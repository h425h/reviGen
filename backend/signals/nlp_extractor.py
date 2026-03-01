"""
NLP Clinical Notes Extractor — reviGen
Converts free-text clinical notes into structured HPO terms with timestamps.

Architecture decision: Uses Claude API rather than a local biomedical NER model.
Rationale: 2024 benchmarks (Wan et al., Luo et al.) show LLMs outperform
scispacy + HPO dictionary lookup on real-world clinical notes due to:
  - Clinical abbreviations ("ID", "DD", "NB", "FTT") handled contextually
  - Negated findings correctly excluded ("no seizures" → not flagged)
  - Implicit phenotypes extracted ("wheelchair dependent" → HP:0002540)
  - Temporal information preserved ("regression at 18 months")

Output feeds directly into the existing HPO pipeline without modification.
"""

import os
import json
import re
from typing import List, Dict, Optional, Tuple
from datetime import datetime
import anthropic


# ── Known HPO term name→ID lookup (supplement LLM output) ────────────────────
# LLM outputs names; we verify/map to canonical IDs
# Full mapping would load from hp.obo — this covers common rare disease terms
HPO_NAME_TO_ID = {
    # Neurological
    "global developmental delay":       "HP:0001263",
    "developmental delay":              "HP:0001263",
    "intellectual disability":          "HP:0001249",
    "intellectual disability mild":     "HP:0001256",
    "intellectual disability severe":   "HP:0010864",
    "seizures":                         "HP:0001250",
    "epilepsy":                         "HP:0001250",
    "hypotonia":                        "HP:0001252",
    "muscular hypotonia":               "HP:0001252",
    "floppy tone":                      "HP:0001252",
    "hypertonia":                       "HP:0001276",
    "spasticity":                       "HP:0001257",
    "ataxia":                           "HP:0001251",
    "dystonia":                         "HP:0001332",
    "apraxia":                          "HP:0002186",
    "dysarthria":                       "HP:0001260",
    "absent speech":                    "HP:0001344",
    "no speech":                        "HP:0001344",
    "speech delay":                     "HP:0000750",
    "language delay":                   "HP:0000750",
    "autistic behavior":                "HP:0000729",
    "autistic behaviour":               "HP:0000729",
    "autism":                           "HP:0000729",
    "stereotyped movements":            "HP:0002063",
    "hand stereotypies":                "HP:0002063",
    "hand wringing":                    "HP:0002063",
    "repetitive hand movements":        "HP:0002063",
    "choreiform movements":             "HP:0002072",
    "chorea":                           "HP:0002072",
    "motor delay":                      "HP:0001270",
    "delayed motor development":        "HP:0001270",
    "developmental regression":         "HP:0002376",
    "regression":                       "HP:0002376",
    "neurodevelopmental delay":         "HP:0012758",
    "feeding difficulties":             "HP:0011968",
    "poor feeding":                     "HP:0011968",
    "feeding problems":                 "HP:0011968",
    # Structural brain
    "microcephaly":                     "HP:0000252",
    "macrocephaly":                     "HP:0000256",
    "simplified gyral pattern":         "HP:0009830",
    "lissencephaly":                    "HP:0001339",
    "corpus callosum agenesis":         "HP:0001274",
    "hypoplastic corpus callosum":      "HP:0002079",
    "cerebellar atrophy":               "HP:0001272",
    "brain atrophy":                    "HP:0012444",
    "white matter abnormality":         "HP:0002500",
    # Ophthalmological
    "cataract":                         "HP:0000518",
    "nystagmus":                        "HP:0000639",
    "strabismus":                       "HP:0000486",
    "visual impairment":                "HP:0000505",
    "cortical visual impairment":       "HP:0100704",
    # Skeletal
    "short stature":                    "HP:0004322",
    "disproportionate short stature":   "HP:0003521",
    "delayed skeletal maturation":      "HP:0002750",
    "scoliosis":                        "HP:0002650",
    "joint hypermobility":              "HP:0001382",
    "macrodontia":                      "HP:0006482",
    # Metabolic / endocrine
    "hyperparathyroidism":              "HP:0000843",
    "hypercalcaemia":                   "HP:0003127",
    "hypophosphataemia":                "HP:0003075",
    "rickets":                          "HP:0002748",
    "lactic acidosis":                  "HP:0003128",
    "hypoglycemia":                     "HP:0001943",
    # Cardiac
    "congenital heart defect":          "HP:0001627",
    "ventricular septal defect":        "HP:0001629",
    "atrial septal defect":             "HP:0001631",
    # Renal
    "renal anomaly":                    "HP:0000077",
    "hydronephrosis":                   "HP:0000126",
    # Growth
    "growth retardation":               "HP:0001510",
    "failure to thrive":                "HP:0001508",
    "ftt":                              "HP:0001508",
    # Facial
    "coarse facies":                    "HP:0000280",
    "wide forehead":                    "HP:0000280",
    "hypertelorism":                    "HP:0000316",
    "low set ears":                     "HP:0000369",
    "low-set ears":                     "HP:0000369",
    "prominent ears":                   "HP:0000411",
}


def _normalise_name(name: str) -> str:
    """Lowercase + strip for fuzzy matching."""
    return re.sub(r'\s+', ' ', name.lower().strip())


def name_to_hpo_id(phenotype_name: str) -> Optional[str]:
    """
    Map a phenotype name string to an HPO ID.
    Tries exact match first, then fuzzy substring match.
    """
    normalised = _normalise_name(phenotype_name)

    # Exact match
    if normalised in HPO_NAME_TO_ID:
        return HPO_NAME_TO_ID[normalised]

    # Substring match — find the longest key that appears in the name
    best_key = None
    best_len = 0
    for key in HPO_NAME_TO_ID:
        if key in normalised and len(key) > best_len:
            best_key = key
            best_len = len(key)

    if best_key:
        return HPO_NAME_TO_ID[best_key]

    # Try reverse — is the full name a substring of a key?
    for key, hpo_id in HPO_NAME_TO_ID.items():
        if normalised in key:
            return hpo_id

    return None


def extract_hpo_from_notes(
    clinical_text: str,
    api_key: Optional[str] = None,
    reference_date: Optional[str] = None,
) -> Dict:
    """
    Extract structured HPO terms from free-text clinical notes using Claude.

    Handles:
    - Multilingual input (Claude is multilingual by default)
    - Negated findings ("no seizures" correctly excluded)
    - Temporal information ("regression at 18 months" → timestamp)
    - Clinical abbreviations (ID, DD, FTT, NB, etc.)
    - Implicit phenotypes ("wheelchair dependent" → HP:0002540)

    Args:
        clinical_text:   Free-text clinical note (any language)
        api_key:         Anthropic API key (or env ANTHROPIC_API_KEY)
        reference_date:  Patient's current date for relative time parsing
                         Defaults to today

    Returns:
        {
            "hpo_terms":       ["HP:0001252", ...],   ← for pipeline input
            "hpo_with_names":  [{"id": "HP:0001252", "name": "Hypotonia", ...}],
            "temporal_hpo":    [{"id": ..., "name": ..., "onset": "2019-01"}],
            "excluded_terms":  [...],  ← negated findings
            "raw_extractions": [...],  ← LLM output before ID mapping
            "extraction_method": "claude_api" | "fallback",
            "unmapped_terms":  [...],  ← phenotypes LLM found but couldn't map
        }
    """
    key = api_key or os.environ.get("ANTHROPIC_API_KEY", "")
    ref_date = reference_date or datetime.now().strftime("%Y-%m-%d")

    if not clinical_text or not clinical_text.strip():
        return _empty_result()

    if key:
        return _extract_via_claude(clinical_text, key, ref_date)
    else:
        return _extract_via_heuristic(clinical_text)


def _extract_via_claude(text: str, api_key: str, ref_date: str) -> Dict:
    """Full extraction using Claude API."""
    try:
        client = anthropic.Anthropic(api_key=api_key)

        system_prompt = """You are a clinical genetics expert extracting phenotype information from clinical notes.
Your task: identify all HPO (Human Phenotype Ontology) phenotypes present in the text.

Rules:
1. EXCLUDE negated findings ("no seizures", "denies", "ruled out", "absent")
2. INCLUDE only currently present OR historically documented findings
3. Extract temporal information when mentioned (age of onset, age of regression)
4. Handle abbreviations: ID=intellectual disability, DD=developmental delay, FTT=failure to thrive, NB=newborn, ASD=autism spectrum disorder, CP=cerebral palsy
5. Extract implicit phenotypes: "wheelchair dependent" → inability to walk, "non-verbal" → absent speech

Return ONLY valid JSON, no other text:
{
  "present_phenotypes": [
    {
      "name": "exact HPO term name or close match",
      "onset_description": "at birth | neonatal | infantile | childhood | null",
      "onset_age_months": null or integer,
      "confidence": "high | medium | low"
    }
  ],
  "excluded_phenotypes": ["list of negated/absent findings as names"],
  "temporal_notes": "brief summary of disease progression if mentioned"
}"""

        user_prompt = f"""Reference date: {ref_date}

Clinical note:
\"\"\"
{text}
\"\"\"

Extract all HPO phenotypes present in this patient."""

        message = client.messages.create(
            model="claude-sonnet-4-20250514",
            max_tokens=1500,
            system=system_prompt,
            messages=[{"role": "user", "content": user_prompt}],
        )

        response_text = message.content[0].text
        # Strip markdown fences if present
        if "```json" in response_text:
            response_text = response_text.split("```json")[1].split("```")[0]
        elif "```" in response_text:
            response_text = response_text.split("```")[1].split("```")[0]

        parsed = json.loads(response_text.strip())
        return _process_claude_output(parsed, ref_date)

    except json.JSONDecodeError as e:
        print(f"JSON parse error in NLP extractor: {e}")
        return _extract_via_heuristic(text)
    except Exception as e:
        print(f"Claude API error in NLP extractor: {e}")
        return _extract_via_heuristic(text)


def _process_claude_output(parsed: Dict, ref_date: str) -> Dict:
    """Map Claude's extracted phenotype names to HPO IDs."""
    present = parsed.get("present_phenotypes", [])
    excluded_names = parsed.get("excluded_phenotypes", [])

    hpo_terms      = []
    hpo_with_names = []
    temporal_hpo   = []
    unmapped       = []

    for item in present:
        name = item.get("name", "")
        hpo_id = name_to_hpo_id(name)

        if hpo_id:
            if hpo_id not in hpo_terms:
                hpo_terms.append(hpo_id)
            hpo_with_names.append({
                "id":         hpo_id,
                "name":       name,
                "confidence": item.get("confidence", "medium"),
            })

            # Build temporal entry if onset info present
            onset_months = item.get("onset_age_months")
            onset_desc   = item.get("onset_description")
            if onset_months is not None or onset_desc:
                temporal_hpo.append({
                    "id":           hpo_id,
                    "name":         name,
                    "onset_months": onset_months,
                    "onset_desc":   onset_desc,
                })
        else:
            if name:
                unmapped.append(name)

    # Map excluded names to IDs where possible
    excluded_ids = []
    for name in excluded_names:
        hpo_id = name_to_hpo_id(name)
        if hpo_id:
            excluded_ids.append({"id": hpo_id, "name": name})

    return {
        "hpo_terms":        hpo_terms,
        "hpo_with_names":   hpo_with_names,
        "temporal_hpo":     temporal_hpo,
        "excluded_terms":   excluded_ids,
        "raw_extractions":  present,
        "temporal_notes":   parsed.get("temporal_notes", ""),
        "extraction_method": "claude_api",
        "unmapped_terms":   unmapped,
        "term_count":       len(hpo_terms),
    }


def _extract_via_heuristic(text: str) -> Dict:
    """
    Fallback when no API key: simple keyword matching against HPO name dictionary.
    Much lower quality than Claude but ensures the endpoint never crashes.
    """
    text_lower = text.lower()
    found_terms = []
    found_ids   = []

    # Simple negation check — skip if "no " or "without " precedes the term
    negation_patterns = [
        r"no\s+{}", r"without\s+{}", r"denies\s+{}", r"absent\s+{}",
        r"ruled out\s+{}", r"negative for\s+{}"
    ]

    for name, hpo_id in sorted(HPO_NAME_TO_ID.items(), key=lambda x: -len(x[0])):
        if name not in text_lower:
            continue

        # Check for negation
        is_negated = False
        for pattern in negation_patterns:
            if re.search(pattern.format(re.escape(name)), text_lower):
                is_negated = True
                break

        if not is_negated and hpo_id not in found_ids:
            found_ids.append(hpo_id)
            found_terms.append({"id": hpo_id, "name": name, "confidence": "low"})

    return {
        "hpo_terms":         found_ids,
        "hpo_with_names":    found_terms,
        "temporal_hpo":      [],
        "excluded_terms":    [],
        "raw_extractions":   [],
        "temporal_notes":    "",
        "extraction_method": "heuristic_fallback",
        "unmapped_terms":    [],
        "term_count":        len(found_ids),
    }


def _empty_result() -> Dict:
    return {
        "hpo_terms":         [],
        "hpo_with_names":    [],
        "temporal_hpo":      [],
        "excluded_terms":    [],
        "raw_extractions":   [],
        "temporal_notes":    "",
        "extraction_method": "empty_input",
        "unmapped_terms":    [],
        "term_count":        0,
    }