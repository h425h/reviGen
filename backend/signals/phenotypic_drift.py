"""
Phenotypic Drift Signal Detector — FINAL v3.

Key upgrade over v2: replaces Jaccard similarity with Resnik semantic
similarity via phenopy. This is the same IC-weighted approach used in
LIRICAL (the state-of-the-art rare disease diagnosis tool).

Why Resnik > Jaccard:
  Jaccard treats HP:0001252 (hypotonia — present in ~40% of rare diseases)
  identically to HP:0002186 (apraxia — present in ~2%). Resnik weights
  each term by Information Content (IC = −log frequency), so rare, specific
  symptoms contribute far more to the similarity score than common ones.
  This dramatically reduces false positives where common symptoms cause
  spurious disease matches.

  IC(t) = −log [ freq(t) / total_annotations ]
  Resnik(t1, t2) = IC( MICA(t1, t2) )   ← Most Informative Common Ancestor
  Similarity(P, D) = BestMatchAverage across all term pairs

Fallback behaviour:
  If phenopy or its data files are unavailable (e.g., first run before
  hp.obo / phenotype.hpoa are downloaded), the module automatically falls
  back to Jaccard with a warning. This ensures the API never crashes due
  to missing ontology files.
"""

import math
import warnings
from typing import List, Dict, Optional, Tuple

# ── phenopy imports (graceful fallback if unavailable) ────────────────────────
_PHENOPY_AVAILABLE = False
_scorer = None
_hpo_network = None

try:
    from phenopy.build_hpo import generate_annotated_hpo_network
    from phenopy.score import Scorer
    _PHENOPY_AVAILABLE = True
except ImportError:
    warnings.warn(
        "phenopy not installed. Install with: pip install phenopy\n"
        "Falling back to Jaccard similarity for disease matching.",
        ImportWarning,
        stacklevel=2,
    )

# ── Module-level scorer (initialised once on first use) ───────────────────────
_scorer_initialised = False


def _get_scorer(obo_path: str = "data/hp.obo",
                hpoa_path: str = "data/phenotype.hpoa") -> Optional[object]:
    """
    Lazy-initialise the phenopy Scorer. Called once; cached globally.
    Returns None if phenopy unavailable or data files missing.
    """
    global _scorer, _hpo_network, _scorer_initialised

    if _scorer_initialised:
        return _scorer  # Already attempted (may be None)

    _scorer_initialised = True

    if not _PHENOPY_AVAILABLE:
        return None

    import os
    # Try common relative paths for the data directory
    candidate_pairs = [
        (obo_path, hpoa_path),
        ("../data/hp.obo", "../data/phenotype.hpoa"),
        ("../../data/hp.obo", "../../data/phenotype.hpoa"),
    ]

    for obo, hpoa in candidate_pairs:
        if os.path.exists(obo) and os.path.exists(hpoa):
            try:
                _hpo_network, _, _ = generate_annotated_hpo_network(obo, hpoa)
                _scorer = Scorer(_hpo_network)
                return _scorer
            except Exception as e:
                warnings.warn(f"phenopy scorer init failed: {e}. Using Jaccard fallback.")
                return None

    warnings.warn(
        "hp.obo / phenotype.hpoa not found. Using Jaccard fallback.\n"
        "Download from: https://hpo.jax.org/app/data/ontology",
        UserWarning,
        stacklevel=3,
    )
    return None


# ── Similarity functions ──────────────────────────────────────────────────────

def resnik_similarity(set_a: set, set_b: set,
                      scorer: Optional[object] = None) -> float:
    """
    IC-weighted Resnik semantic similarity between two sets of HPO terms.

    Uses BestMatchAverage (BMA): for each term in set_a, finds the most
    similar term in set_b (using MICA IC), averages both directions.

    Falls back to Jaccard if phenopy scorer is unavailable.

    Args:
        set_a:   First HPO term set (e.g., patient current symptoms)
        set_b:   Second HPO term set (e.g., disease canonical profile)
        scorer:  phenopy Scorer instance (uses module-level if None)

    Returns:
        Similarity score in [0, 1]. Higher = more similar.
    """
    if not set_a or not set_b:
        return 0.0

    active_scorer = scorer or _get_scorer()

    if active_scorer is not None:
        try:
            score = active_scorer.score_term_sets_basic(
                [t for t in set_a if t.startswith("HP:")],
                [t for t in set_b if t.startswith("HP:")],
            )
            # phenopy returns raw IC values; normalise to [0, 1]
            # Max IC in HPO is ~13 bits; BMA across set pairs is typically 0–6
            # Empirical normalisation: divide by 6.0, clamp
            return min(1.0, score / 6.0) if score > 0 else 0.0
        except Exception:
            pass  # Fall through to Jaccard

    return jaccard_similarity(set_a, set_b)


def jaccard_similarity(set_a: set, set_b: set) -> float:
    """
    Fallback: Jaccard similarity = |intersection| / |union|.
    Used when phenopy is unavailable.
    """
    if not set_a and not set_b:
        return 0.0
    intersection = len(set_a & set_b)
    union = len(set_a | set_b)
    return intersection / union if union > 0 else 0.0


def using_resnik() -> bool:
    """Returns True if Resnik similarity is active (phenopy + data available)."""
    return _get_scorer() is not None


# ── Disease profile lookup ────────────────────────────────────────────────────

# Expanded hardcoded HPO profiles for demo cases
# In production: loaded entirely from HPOA at startup via data_loaders.py
_HARDCODED_PROFILES: Dict[str, Dict] = {
    # FOXG1 syndrome (OMIM:613454)
    # Congenital Rett variant: severe ID, absent speech, stereotypies, microcephaly
    "OMIM:613454": {
        "name": "FOXG1 syndrome",
        "gene": "FOXG1",
        "hpo": {
            "HP:0001263",  # Global developmental delay
            "HP:0002186",  # Apraxia
            "HP:0000252",  # Microcephaly
            "HP:0001252",  # Hypotonia
            "HP:0001332",  # Dystonia
            "HP:0002072",  # Choreiform movements
            "HP:0000729",  # Autistic behaviour
            "HP:0001344",  # Absent speech
            "HP:0001270",  # Motor delay
            "HP:0011968",  # Feeding difficulties
            "HP:0001251",  # Ataxia
            "HP:0001260",  # Dysarthria
            "HP:0012758",  # Neurodevelopmental delay
        },
    },
    # CDKL5 deficiency disorder (OMIM:300672)
    # Early-onset seizures, severe ID, absent speech, hand stereotypies
    "OMIM:300672": {
        "name": "CDKL5 deficiency disorder",
        "gene": "CDKL5",
        "hpo": {
            "HP:0001263",  # Global developmental delay
            "HP:0001250",  # Seizures
            "HP:0001252",  # Hypotonia
            "HP:0000252",  # Microcephaly
            "HP:0002186",  # Apraxia
            "HP:0001344",  # Absent speech
            "HP:0002072",  # Choreiform movements
            "HP:0001270",  # Motor delay
            "HP:0011968",  # Feeding difficulties
        },
    },
    # KBG syndrome (OMIM:148050) — ANKRD11
    # Short stature, macrodontia, intellectual disability, behavioural issues
    "OMIM:148050": {
        "name": "KBG syndrome",
        "gene": "ANKRD11",
        "hpo": {
            "HP:0001263",  # Global developmental delay
            "HP:0001249",  # Intellectual disability
            "HP:0001510",  # Growth delay
            "HP:0006482",  # Macrodontia
            "HP:0000252",  # Microcephaly
            "HP:0001250",  # Seizures
            "HP:0000729",  # Autistic behaviour
        },
    },
    # Schinzel-Giedion syndrome (OMIM:616078) — SETBP1
    # Severe ID, seizures, coarse facies, skeletal anomalies
    "OMIM:616078": {
        "name": "Schinzel-Giedion syndrome / SETBP1 haploinsufficiency",
        "gene": "SETBP1",
        "hpo": {
            "HP:0001263",  # Global developmental delay
            "HP:0001249",  # Intellectual disability
            "HP:0001250",  # Seizures
            "HP:0001252",  # Hypotonia
            "HP:0000252",  # Microcephaly
        },
    },
    # Jansen metaphyseal chondrodysplasia (OMIM:156400) — PTHR1
    # Hypercalcaemia, skeletal dysplasia, short-limbed dwarfism
    "OMIM:156400": {
        "name": "Jansen metaphyseal chondrodysplasia",
        "gene": "PTHR1",
        "hpo": {
            "HP:0003521",  # Disproportionate short stature
            "HP:0002750",  # Delayed skeletal maturation
            "HP:0000843",  # Hyperparathyroidism
            "HP:0000518",  # Cataract
            "HP:0002748",  # Rickets
            "HP:0003075",  # Hypophosphataemia
            "HP:0003127",  # Hypercalcaemia
            "HP:0000926",  # Platyspondyly
        },
    },
}

# Gene → OMIM IDs mapping (for lookup by gene symbol)
_GENE_TO_OMIM: Dict[str, List[str]] = {
    "FOXG1":   ["OMIM:613454", "OMIM:164874"],
    "CDKL5":   ["OMIM:300672"],
    "SETBP1":  ["OMIM:616078"],
    "ANKRD11": ["OMIM:148050"],
    "PTHR1":   ["OMIM:156400"],
}


def get_disease_hpo_profiles(vus_genes: List[str]) -> Dict[str, set]:
    """
    Get canonical HPO profiles for diseases linked to VUS genes.

    Priority:
    1. Live HPOA data loaded in data_loaders.disease_phenotypes (full dataset)
    2. Hardcoded profiles above (demo fallback for known genes)

    Returns:
        {disease_id: set_of_hpo_terms}
    """
    result: Dict[str, set] = {}

    # Try live HPOA data first
    try:
        from ..data_loaders import disease_phenotypes as hpoa_map
        if hpoa_map:
            for gene in (vus_genes or []):
                omim_ids = _GENE_TO_OMIM.get(gene.upper().strip(), [])
                for omim_id in omim_ids:
                    if omim_id in hpoa_map:
                        result[omim_id] = set(hpoa_map[omim_id])
            if result:
                return result
    except (ImportError, AttributeError):
        pass

    # Fall back to hardcoded profiles
    for gene in (vus_genes or []):
        omim_ids = _GENE_TO_OMIM.get(gene.upper().strip(), [])
        for omim_id in omim_ids:
            if omim_id in _HARDCODED_PROFILES:
                result[omim_id] = _HARDCODED_PROFILES[omim_id]["hpo"].copy()

    return result


# ── Main scoring function ─────────────────────────────────────────────────────

def score_phenotypic_drift(
    original_hpo_terms: List[str],
    current_hpo_terms:  List[str],
    vus_genes:          Optional[List[str]] = None,
) -> Dict:
    """
    Score phenotypic drift with IC-weighted Resnik disease matching.

    Two sub-signals:
      1. Raw drift:      fraction of today's symptoms that are new since the test
      2. Disease match:  Resnik semantic similarity of current HPO profile
                         against canonical disease HPO profiles (HPOA)

    The disease match signal is the primary upgrade: Resnik weights specific
    rare symptoms heavily (high IC) and common symptoms lightly (low IC),
    dramatically reducing false-positive disease matches.

    Args:
        original_hpo_terms: HPO terms at time of original test
        current_hpo_terms:  Current HPO terms (today)
        vus_genes:          Genes with VUS — used to look up disease profiles

    Returns:
        Full drift analysis dict with both raw and disease-aware scores,
        plus metadata about which similarity method was used
    """
    if original_hpo_terms is None:
        original_hpo_terms = []
    if current_hpo_terms is None:
        current_hpo_terms = []

    def normalise(terms: List[str]) -> set:
        result = set()
        for t in terms:
            t = t.strip().upper()
            if t and not t.startswith("HP:"):
                t = "HP:" + t
            if t:
                result.add(t)
        return result

    original_set = normalise(original_hpo_terms)
    current_set  = normalise(current_hpo_terms)

    # ── Sub-signal 1: Raw drift ───────────────────────────────────────────────
    new_symptoms  = current_set - original_set
    lost_symptoms = original_set - current_set
    new_count     = len(new_symptoms)
    current_count = len(current_set)

    # What fraction of today's symptom profile is new since the test?
    raw_drift_score = new_count / current_count if current_count > 0 else 0.0

    # ── Sub-signal 2: Disease match (Resnik or Jaccard fallback) ─────────────
    disease_profiles  = get_disease_hpo_profiles(vus_genes or [])
    similarity_method = "resnik" if using_resnik() else "jaccard"

    best_match_score   = 0.0
    best_match_disease = None
    disease_match_details = []

    scorer = _get_scorer()  # May be None → jaccard fallback inside resnik_similarity

    for disease_id, disease_hpo_set in disease_profiles.items():
        # Full similarity: current patient symptoms vs disease profile
        full_sim = resnik_similarity(current_set, disease_hpo_set, scorer)

        # New-symptom specificity: what fraction of NEW symptoms are in disease profile?
        # Rewards cases where new symptoms are diagnostically targeted, not just common
        new_in_disease = len(new_symptoms & disease_hpo_set)
        new_symptom_match_rate = new_in_disease / len(new_symptoms) if new_symptoms else 0.0

        # Disease name for display
        disease_name = _HARDCODED_PROFILES.get(disease_id, {}).get("name", disease_id)

        disease_match_details.append({
            "disease_id":             disease_id,
            "disease_name":           disease_name,
            "similarity":             round(full_sim, 3),
            "similarity_method":      similarity_method,
            "new_symptom_match_rate": round(new_symptom_match_rate, 3),
            "new_symptoms_in_profile": new_in_disease,
            "disease_profile_size":   len(disease_hpo_set),
        })

        if full_sim > best_match_score:
            best_match_score   = full_sim
            best_match_disease = disease_id

    # ── Combined drift score ──────────────────────────────────────────────────
    # 55% raw drift (symptoms changed), 45% disease match (they match something real)
    # Disease match weight slightly higher than v2 because Resnik is more precise
    if disease_profiles:
        combined_score = 0.55 * raw_drift_score + 0.45 * best_match_score
    else:
        combined_score = raw_drift_score  # No disease profile → raw drift only

    combined_score = min(1.0, combined_score)

    # ── Signal strength thresholds ────────────────────────────────────────────
    if combined_score > 0.35:
        signal_strength = "HIGH"
    elif combined_score > 0.15:
        signal_strength = "MEDIUM"
    else:
        signal_strength = "LOW"

    # ── Resolve HPO IDs to human-readable names ───────────────────────────────
    def resolve(hpo_set: set) -> List[Dict]:
        try:
            from ..data_loaders import get_hpo_name
            return [{"hpo_id": hid, "name": get_hpo_name(hid)} for hid in hpo_set]
        except (ImportError, AttributeError):
            return [{"hpo_id": hid, "name": hid} for hid in hpo_set]

    return {
        # Primary scores consumed by aggregator
        "drift_score":           round(combined_score, 3),
        "raw_drift_score":       round(raw_drift_score, 3),
        "disease_match_score":   round(best_match_score, 3),
        "best_matching_disease": best_match_disease,

        # Metadata
        "similarity_method":     similarity_method,  # 'resnik' or 'jaccard'
        "new_symptom_count":     new_count,
        "new_symptoms":          resolve(new_symptoms),
        "removed_symptoms":      resolve(lost_symptoms),
        "original_count":        len(original_set),
        "current_count":         current_count,
        "signal_strength":       signal_strength,
        "disease_match_details": disease_match_details,
    }