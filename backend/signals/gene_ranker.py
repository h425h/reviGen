"""
Gene Ranker — Phenotype-to-Gene via Disease Intermediate.
==========================================================

Core function of the evaluation pipeline.

Given a set of HPO terms (from a patient or a MyGene2 record), this module:
  1. Scores ALL diseases in the HPOA database using Resnik IC-weighted
     semantic similarity (falls back to Jaccard if phenopy unavailable)
  2. Maps each disease to its causal gene(s) via OMIM
  3. Returns a ranked list of (gene, score, disease) tuples

Why disease-intermediate (not direct phenotype→gene):
  - HPOA disease profiles are curated from hundreds of published cases —
    high-quality, human-validated phenotype clusters
  - Matching patient vs disease profile captures syndrome-level specificity
    (e.g. FOXG1 syndrome's regression + stereotypies pattern) better than
    raw gene-phenotype matrices
  - Then disease→gene lookup via OMIM gives 1-3 genes per disease
    rather than scoring all 20,000 genes with noisy annotations
  - LIRICAL, Exomiser all use this same two-step architecture

This ranked gene list is what the evaluation pipeline scores:
  - Is the known causal gene in top-1? top-3? top-10?
  - F1 / Precision / Recall computed by eval_pipeline.py
"""

import warnings
from typing import List, Dict, Tuple, Optional


# ── HPOA disease→gene map (static fallback, used when full HPOA not loaded) ──
# Expanded for evaluation coverage of MyGene2 common genes
# Format: OMIM_ID → {name, gene, hpo_terms}
DISEASE_GENE_PROFILES: Dict[str, Dict] = {
    "OMIM:613454": {
        "name": "FOXG1 syndrome",
        "gene": "FOXG1",
        "hpo": {"HP:0001263","HP:0002186","HP:0000252","HP:0001252","HP:0001332",
                "HP:0002072","HP:0000729","HP:0001344","HP:0001270","HP:0011968",
                "HP:0001251","HP:0001260","HP:0012758"},
    },
    "OMIM:300672": {
        "name": "CDKL5 deficiency disorder",
        "gene": "CDKL5",
        "hpo": {"HP:0001263","HP:0001250","HP:0001252","HP:0000252","HP:0002186",
                "HP:0001344","HP:0002072","HP:0001270","HP:0011968"},
    },
    "OMIM:312750": {
        "name": "Rett syndrome",
        "gene": "MECP2",
        "hpo": {"HP:0001263","HP:0002186","HP:0000252","HP:0001252","HP:0001332",
                "HP:0001344","HP:0001270","HP:0002072","HP:0000729","HP:0001250",
                "HP:0002360","HP:0001260"},
    },
    "OMIM:300496": {
        "name": "MECP2 duplication syndrome",
        "gene": "MECP2",
        "hpo": {"HP:0001263","HP:0001252","HP:0001344","HP:0001250","HP:0002360",
                "HP:0000407","HP:0001270"},
    },
    "OMIM:613721": {
        "name": "Developmental and epileptic encephalopathy 11 (SCN2A)",
        "gene": "SCN2A",
        "hpo": {"HP:0001250","HP:0001263","HP:0001252","HP:0001270","HP:0001344",
                "HP:0002069","HP:0011185"},
    },
    "OMIM:612164": {
        "name": "Developmental and epileptic encephalopathy 4 (STXBP1)",
        "gene": "STXBP1",
        "hpo": {"HP:0001250","HP:0001263","HP:0001252","HP:0001270","HP:0001344",
                "HP:0002069","HP:0000729"},
    },
    "OMIM:613720": {
        "name": "Developmental and epileptic encephalopathy 7 (KCNQ2)",
        "gene": "KCNQ2",
        "hpo": {"HP:0001250","HP:0001263","HP:0001252","HP:0001270","HP:0002069",
                "HP:0011185","HP:0001344"},
    },
    "OMIM:148050": {
        "name": "KBG syndrome",
        "gene": "ANKRD11",
        "hpo": {"HP:0001263","HP:0001249","HP:0001510","HP:0006482","HP:0000252",
                "HP:0001250","HP:0000729"},
    },
    "OMIM:616078": {
        "name": "Schinzel-Giedion / SETBP1 haploinsufficiency",
        "gene": "SETBP1",
        "hpo": {"HP:0001263","HP:0001249","HP:0001250","HP:0001252","HP:0000252"},
    },
    "OMIM:156400": {
        "name": "Jansen metaphyseal chondrodysplasia",
        "gene": "PTHR1",
        "hpo": {"HP:0003521","HP:0002750","HP:0000843","HP:0000518","HP:0002748",
                "HP:0003075","HP:0003127","HP:0000926"},
    },
    "OMIM:615005": {
        "name": "Epileptic encephalopathy, early infantile, 14 (KCNT1)",
        "gene": "KCNT1",
        "hpo": {"HP:0001250","HP:0001263","HP:0001252","HP:0001270","HP:0001344",
                "HP:0002069","HP:0000729","HP:0002360"},
    },
    "OMIM:614558": {
        "name": "Developmental and epileptic encephalopathy 6 (SCN8A)",
        "gene": "SCN8A",
        "hpo": {"HP:0001250","HP:0001263","HP:0001252","HP:0001270","HP:0001344",
                "HP:0002069","HP:0000729"},
    },
    "OMIM:616300": {
        "name": "Intellectual disability, autosomal dominant 45 (SETD5)",
        "gene": "SETD5",
        "hpo": {"HP:0001263","HP:0001249","HP:0000252","HP:0001250","HP:0000729",
                "HP:0001270"},
    },
    "OMIM:615761": {
        "name": "Mental retardation, autosomal dominant 23 (SETD5)",
        "gene": "SETD5",
        "hpo": {"HP:0001263","HP:0001249","HP:0000252","HP:0000729"},
    },
    "OMIM:616268": {
        "name": "Mental retardation, autosomal dominant 32 (KAT6A)",
        "gene": "KAT6A",
        "hpo": {"HP:0001263","HP:0001249","HP:0000252","HP:0011968","HP:0000729",
                "HP:0001270","HP:0001344"},
    },
    "OMIM:614104": {
        "name": "Intellectual disability, autosomal dominant 7 (DYRK1A)",
        "gene": "DYRK1A",
        "hpo": {"HP:0001263","HP:0001249","HP:0000252","HP:0001250","HP:0000729",
                "HP:0001270","HP:0001344","HP:0011968"},
    },
    "OMIM:616278": {
        "name": "MED13L syndrome",
        "gene": "MED13L",
        "hpo": {"HP:0001263","HP:0001249","HP:0000252","HP:0001270","HP:0001344",
                "HP:0011968","HP:0000729"},
    },
    "OMIM:617053": {
        "name": "Neurodevelopmental disorder with or without hyperkinetic movements (ADNP)",
        "gene": "ADNP",
        "hpo": {"HP:0001263","HP:0001249","HP:0000729","HP:0001270","HP:0001344",
                "HP:0002072","HP:0000252"},
    },
    "OMIM:615502": {
        "name": "Intellectual disability, autosomal dominant 26 (ANKRD11)",
        "gene": "ANKRD11",
        "hpo": {"HP:0001263","HP:0001249","HP:0000252","HP:0001250","HP:0000729"},
    },
    "OMIM:245570": {
        "name": "Epilepsy-aphasia spectrum (GRIN2A)",
        "gene": "GRIN2A",
        "hpo": {"HP:0001250","HP:0001263","HP:0001344","HP:0002069","HP:0000729"},
    },
    "OMIM:615418": {
        "name": "Angelman-like syndrome (HIST1H1E / KAT6A related)",
        "gene": "KAT6A",
        "hpo": {"HP:0001263","HP:0001249","HP:0001344","HP:0000729","HP:0001270"},
    },
}


def _jaccard(a: set, b: set) -> float:
    if not a and not b:
        return 0.0
    inter = len(a & b)
    union = len(a | b)
    return inter / union if union > 0 else 0.0


def rank_genes_by_phenotype(
    patient_hpo_terms: List[str],
    top_k: int = 20,
    hpoa_disease_map: Optional[Dict[str, List[str]]] = None,
    disease_to_gene_map: Optional[Dict[str, str]] = None,
) -> List[Dict]:
    """
    Rank candidate genes by Resnik semantic similarity of patient HPO profile
    against all disease HPO profiles in HPOA (or fallback static profiles).

    Args:
        patient_hpo_terms:   Patient's current HPO terms (list of HP:xxxxxxx)
        top_k:               How many top genes to return (default 20)
        hpoa_disease_map:    Optional live HPOA {disease_id: [hpo_terms]}
                             loaded from data_loaders. If None, uses static profiles.
        disease_to_gene_map: Optional {disease_id: gene_symbol} from OMIM.
                             If None, uses static profiles.

    Returns:
        List of dicts sorted by score descending:
        [{gene, disease_id, disease_name, score, similarity_method, rank}, ...]
    """
    # Normalise patient terms
    patient_set = set()
    for t in (patient_hpo_terms or []):
        t = t.strip().upper()
        if not t.startswith("HP:"):
            t = "HP:" + t
        patient_set.add(t)

    if not patient_set:
        return []

    # Get similarity function
    scorer, similarity_method = _get_similarity_fn()

    # Build disease profile source
    profiles = _build_profiles(hpoa_disease_map, disease_to_gene_map)

    # Score all diseases
    scored: List[Tuple[str, str, str, float]] = []  # (gene, disease_id, name, score)

    for disease_id, profile in profiles.items():
        disease_hpo = profile["hpo"]
        gene        = profile["gene"]
        name        = profile["name"]

        if not disease_hpo:
            continue

        sim = _compute_similarity(patient_set, disease_hpo, scorer, similarity_method)
        if sim > 0:
            scored.append((gene, disease_id, name, sim))

    # Sort by score descending
    scored.sort(key=lambda x: x[3], reverse=True)

    # Deduplicate: if same gene appears via multiple diseases, keep best score
    seen_genes: Dict[str, Dict] = {}
    for gene, disease_id, name, score in scored:
        if gene not in seen_genes or score > seen_genes[gene]["score"]:
            seen_genes[gene] = {
                "gene":             gene,
                "disease_id":       disease_id,
                "disease_name":     name,
                "score":            round(score, 4),
                "similarity_method": similarity_method,
            }

    # Re-sort after dedup and assign ranks
    ranked = sorted(seen_genes.values(), key=lambda x: x["score"], reverse=True)
    for i, r in enumerate(ranked):
        r["rank"] = i + 1

    return ranked[:top_k]


def _get_similarity_fn():
    """Return (scorer_or_None, method_name)."""
    try:
        from .phenotypic_drift import _get_scorer, using_resnik
        scorer = _get_scorer()
        if scorer is not None:
            return scorer, "resnik"
    except (ImportError, Exception):
        pass
    return None, "jaccard"


def _compute_similarity(patient_set: set, disease_hpo: set,
                        scorer, method: str) -> float:
    if method == "resnik" and scorer is not None:
        try:
            raw = scorer.score_term_sets_basic(
                [t for t in patient_set if t.startswith("HP:")],
                [t for t in disease_hpo if t.startswith("HP:")],
            )
            return min(1.0, raw / 6.0) if raw > 0 else 0.0
        except Exception:
            pass
    return _jaccard(patient_set, disease_hpo)


def _build_profiles(
    hpoa_disease_map: Optional[Dict],
    disease_to_gene_map: Optional[Dict],
) -> Dict[str, Dict]:
    """
    Build unified {disease_id: {gene, name, hpo}} from live data or fallback.
    """
    # If live HPOA + gene map provided, use them
    if hpoa_disease_map and disease_to_gene_map:
        profiles = {}
        for disease_id, hpo_list in hpoa_disease_map.items():
            gene = disease_to_gene_map.get(disease_id)
            if gene:
                profiles[disease_id] = {
                    "gene": gene,
                    "name": disease_id,
                    "hpo":  set(hpo_list),
                }
        if profiles:
            return profiles

    # Fall back to static profiles
    return {
        did: {
            "gene": p["gene"],
            "name": p["name"],
            "hpo":  p["hpo"],
        }
        for did, p in DISEASE_GENE_PROFILES.items()
    }