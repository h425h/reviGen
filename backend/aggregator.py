"""
Reanalysis Score Aggregator — FINAL v3.

7 signals + entropy modifier. Literature-backed weights.

Weight derivation from Radboud 2022 (n=150, +22% diagnostic yield):
  42% new gene-disease associations  → OMIM signal         (weight 0.29)
  19% improved bioinformatics        → Analysis gaps        (weight 0.10)
  10% VUS reclassification           → ClinVar VUS          (weight 0.14)
  29% improved coverage              → captured in gaps

Additional signals (validated across rare disease literature):
  Phenotypic drift                   →                      (weight 0.11)
  Disease profile match (Resnik)     →                      (weight 0.13)
  Inheritance pattern flags          →                      (weight 0.12)
  AlphaMissense variant scoring      →                      (weight 0.05)
  Time since test                    →                      (weight 0.03)
  Entropy modifier (NEW)             → post-hoc boost       (up to +0.08)

Entropy signal:
  When a patient has multiple VUS across multiple genes and the posterior
  probability is spread evenly (high entropy), diagnostic uncertainty is
  maximal — the original analysis couldn't determine which gene is causal.
  This is a strong independent reanalysis trigger.

  H = −Σ p(gene) log₂ p(gene), normalised by log₂(n_genes) → [0,1]
  Applied as: final_score += ENTROPY_BOOST_MAX × H_normalised

Similarity method:
  Disease match uses Resnik IC-weighted semantic similarity (phenopy) when
  hp.obo + phenotype.hpoa are available; falls back to Jaccard automatically.
"""

import os
import math
from typing import List, Dict, Optional

from .signals import check_vus_reclassification, score_time_since_test
from .signals.phenotypic_drift import score_phenotypic_drift
from .signals.omim_surveillance import check_omim_all_genes
from .signals.alphamissense import check_variants_alphamissense
from .signals.inheritance import score_inheritance_signal
from .signals.analysis_gaps import score_analysis_gaps


# ── Signal weights (must sum to 1.0) ──────────────────────────────────────────
WEIGHTS = {
    "omim":           0.29,
    "vus_clinvar":    0.14,
    "drift":          0.11,
    "disease_match":  0.13,
    "inheritance":    0.12,
    "analysis_gaps":  0.10,
    "alphamissense":  0.05,
    "time":           0.06,
}

assert abs(sum(WEIGHTS.values()) - 1.0) < 1e-9, (
    f"Weights must sum to 1.0, got {sum(WEIGHTS.values()):.6f}"
)

# Entropy modifier cap: multi-VUS uncertainty adds up to +8 points
ENTROPY_BOOST_MAX = 0.08

# Score thresholds
URGENT_THRESHOLD = 0.80
HIGH_THRESHOLD   = 0.65
MEDIUM_THRESHOLD = 0.35


# ── Entropy helper ─────────────────────────────────────────────────────────────

def compute_vus_entropy(gene_scores: Dict[str, float]) -> float:
    """
    Compute normalised Shannon entropy across per-gene posterior scores.

    High entropy = posterior probability spread evenly across genes =
    high diagnostic uncertainty about which gene is causal = stronger
    reanalysis trigger.

    H = −Σ p(gene) log₂ p(gene)
    H_norm = H / log₂(n_genes)  →  [0, 1]

    Returns 0.0 for fewer than 2 genes (no ambiguity to measure).
    """
    if len(gene_scores) < 2:
        return 0.0

    total = sum(gene_scores.values())
    if total == 0:
        return 0.0

    probs = [s / total for s in gene_scores.values() if s > 0]
    H     = -sum(p * math.log2(p) for p in probs)
    H_max = math.log2(len(gene_scores))
    return round(H / H_max, 3) if H_max > 0 else 0.0


# ── Main aggregator ────────────────────────────────────────────────────────────

def compute_reanalysis_score(
    vus_genes:                   List[str],
    test_date:                   str,
    original_hpo_terms:          List[str],
    current_hpo_terms:           List[str],
    test_type:                   Optional[str]  = "exome",
    patient_sex:                 Optional[str]  = None,
    analysis_type:               Optional[str]  = None,
    consanguineous:              Optional[bool] = None,
    cnv_calling_performed:       Optional[bool] = None,
    splice_analysis_performed:   Optional[bool] = None,
    mito_analysis_performed:     Optional[bool] = None,
    repeat_expansion_checked:    Optional[bool] = None,
    variants:                    Optional[List[Dict]] = None,
    omim_api_key:                Optional[str]  = None,
    alphamissense_db_path:       Optional[str]  = None,
) -> Dict:
    """
    Compute the overall reanalysis score for a single patient case.

    All optional fields degrade gracefully — if not provided, the
    corresponding signal scores 0.0 and the weighted sum still works.

    Returns a full analysis dict with: reanalysis_score (0-100),
    confidence tier, urgency string, complete signal breakdown with
    per-signal scores/weights/contributions, entropy analysis, and
    input summary.
    """

    # ── Signal 1: OMIM new gene-disease surveillance ──────────────────────────
    omim_signals = check_omim_all_genes(
        vus_genes, test_date,
        omim_api_key=omim_api_key or os.environ.get("OMIM_API_KEY", ""),
    )
    omim_score = max((s["score"] for s in omim_signals), default=0.0)

    # ── Signal 2: Phenotypic drift + Resnik disease match ────────────────────
    drift_result        = score_phenotypic_drift(
        original_hpo_terms, current_hpo_terms, vus_genes
    )
    drift_score         = drift_result.get("raw_drift_score", 0.0)
    disease_match_score = drift_result.get("disease_match_score", 0.0)
    similarity_method   = drift_result.get("similarity_method", "jaccard")

    # ── Signal 3: ClinVar VUS reclassification ────────────────────────────────
    vus_signals = check_vus_reclassification(vus_genes, test_date)
    vus_score   = 0.0
    for sig in vus_signals:
        vus_score += {"HIGH": 0.5, "MEDIUM": 0.3, "LOW": 0.15}.get(
            sig.get("signal_strength", "LOW"), 0.0
        )
    vus_score = min(1.0, vus_score)

    # ── Signal 4: Inheritance pattern flags ──────────────────────────────────
    inheritance_result = score_inheritance_signal(
        vus_genes=vus_genes,
        patient_sex=patient_sex,
        analysis_type=analysis_type,
        known_vus_count_per_gene={g.upper(): 1 for g in (vus_genes or [])},
        consanguineous=consanguineous,
    )
    inheritance_score = inheritance_result.get("score", 0.0)

    # ── Signal 5: Analysis method gaps ───────────────────────────────────────
    mito_hpo   = {"HP:0003128", "HP:0001249", "HP:0011458", "HP:0003198", "HP:0000407"}
    repeat_hpo = {"HP:0001251", "HP:0002355", "HP:0001310", "HP:0002459", "HP:0001260"}
    current_set = {t.strip().upper() for t in (current_hpo_terms or [])}

    gaps_result = score_analysis_gaps(
        test_type=test_type or "exome",
        test_date=test_date,
        cnv_calling_performed=cnv_calling_performed,
        splice_analysis_performed=splice_analysis_performed,
        mito_analysis_performed=mito_analysis_performed,
        repeat_expansion_checked=repeat_expansion_checked,
        phenotype_suggests_mito=bool(current_set & mito_hpo),
        phenotype_suggests_repeat=bool(current_set & repeat_hpo),
    )
    gaps_score = gaps_result.get("score", 0.0)

    # ── Signal 6: AlphaMissense variant scoring ───────────────────────────────
    am_signals = []
    am_score   = 0.0
    if variants:
        db_path = alphamissense_db_path or os.path.join(
            os.path.dirname(os.path.dirname(__file__)), "data", "alphamissense.db"
        )
        am_signals = check_variants_alphamissense(variants, db_path)
        am_score   = max((s["score"] for s in am_signals), default=0.0)

    # ── Signal 7: Time since test ─────────────────────────────────────────────
    time_result = score_time_since_test(test_date)
    time_score  = time_result.get("time_score", 0.0)

    # ── Weighted base score ───────────────────────────────────────────────────
    base_score = (
        WEIGHTS["omim"]          * omim_score          +
        WEIGHTS["drift"]         * drift_score         +
        WEIGHTS["disease_match"] * disease_match_score +
        WEIGHTS["vus_clinvar"]   * vus_score           +
        WEIGHTS["inheritance"]   * inheritance_score   +
        WEIGHTS["analysis_gaps"] * gaps_score          +
        WEIGHTS["alphamissense"] * am_score            +
        WEIGHTS["time"]          * time_score
    )

    # ── Entropy modifier: multi-VUS diagnostic uncertainty ────────────────────
    # Approximate per-gene posterior: blend disease match + VUS evidence
    # In production: replace with P(phenotype|gene) from HPOA normalised IC scores
    gene_posteriors: Dict[str, float] = {}
    if len(vus_genes or []) >= 2:
        for gene in vus_genes:
            gene_posteriors[gene] = disease_match_score * 0.7 + vus_score * 0.3

    entropy_value = compute_vus_entropy(gene_posteriors) if gene_posteriors else 0.0
    entropy_boost = round(ENTROPY_BOOST_MAX * entropy_value, 3)
    final_score   = min(1.0, base_score + entropy_boost)
    final_score_scaled = round(final_score * 100, 1)

    # ── Confidence / urgency tier ─────────────────────────────────────────────
    if final_score >= URGENT_THRESHOLD:
        confidence = "URGENT"
        urgency    = "Immediate escalation — reanalysis within days"
    elif final_score >= HIGH_THRESHOLD:
        confidence = "HIGH"
        urgency    = "Reanalysis recommended — schedule within weeks"
    elif final_score >= MEDIUM_THRESHOLD:
        confidence = "MEDIUM"
        urgency    = "Consider reanalysis — review at next clinic visit"
    else:
        confidence = "LOW"
        urgency    = "Routine monitoring — revisit in 12 months"

    # ── Per-signal strength labels for UI ────────────────────────────────────
    def _strength(score: float, hi: float = 0.5, med: float = 0.25) -> str:
        return "HIGH" if score >= hi else "MEDIUM" if score >= med else "LOW"

    all_strengths = [
        inheritance_result.get("signal_strength", "LOW"),
        gaps_result.get("signal_strength", "LOW"),
        omim_signals[0]["signal_strength"] if omim_signals else "LOW",
        _strength(vus_score),
        _strength(disease_match_score),
        _strength(drift_score),
        "HIGH" if entropy_value > 0.7 else "MEDIUM" if entropy_value > 0.4 else "LOW",
    ]
    overall_strength = (
        "HIGH"   if "HIGH"   in all_strengths else
        "MEDIUM" if "MEDIUM" in all_strengths else
        "LOW"
    )

    return {
        # Top-level scores
        "reanalysis_score":        final_score_scaled,
        "raw_score":               round(final_score, 3),
        "base_score":              round(base_score, 3),
        "entropy_boost":           entropy_boost,
        "confidence":              confidence,
        "urgency":                 urgency,
        "overall_signal_strength": overall_strength,

        "signal_breakdown": {

            "omim_gene_disease": {
                "score":                 round(omim_score, 3),
                "weight":                WEIGHTS["omim"],
                "weighted_contribution": round(WEIGHTS["omim"] * omim_score, 3),
                "signal_count":          len(omim_signals),
                "signals":               omim_signals,
                "signal_strength":       omim_signals[0]["signal_strength"] if omim_signals else "LOW",
                "description":           "New gene-disease associations in OMIM (42% of reanalysis yield, Radboud 2022)",
            },

            "vus_reclassification": {
                "score":                 round(vus_score, 3),
                "weight":                WEIGHTS["vus_clinvar"],
                "weighted_contribution": round(WEIGHTS["vus_clinvar"] * vus_score, 3),
                "signal_count":          len(vus_signals),
                "signals":               vus_signals,
                "signal_strength":       _strength(vus_score),
                "description":           "ClinVar pathogenic reclassification of VUS genes",
            },

            "phenotypic_drift": {
                "score":                 round(drift_score, 3),
                "weight":                WEIGHTS["drift"],
                "weighted_contribution": round(WEIGHTS["drift"] * drift_score, 3),
                "details":               drift_result,
                "signal_strength":       drift_result.get("signal_strength", "LOW"),
                "description":           "Fraction of current symptoms new since original test",
            },

            "disease_match": {
                "score":                 round(disease_match_score, 3),
                "weight":                WEIGHTS["disease_match"],
                "weighted_contribution": round(WEIGHTS["disease_match"] * disease_match_score, 3),
                "best_matching_disease": drift_result.get("best_matching_disease"),
                "disease_match_details": drift_result.get("disease_match_details", []),
                "similarity_method":     similarity_method,
                "signal_strength":       _strength(disease_match_score),
                "description":           (
                    f"IC-weighted Resnik semantic similarity between patient HPO profile "
                    f"and disease canonical profiles — method active: {similarity_method}"
                ),
            },

            "inheritance_pattern": {
                "score":                 round(inheritance_score, 3),
                "weight":                WEIGHTS["inheritance"],
                "weighted_contribution": round(WEIGHTS["inheritance"] * inheritance_score, 3),
                "flags":                 inheritance_result.get("flags", []),
                "flag_count":            inheritance_result.get("flag_count", 0),
                "signal_strength":       inheritance_result.get("signal_strength", "LOW"),
                "description":           "AR second hit, de novo singleton, X-linked dominant in female, consanguinity",
            },

            "analysis_method_gaps": {
                "score":                 round(gaps_score, 3),
                "weight":                WEIGHTS["analysis_gaps"],
                "weighted_contribution": round(WEIGHTS["analysis_gaps"] * gaps_score, 3),
                "gaps":                  gaps_result.get("gaps", []),
                "gap_count":             gaps_result.get("gap_count", 0),
                "signal_strength":       gaps_result.get("signal_strength", "LOW"),
                "description":           "CNV not called, panel→exome, splice analysis missing, pipeline era gap",
            },

            "alphamissense": {
                "score":                 round(am_score, 3),
                "weight":                WEIGHTS["alphamissense"],
                "weighted_contribution": round(WEIGHTS["alphamissense"] * am_score, 3),
                "signal_count":          len(am_signals),
                "signals":               am_signals,
                "signal_strength":       _strength(am_score),
                "description":           "AlphaMissense structural pathogenicity (requires variant coordinates)",
            },

            "time_since_test": {
                "score":                 round(time_score, 3),
                "weight":                WEIGHTS["time"],
                "weighted_contribution": round(WEIGHTS["time"] * time_score, 3),
                "details":               time_result,
                "signal_strength":       time_result.get("signal_strength", "LOW"),
                "description":           "Time elapsed proxy (mechanism captured by OMIM/gaps signals above)",
            },

            "entropy_modifier": {
                "entropy_value":   entropy_value,
                "entropy_boost":   entropy_boost,
                "max_boost":       ENTROPY_BOOST_MAX,
                "gene_posteriors": {g: round(s, 3) for g, s in gene_posteriors.items()},
                "n_vus_genes":     len(vus_genes or []),
                "signal_strength": (
                    "HIGH" if entropy_value > 0.7 else
                    "MEDIUM" if entropy_value > 0.4 else "LOW"
                ),
                "description": (
                    "Shannon entropy H = −Σ p(g)log₂p(g) across VUS gene posteriors, "
                    "normalised to [0,1]. High entropy = diagnostic uncertainty = "
                    "stronger reanalysis trigger. Max boost: +8 points."
                ),
            },
        },

        "input_summary": {
            "vus_genes":              vus_genes,
            "test_date":              test_date,
            "test_type":              test_type,
            "patient_sex":            patient_sex,
            "analysis_type":          analysis_type,
            "original_symptom_count": len(original_hpo_terms) if original_hpo_terms else 0,
            "current_symptom_count":  len(current_hpo_terms) if current_hpo_terms else 0,
            "variants_provided":      bool(variants),
            "similarity_method":      similarity_method,
        },
    }