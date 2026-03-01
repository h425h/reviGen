"""
Analysis Method Gap Signal — NEW SIGNAL.

NEJM 2021 (Eldomery et al.): 10.8% of diagnoses from existing genomic data
required additional analytic methods — primarily CNV calling on exome data.

Standard exome sequencing poorly detects:
  - Copy number variants (CNVs) / large deletions/duplications
  - Deep intronic splice variants (>300bp from exon)  
  - Tandem repeat expansions (Huntington's, FMR1, SCA disorders)
  - Complex structural variants (inversions, translocations)
  - Mitochondrial DNA variants (low heteroplasmy)
  - Long-range regulatory variants (promoter, enhancer)

Pipeline version gaps:
  - BWA-MEM → BWA-MEM2: ~15% improvement in alignment sensitivity at SVs
  - GATK 3.x → 4.x: improved genotyping of homopolymer regions, indels
  - DeepVariant (2018) vs older callers: ~5% increase in SNV/indel accuracy
  - gnomAD v2 → v4: 2x more samples, better AF estimates for interpretation

This signal scores cases where the original analysis method had known
limitations that leave a significant portion of the genome uninterrogated.
"""

from typing import Optional, Dict, List
from datetime import datetime


# Pipeline era definitions (approximate dates of major tool releases)
PIPELINE_ERAS = {
    # BWA
    "bwa_mem_v1": {"released": "2013-01-01", "issues": ["SV detection limited", "STR regions"]},
    "bwa_mem2":   {"released": "2019-01-01", "issues": []},
    
    # GATK
    "gatk_3x":    {"released": "2012-01-01", "issues": ["homopolymer indels", "strand bias artefacts"]},
    "gatk_4x":    {"released": "2018-01-01", "issues": ["minor VQSR tuning"]},
    
    # DeepVariant
    "deepvariant_v1": {"released": "2018-01-01", "issues": ["limited training data"]},
    "deepvariant_v1_6": {"released": "2023-01-01", "issues": []},
}

# Analysis completeness checklist
ANALYSIS_METHODS = {
    "cnv_calling": {
        "description": "Copy number variant detection (ExomeDepth, GATK-gCNV, CNVkit)",
        "missing_score": 0.70,
        "strength": "HIGH",
        "reasoning": (
            "CNV calling was NOT performed on the original data. NEJM 2021: 10.8% of "
            "positive reanalyses required CNV detection. Large exon deletions/duplications "
            "can be called retrospectively from existing exome BAM files using "
            "ExomeDepth, GATK gCNV, or CNVkit without new sequencing."
        ),
    },
    "splice_ai_scoring": {
        "description": "Deep intronic splice variant detection (SpliceAI, CADD-Splice)",
        "missing_score": 0.55,
        "strength": "MEDIUM",
        "reasoning": (
            "Splice variant analysis not performed. SpliceAI (delta score ≥0.2) detects "
            "cryptic splice sites up to 50bp from exon boundaries that standard annotation "
            "misses. Can be applied to existing VCF without re-sequencing."
        ),
    },
    "mito_analysis": {
        "description": "Mitochondrial DNA variant calling (Haplogrep, Mutect2 in MT mode)",
        "missing_score": 0.50,
        "strength": "MEDIUM",
        "reasoning": (
            "Mitochondrial DNA not analysed or analysed with insufficient sensitivity. "
            "Low-level heteroplasmy (5-30%) requires specialised callers. "
            "MT analysis on existing exome data adds diagnostic yield for "
            "neuromuscular, encephalopathy, and metabolic presentations."
        ),
    },
    "repeat_expansion": {
        "description": "Tandem repeat expansion detection (ExpansionHunter, STRipy)",
        "missing_score": 0.60,
        "strength": "HIGH",
        "reasoning": (
            "Repeat expansion analysis not performed. Standard variant callers miss "
            "pathogenic expansions (FMR1→Fragile X, C9orf72→ALS/FTD, HTT→Huntington's, "
            "ATXN→spinocerebellar ataxias). ExpansionHunter can run on existing aligned reads."
        ),
    },
}

# Sequencing modality upgrade potential
MODALITY_UPGRADES = {
    "gene_panel_to_exome": {
        "score": 0.70,
        "strength": "HIGH",
        "reasoning": (
            "Original test was a GENE PANEL. Whole exome sequencing covers "
            "~20,000 genes vs panel's limited scope. Genes not on the original "
            "panel cannot be assessed — exome reanalysis may reveal causal variants "
            "in genes not yet associated with this phenotype at time of testing."
        ),
    },
    "exome_to_genome": {
        "score": 0.55,
        "strength": "MEDIUM",
        "reasoning": (
            "Whole genome sequencing has 7-10% additional diagnostic yield vs exome. "
            "Covers: regulatory regions, deep intronic variants, pseudogenes (e.g. SMN1/SMN2), "
            "structural variants, and GC-rich regions poorly captured by exome capture kits."
        ),
    },
}


def score_analysis_gaps(
    test_type: str,                              # 'exome', 'gene panel', 'genome', 'other'
    test_date: str,                              # YYYY-MM-DD
    cnv_calling_performed: Optional[bool] = None,
    splice_analysis_performed: Optional[bool] = None,
    mito_analysis_performed: Optional[bool] = None,
    repeat_expansion_checked: Optional[bool] = None,
    pipeline_version: Optional[str] = None,     # e.g., 'gatk_3x', 'bwa_mem_v1'
    phenotype_suggests_repeat: bool = False,     # e.g., ataxia, tremor, cognitive decline
    phenotype_suggests_mito: bool = False,       # e.g., myopathy, encephalopathy, lactic acidosis
) -> Dict:
    """
    Score analysis method gaps that may explain a negative result.
    
    Args:
        test_type:                  Sequencing modality used
        test_date:                  Date of original test
        cnv_calling_performed:      Was CNV calling run? None = unknown
        splice_analysis_performed:  Was deep splice analysis done?
        mito_analysis_performed:    Was MT DNA analysed?
        repeat_expansion_checked:   Were STR expansions assessed?
        pipeline_version:           Known pipeline tool version
        phenotype_suggests_repeat:  Does phenotype fit repeat expansion disorders?
        phenotype_suggests_mito:    Does phenotype fit mitochondrial disease?
    
    Returns:
        Signal dict with gap flags and scores
    """
    flags      = []
    max_score  = 0.0
    strengths  = []

    # Normalize test type
    test_type_norm = test_type.lower().strip() if test_type else ''

    # --- Gap 1: Gene panel → should have been exome ---
    if 'panel' in test_type_norm:
        info = MODALITY_UPGRADES["gene_panel_to_exome"]
        flags.append({
            "gap":       "GENE_PANEL_SCOPE_LIMITED",
            "score":     info["score"],
            "strength":  info["strength"],
            "reasoning": info["reasoning"],
        })
        max_score = max(max_score, info["score"])
        strengths.append(info["strength"])

    # --- Gap 2: Exome → genome upgrade ---
    elif 'exome' in test_type_norm:
        info = MODALITY_UPGRADES["exome_to_genome"]
        flags.append({
            "gap":       "EXOME_GENOME_UPGRADE_AVAILABLE",
            "score":     info["score"],
            "strength":  info["strength"],
            "reasoning": info["reasoning"],
        })
        max_score = max(max_score, info["score"])
        strengths.append(info["strength"])

    # --- Gap 3: CNV calling ---
    # If not performed OR unknown (None = we don't know → flag it)
    if cnv_calling_performed is False or (
        cnv_calling_performed is None and 'genome' not in test_type_norm
    ):
        info = ANALYSIS_METHODS["cnv_calling"]
        flags.append({
            "gap":       "CNV_CALLING_NOT_PERFORMED",
            "score":     info["missing_score"],
            "strength":  info["strength"],
            "reasoning": info["reasoning"],
        })
        max_score = max(max_score, info["missing_score"])
        strengths.append(info["strength"])

    # --- Gap 4: Splice analysis ---
    if splice_analysis_performed is False:
        info = ANALYSIS_METHODS["splice_ai_scoring"]
        flags.append({
            "gap":       "SPLICE_ANALYSIS_NOT_PERFORMED",
            "score":     info["missing_score"],
            "strength":  info["strength"],
            "reasoning": info["reasoning"],
        })
        max_score = max(max_score, info["missing_score"])
        strengths.append(info["strength"])

    # --- Gap 5: Mitochondrial analysis (especially relevant to phenotype) ---
    if phenotype_suggests_mito and (
        mito_analysis_performed is False or mito_analysis_performed is None
    ):
        info = ANALYSIS_METHODS["mito_analysis"]
        flags.append({
            "gap":       "MITO_ANALYSIS_MISSING",
            "score":     info["missing_score"],
            "strength":  info["strength"],
            "reasoning": info["reasoning"],
        })
        max_score = max(max_score, info["missing_score"])
        strengths.append(info["strength"])

    # --- Gap 6: Repeat expansion (especially relevant to phenotype) ---
    if phenotype_suggests_repeat and not repeat_expansion_checked:
        info = ANALYSIS_METHODS["repeat_expansion"]
        flags.append({
            "gap":       "REPEAT_EXPANSION_NOT_CHECKED",
            "score":     info["missing_score"],
            "strength":  info["strength"],
            "reasoning": info["reasoning"],
        })
        max_score = max(max_score, info["missing_score"])
        strengths.append(info["strength"])

    # --- Gap 7: Old pipeline version ---
    if pipeline_version and pipeline_version in PIPELINE_ERAS:
        era = PIPELINE_ERAS[pipeline_version]
        issues = era.get("issues", [])
        if issues:
            flags.append({
                "gap":       "OLD_PIPELINE_VERSION",
                "score":     0.40,
                "strength":  "MEDIUM",
                "reasoning": (
                    f"Pipeline version '{pipeline_version}' (released ~{era['released'][:4]}) "
                    f"has known limitations: {', '.join(issues)}. "
                    f"Reprocessing with current tools may reveal previously missed variants."
                ),
            })
            max_score = max(max_score, 0.40)
            strengths.append("MEDIUM")

    # --- Gap 8: Old test — bioinformatics pipeline improvement signal ---
    # Even without knowing the exact pipeline, older tests used older tools
    try:
        test_dt = datetime.strptime(test_date, "%Y-%m-%d")
        years_old = (datetime.now() - test_dt).days / 365.25
        if years_old >= 5 and not pipeline_version:
            flags.append({
                "gap":       "PIPELINE_ERA_GAP",
                "score":     0.45,
                "strength":  "MEDIUM",
                "reasoning": (
                    f"Test performed {years_old:.1f} years ago. Bioinformatics pipelines "
                    f"have substantially improved since then (GATK 3→4, DeepVariant, "
                    f"improved alignment algorithms). Radboud 2022: 19% of new diagnoses "
                    f"attributed to improved bioinformatics pipeline alone."
                ),
            })
            max_score = max(max_score, 0.45)
            strengths.append("MEDIUM")
    except (ValueError, TypeError):
        pass

    # Overall strength
    if "HIGH" in strengths:
        overall_strength = "HIGH"
    elif "MEDIUM" in strengths:
        overall_strength = "MEDIUM"
    elif flags:
        overall_strength = "LOW"
    else:
        overall_strength = "LOW"
        max_score = 0.0

    return {
        "signal_type":     "analysis_method_gaps",
        "score":           round(max_score, 3),
        "signal_strength": overall_strength,
        "gaps":            flags,
        "gap_count":       len(flags),
        "reasoning": (
            f"{len(flags)} analysis gaps identified: "
            + "; ".join(f["gap"] for f in flags)
        ) if flags else "No analysis gaps identified.",
    }