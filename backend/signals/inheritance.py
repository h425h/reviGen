"""
Inheritance Pattern & Family Structure Signal — NEW SIGNAL.

Entirely absent from current reviGen implementation. This signal addresses
several major sources of missed diagnoses:

1. AR disease + single VUS = probable SECOND HIT MISSED
   - Autosomal recessive disease requires 2 pathogenic variants (in trans)
   - If patient has 1 known pathogenic/VUS in an AR gene, the second hit
     (deep intronic splice variant, CNV, maternal deletion) is frequently missed
   - Singleton analysis (no parents) cannot phase variants → both hits ambiguous

2. X-linked variants in females
   - Standard pipelines assume X-linked = hemizygous males only
   - MECP2 het LoF → Rett syndrome in females (classic X-linked dominant)
   - Pipelines often label heterozygous X-linked variants as 'carrier status'
     and exclude from pathogenic assessment

3. Consanguinity / founder populations
   - gnomAD AF computed in outbred populations
   - Variant with AF=0.01 in gnomAD is NOT evidence of benignity for AR disease
     in consanguineous family (Bedouin, Ashkenazi, Finnish, etc.)

4. Singleton vs. trio analysis
   - Trio (proband + both parents) enables de novo variant calling
   - De novo variants in neurodevelopmental disease genes are highly
     enriched for pathogenicity (odds ratio ~10-50x vs inherited variants)
   - 'Negative' singleton analysis may simply have missed a de novo variant
     that trio analysis would confirm

OMIM inheritance modes:
  AR  = Autosomal recessive
  AD  = Autosomal dominant  
  XLR = X-linked recessive
  XLD = X-linked dominant
  MT  = Mitochondrial
"""

from typing import List, Dict, Optional

# Known inheritance modes per gene (OMIM-derived, abbreviated for demo)
# In production: query OMIM API or use a full OMIM gene-phenotype table
GENE_INHERITANCE_MAP = {
    # AR genes (recessive — single VUS is a major trigger)
    "FOXG1":   "AD",    # De novo dominant; important for singleton flag
    "CDKL5":   "XLD",   # X-linked dominant; het females affected
    "SETBP1":  "AD",
    "ANKRD11": "AD",
    "PTHR1":   "AD",
    "STXBP1":  "AD",
    "KCNQ2":   "AD",
    "SCN2A":   "AD",
    "DYRK1A":  "AD",
    "MED13L":  "AD",
    "SETD5":   "AD",
    "KAT6A":   "AD",
    "MECP2":   "XLD",   # X-linked dominant; het females → Rett syndrome
    "PCDH19":  "XLD",   # X-linked; het females affected, hemizygous males spared
    "DHTKD1":  "AR",    # 2-oxoadipic acidemia; second hit often missed
    "PAH":     "AR",    # Phenylketonuria
    "HEX A":   "AR",    # Tay-Sachs
    "HEXA":    "AR",
    "CFTR":    "AR",    # Cystic fibrosis
    "SMN1":    "AR",    # SMA
    "FANCA":   "AR",    # Fanconi anemia
    "BRCA2":   "AR",    # Fanconi anemia (biallelic); AD breast cancer (mono)
    "GBA":     "AR",    # Gaucher disease
    "ARSA":    "AR",    # Metachromatic leukodystrophy
    "CLN3":    "AR",    # Batten disease
    "ASPA":    "AR",    # Canavan disease
    "GALC":    "AR",    # Krabbe disease
    "SMPD1":   "AR",    # Niemann-Pick type A/B
    "NPC1":    "AR",    # Niemann-Pick type C
    "ABCD1":   "XLR",   # Adrenoleukodystrophy (X-linked recessive)
    "F8":      "XLR",   # Haemophilia A
    "F9":      "XLR",   # Haemophilia B
    "DMD":     "XLR",   # Duchenne/Becker muscular dystrophy
    "RPGR":    "XLR",   # Retinitis pigmentosa GTPase regulator
}


def score_inheritance_signal(
    vus_genes: List[str],
    patient_sex: Optional[str] = None,           # 'M', 'F', or None
    analysis_type: Optional[str] = None,          # 'singleton', 'duo', 'trio'
    known_vus_count_per_gene: Optional[Dict[str, int]] = None,  # {gene: count}
    consanguineous: Optional[bool] = None,        # True/False/None
    vus_chromosomes: Optional[Dict[str, str]] = None,  # {gene: 'X'|'auto'}
) -> Dict:
    """
    Score inheritance-based reanalysis triggers.
    
    Args:
        vus_genes:               Genes with VUS on original test
        patient_sex:             'M', 'F', or None if unknown
        analysis_type:           Type of analysis performed: 'singleton'|'duo'|'trio'
        known_vus_count_per_gene: How many VUS/pathogenic found per gene
        consanguineous:          Whether family is consanguineous
        vus_chromosomes:         Chromosome location per gene (for X-linked check)
    
    Returns:
        Signal dict with inheritance-based scores and flags
    """
    flags      = []
    max_score  = 0.0
    strengths  = []

    vus_count = known_vus_count_per_gene or {}

    for gene in (vus_genes or []):
        gene_upper = gene.upper().strip()
        inheritance = GENE_INHERITANCE_MAP.get(gene_upper, "UNKNOWN")

        # --- Flag 1: AR gene with single VUS ---
        # Missing second hit is the most common cause of 'negative' in AR disease
        if inheritance == "AR":
            n_vus = vus_count.get(gene_upper, 1)  # Assume at least 1 if not specified
            if n_vus == 1:
                flags.append({
                    "flag":      "AR_SINGLE_HIT",
                    "gene":      gene,
                    "score":     0.80,
                    "strength":  "HIGH",
                    "reasoning": (
                        f"{gene} causes autosomal recessive disease. Patient has only 1 VUS. "
                        f"Second hit (deep intronic variant, CNV, large deletion) is frequently "
                        f"missed on standard exome sequencing. Trio analysis and CNV calling strongly recommended."
                    ),
                })
                max_score = max(max_score, 0.80)
                strengths.append("HIGH")

            elif n_vus >= 2:
                # Two VUS in same AR gene — possible compound heterozygosity
                flags.append({
                    "flag":      "AR_COMPOUND_HET_CANDIDATE",
                    "gene":      gene,
                    "score":     0.85,
                    "strength":  "HIGH",
                    "reasoning": (
                        f"{gene} causes AR disease. Patient has {n_vus} VUS in this gene. "
                        f"If variants are in TRANS (different alleles, one per parent), "
                        f"this constitutes a complete AR genotype. Parental phase confirmation needed."
                    ),
                })
                max_score = max(max_score, 0.85)
                strengths.append("HIGH")

        # --- Flag 2: X-linked gene in female ---
        # Heterozygous X-linked variants in females are often dismissed as 'carriers'
        # but cause disease in X-linked dominant conditions (MECP2, CDKL5, PCDH19)
        if inheritance in ("XLD", "XLR") and patient_sex == "F":
            flags.append({
                "flag":      "XLINKED_FEMALE",
                "gene":      gene,
                "score":     0.75,
                "strength":  "HIGH",
                "reasoning": (
                    f"{gene} has X-linked inheritance. Female patient with heterozygous variant "
                    f"may be affected if disease is X-linked dominant (e.g. MECP2→Rett, CDKL5→CDD, PCDH19→EIEE). "
                    f"Standard pipelines often label these 'carrier' — clinical correlation needed."
                ),
            })
            max_score = max(max_score, 0.75)
            strengths.append("HIGH")

        # --- Flag 3: De novo AD gene with singleton analysis ---
        # De novo variants are missed in singleton analysis; trio unlocks this
        if inheritance == "AD" and analysis_type == "singleton":
            flags.append({
                "flag":      "AD_DE_NOVO_MISSED",
                "gene":      gene,
                "score":     0.65,
                "strength":  "MEDIUM",
                "reasoning": (
                    f"{gene} causes autosomal dominant disease, commonly via de novo mutations. "
                    f"Singleton analysis cannot confirm de novo status. Trio analysis "
                    f"(proband + both parents) dramatically increases diagnostic yield for "
                    f"neurodevelopmental AD genes."
                ),
            })
            max_score = max(max_score, 0.65)
            strengths.append("MEDIUM")

    # --- Flag 4: Consanguinity modifier ---
    if consanguineous is True and any(
        GENE_INHERITANCE_MAP.get(g.upper(), "") == "AR" for g in (vus_genes or [])
    ):
        flags.append({
            "flag":      "CONSANGUINITY_AR_MODIFIER",
            "gene":      "MULTIPLE",
            "score":     0.60,
            "strength":  "MEDIUM",
            "reasoning": (
                "Consanguineous family with AR gene VUS. Population allele frequencies "
                "in gnomAD (outbred populations) may underestimate pathogenicity risk. "
                "Homozygous variant in AR gene in consanguineous family is strong "
                "pathogenicity indicator even at moderate population AF."
            ),
        })
        max_score = max(max_score, 0.60)
        strengths.append("MEDIUM")

    # --- Flag 5: Singleton analysis overall ---
    if analysis_type == "singleton" and not any(f["flag"] == "AD_DE_NOVO_MISSED" for f in flags):
        flags.append({
            "flag":      "SINGLETON_LIMITATION",
            "gene":      "ALL",
            "score":     0.40,
            "strength":  "MEDIUM",
            "reasoning": (
                "Original analysis was singleton (proband only). Trio analysis enables: "
                "de novo variant identification, variant phasing (cis vs trans for AR), "
                "and segregation analysis. Significantly increases diagnostic yield."
            ),
        })
        max_score = max(max_score, 0.40)
        strengths.append("MEDIUM")

    # Determine overall signal strength
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
        "signal_type":     "inheritance_pattern",
        "score":           round(max_score, 3),
        "signal_strength": overall_strength,
        "flags":           flags,
        "flag_count":      len(flags),
        "reasoning":       (
            f"{len(flags)} inheritance-based flags raised: "
            + "; ".join(f["flag"] for f in flags)
        ) if flags else "No inheritance-based flags.",
    }