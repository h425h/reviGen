"""
VUS Reclassification Signal Detector — IMPROVED.

Key fix: the original logic checked if ClinVar had a pathogenic entry
*after* the test date. This was broken because most FOXG1/PTHR1 pathogenic
entries in ClinVar predate 2019.

Corrected logic:
1. Does ClinVar confirm this gene has pathogenic variants AT ALL?
2. What is the review quality/consensus strength?
3. Have new entries appeared since the test date?

This mirrors real clinical VUS reanalysis: if a gene is now confirmed
pathogenic in ClinVar with expert evidence, a VUS from years ago warrants
re-examination regardless of when ClinVar was updated.
"""
from datetime import datetime
from typing import List, Dict
import pandas as pd
from ..data_loaders import get_clinvar_data


REVIEW_QUALITY = {
    'practice guideline':                                   5,
    'reviewed by expert panel':                             4,
    'criteria provided, multiple submitters, no conflicts': 3,
    'criteria provided, single submitter':                  2,
    'no assertion criteria provided':                       1,
    'no classification provided':                           0,
}


def _review_rank(review_status: str) -> int:
    rs = str(review_status).lower().strip()
    for key, val in REVIEW_QUALITY.items():
        if key in rs:
            return val
    return 0


def check_vus_reclassification(vus_genes: List[str], test_date: str) -> List[Dict]:
    """
    Check whether VUS genes now have strong pathogenic evidence in ClinVar.

    Three-tier scoring:
    - HIGH:   Expert panel / multiple submitters confirm pathogenic
    - MEDIUM: Single submitter pathogenic, or new entries added after test date
    - LOW:    Only likely pathogenic with weak review quality
    """
    signals = []

    if not vus_genes:
        return signals

    try:
        test_datetime = datetime.strptime(test_date, '%Y-%m-%d')
    except (ValueError, TypeError):
        return signals

    clinvar_df = get_clinvar_data()
    if clinvar_df.empty:
        return signals

    for gene in vus_genes:
        gene_upper = gene.upper().strip()

        gene_mask = clinvar_df['GeneSymbol'].str.upper().str.strip() == gene_upper
        gene_variants = clinvar_df[gene_mask].copy()

        if gene_variants.empty:
            continue

        # Find all pathogenic/likely pathogenic entries (exclude uncertain/benign)
        path_mask = (
            gene_variants['ClinicalSignificance'].str.contains('pathogenic', case=False, na=False) &
            ~gene_variants['ClinicalSignificance'].str.contains('uncertain|benign', case=False, na=False)
        )
        pathogenic_entries = gene_variants[path_mask].copy()

        if pathogenic_entries.empty:
            continue

        # Rank entries by review quality, take best
        pathogenic_entries['review_rank'] = pathogenic_entries['ReviewStatus'].apply(_review_rank)
        best_entry = pathogenic_entries.sort_values('review_rank', ascending=False).iloc[0]

        clin_sig      = str(best_entry.get('ClinicalSignificance', ''))
        review_status = str(best_entry.get('ReviewStatus', ''))
        review_rank   = int(best_entry['review_rank'])
        phenotype     = str(best_entry.get('PhenotypeList', 'Unknown'))

        classification = ('Likely Pathogenic' if 'likely pathogenic' in clin_sig.lower()
                          else 'Pathogenic')

        # Count entries added after test date
        new_entries_count = 0
        last_eval_str = str(best_entry.get('LastEvaluated', ''))

        for _, row in pathogenic_entries.iterrows():
            last_eval = str(row.get('LastEvaluated', ''))
            if pd.notna(last_eval) and last_eval not in ('', 'nan', '-'):
                for fmt in ['%Y-%m-%d', '%b %d, %Y', '%Y/%m/%d', '%m/%d/%Y']:
                    try:
                        dt = datetime.strptime(last_eval.strip(), fmt)
                        if dt > test_datetime:
                            new_entries_count += 1
                        break
                    except ValueError:
                        continue

        # Assign signal strength based on review quality
        if review_rank >= 3:
            signal_strength = 'HIGH'
        elif review_rank == 2 and classification == 'Pathogenic':
            signal_strength = 'MEDIUM'
        else:
            signal_strength = 'LOW'

        # Boost if new evidence appeared after test
        if new_entries_count > 0 and signal_strength == 'MEDIUM':
            signal_strength = 'HIGH'

        signals.append({
            'gene':                     gene,
            'old_classification':       'Uncertain Significance',
            'new_classification':       classification,
            'signal_strength':          signal_strength,
            'review_status':            review_status,
            'review_quality_rank':      review_rank,
            'associated_disease':       phenotype[:200],
            'last_evaluated':           last_eval_str,
            'new_entries_after_test':   new_entries_count,
            'total_pathogenic_entries': len(pathogenic_entries),
            'reasoning': (
                f"ClinVar has {len(pathogenic_entries)} pathogenic "
                f"{'entries' if len(pathogenic_entries) > 1 else 'entry'} for {gene} "
                f"(review: {review_status}). "
                f"{new_entries_count} new {'entries' if new_entries_count != 1 else 'entry'} "
                f"after test date."
            )
        })

    return signals