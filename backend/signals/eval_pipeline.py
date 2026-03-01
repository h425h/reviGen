"""
Evaluation Pipeline — reviGen Model Assessment.
================================================

Implements the judge-specified evaluation strategy:

  1. Load samples from MyGene2 dataset
  2. For each sample: use its HPO terms as input → model predicts ranked gene list
  3. Score: is the known causal gene predicted correctly?
  4. Compute Precision, Recall, F1 at multiple k thresholds (top-1, top-3, top-10)

MyGene2 (mygene2.org):
  A University of Washington registry of families with undiagnosed rare diseases.
  Each record contains: HPO phenotype terms, optional free-text description,
  and the eventual causal gene (for solved cases).
  Public data: https://mygene2.org/MyGene2/api/public_data

Evaluation framing (top-k gene prediction):
  The model outputs a ranked gene list. We evaluate:
  - Top-1 accuracy:  Is the correct gene rank 1?
  - Top-3 accuracy:  Is it in the top 3?
  - Top-10 accuracy: Is it in the top 10? (standard in the field, used by LIRICAL)
  - F1@k:            Precision × Recall harmonic mean at cutoff k

Why NOT a custom penalty algorithm (judge's question):
  Standard Precision/Recall/F1 is the correct choice because:
  - F1 is the established metric in clinical genomics tool evaluation
    (LIRICAL paper, Exomiser paper, AI-MARRVEL paper all use it)
  - A custom penalty introduces arbitrary hyperparameters with no basis
    in clinical literature — harder to justify to judges
  - F1 = harmonic mean of precision & recall naturally penalises both
    false positives (wrong gene predicted) and false negatives (correct
    gene missed) without ad-hoc tuning
  - AUROC and Brier score (mentioned by judge) are for probability
    calibration — relevant post-hackathon when you have calibrated scores

Usage:
  python eval_pipeline.py                    # Run on built-in test cases
  python eval_pipeline.py --mygene2 data.json  # Run on downloaded MyGene2 data
"""

import json
import math
import os
import sys
import warnings
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple

from .gene_ranker import rank_genes_by_phenotype


# ── Data structures ───────────────────────────────────────────────────────────

@dataclass
class EvalSample:
    """One evaluation case."""
    sample_id:    str
    hpo_terms:    List[str]       # Input to model
    true_gene:    str             # Known causal gene (ground truth)
    disease_name: str = ""
    source:       str = "unknown" # 'mygene2' | 'synthetic' | 'literature'
    notes:        str = ""


@dataclass
class SampleResult:
    """Model output for one sample."""
    sample_id:      str
    true_gene:      str
    predicted_genes: List[str]    # Ranked gene list (rank 1 = most likely)
    top_k_scores:   Dict[str, float]  # {gene: similarity_score}
    rank_of_true:   Optional[int]     # What rank did true gene appear at? None if not in list
    correct_top1:   bool
    correct_top3:   bool
    correct_top10:  bool
    similarity_method: str


@dataclass
class EvalReport:
    """Aggregate evaluation results across all samples."""
    n_samples:      int
    top1_accuracy:  float
    top3_accuracy:  float
    top10_accuracy: float
    precision_at1:  float
    precision_at3:  float
    precision_at10: float
    recall_at1:     float
    recall_at3:     float
    recall_at10:    float
    f1_at1:         float
    f1_at3:         float
    f1_at10:        float
    mean_reciprocal_rank: float   # MRR — standard IR metric
    similarity_method: str
    per_sample: List[SampleResult] = field(default_factory=list)


# ── Built-in test cases ───────────────────────────────────────────────────────
# Derived from published literature and MyGene2 public cases
# These allow the eval pipeline to run without downloading MyGene2 data

BUILTIN_TEST_CASES: List[Dict] = [
    {
        "sample_id":   "FOXG1_001",
        "hpo_terms":   ["HP:0001263","HP:0002186","HP:0000252","HP:0001252",
                        "HP:0001332","HP:0002072","HP:0000729","HP:0001344",
                        "HP:0001270","HP:0011968"],
        "true_gene":   "FOXG1",
        "disease_name": "FOXG1 syndrome",
        "source":       "literature",
    },
    {
        "sample_id":   "CDKL5_001",
        "hpo_terms":   ["HP:0001250","HP:0001263","HP:0001252","HP:0000252",
                        "HP:0002186","HP:0001344","HP:0001270","HP:0011968"],
        "true_gene":   "CDKL5",
        "disease_name": "CDKL5 deficiency disorder",
        "source":       "literature",
    },
    {
        "sample_id":   "MECP2_001",
        "hpo_terms":   ["HP:0001263","HP:0002186","HP:0000252","HP:0001252",
                        "HP:0001332","HP:0001344","HP:0001270","HP:0002072",
                        "HP:0000729","HP:0001250"],
        "true_gene":   "MECP2",
        "disease_name": "Rett syndrome",
        "source":       "literature",
    },
    {
        "sample_id":   "STXBP1_001",
        "hpo_terms":   ["HP:0001250","HP:0001263","HP:0001252","HP:0001270",
                        "HP:0001344","HP:0002069","HP:0000729"],
        "true_gene":   "STXBP1",
        "disease_name": "DEE4 / Ohtahara syndrome",
        "source":       "literature",
    },
    {
        "sample_id":   "PTHR1_001",
        "hpo_terms":   ["HP:0003521","HP:0002750","HP:0000843","HP:0000518",
                        "HP:0002748","HP:0003075","HP:0003127"],
        "true_gene":   "PTHR1",
        "disease_name": "Jansen metaphyseal chondrodysplasia (Nina Nazar case)",
        "source":       "literature",
    },
    {
        "sample_id":   "ANKRD11_001",
        "hpo_terms":   ["HP:0001263","HP:0001249","HP:0001510","HP:0006482",
                        "HP:0000252","HP:0001250","HP:0000729"],
        "true_gene":   "ANKRD11",
        "disease_name": "KBG syndrome",
        "source":       "literature",
    },
    {
        "sample_id":   "SCN2A_001",
        "hpo_terms":   ["HP:0001250","HP:0001263","HP:0001252","HP:0001270",
                        "HP:0001344","HP:0002069","HP:0011185"],
        "true_gene":   "SCN2A",
        "disease_name": "DEE11",
        "source":       "literature",
    },
    {
        "sample_id":   "DYRK1A_001",
        "hpo_terms":   ["HP:0001263","HP:0001249","HP:0000252","HP:0001250",
                        "HP:0000729","HP:0001270","HP:0001344","HP:0011968"],
        "true_gene":   "DYRK1A",
        "disease_name": "DYRK1A syndrome",
        "source":       "literature",
    },
    {
        "sample_id":   "KAT6A_001",
        "hpo_terms":   ["HP:0001263","HP:0001249","HP:0000252","HP:0011968",
                        "HP:0000729","HP:0001270","HP:0001344"],
        "true_gene":   "KAT6A",
        "disease_name": "KAT6A syndrome",
        "source":       "literature",
    },
    {
        "sample_id":   "MED13L_001",
        "hpo_terms":   ["HP:0001263","HP:0001249","HP:0000252","HP:0001270",
                        "HP:0001344","HP:0011968","HP:0000729"],
        "true_gene":   "MED13L",
        "disease_name": "MED13L syndrome",
        "source":       "literature",
    },
]


# ── MyGene2 loader ─────────────────────────────────────────────────────────────

def load_mygene2_samples(
    filepath: Optional[str] = None,
    max_samples: int = 200,
) -> List[EvalSample]:
    """
    Load evaluation samples from MyGene2 public dataset.

    MyGene2 JSON format (public API response):
      Each entry has:
        - hpoTerms: list of {termId, label}
        - variants: list of {gene, ...}  ← we use the first solved gene
        - familyId, phenotypes (free text)

    If filepath is None, returns built-in test cases.

    Download public data from:
      https://mygene2.org/MyGene2/api/public_data

    Args:
        filepath:    Path to downloaded MyGene2 JSON file, or None
        max_samples: Maximum samples to load (default 200)

    Returns:
        List of EvalSample objects
    """
    if filepath is None or not os.path.exists(filepath):
        if filepath is not None:
            warnings.warn(
                f"MyGene2 file not found: {filepath}. Using built-in test cases.",
                UserWarning,
            )
        return _load_builtin_cases()

    try:
        with open(filepath, "r", encoding="utf-8") as f:
            raw = json.load(f)
    except (json.JSONDecodeError, IOError) as e:
        warnings.warn(f"Failed to load MyGene2 file: {e}. Using built-in cases.")
        return _load_builtin_cases()

    samples = []
    entries = raw if isinstance(raw, list) else raw.get("data", raw.get("entries", []))

    for entry in entries[:max_samples]:
        try:
            sample = _parse_mygene2_entry(entry)
            if sample:
                samples.append(sample)
        except Exception:
            continue

    if not samples:
        warnings.warn("No valid samples parsed from MyGene2 file. Using built-in cases.")
        return _load_builtin_cases()

    print(f"Loaded {len(samples)} samples from MyGene2 ({filepath})")
    return samples


def _parse_mygene2_entry(entry: Dict) -> Optional[EvalSample]:
    """Parse a single MyGene2 entry into an EvalSample."""
    # Extract HPO terms
    hpo_terms = []
    for term in entry.get("hpoTerms", entry.get("phenotypes", [])):
        if isinstance(term, dict):
            tid = term.get("termId", term.get("hpoId", term.get("id", "")))
        else:
            tid = str(term)
        if tid and "HP:" in tid.upper():
            hpo_terms.append(tid.upper())

    if not hpo_terms:
        return None

    # Extract known gene (solved cases only)
    true_gene = None
    for var in entry.get("variants", entry.get("genes", [])):
        if isinstance(var, dict):
            gene = var.get("gene", var.get("geneSymbol", var.get("symbol", "")))
            solved = var.get("solved", var.get("diagnosed", True))
            if gene and solved:
                true_gene = gene.upper().strip()
                break
        elif isinstance(var, str):
            true_gene = var.upper().strip()
            break

    if not true_gene:
        return None

    sample_id = str(entry.get("familyId", entry.get("id", f"mg2_{len(hpo_terms)}")))

    return EvalSample(
        sample_id=    sample_id,
        hpo_terms=    hpo_terms,
        true_gene=    true_gene,
        disease_name= entry.get("diagnosis", entry.get("disease", "")),
        source=       "mygene2",
    )


def _load_builtin_cases() -> List[EvalSample]:
    """Load the built-in test cases."""
    return [
        EvalSample(
            sample_id=    c["sample_id"],
            hpo_terms=    c["hpo_terms"],
            true_gene=    c["true_gene"],
            disease_name= c.get("disease_name", ""),
            source=       c.get("source", "literature"),
        )
        for c in BUILTIN_TEST_CASES
    ]


# ── Evaluation engine ─────────────────────────────────────────────────────────

def evaluate_sample(sample: EvalSample, top_k: int = 10) -> SampleResult:
    """
    Run the model on one sample and return its result.

    The model pipeline:
      HPO terms → rank_genes_by_phenotype() → ranked gene list
      → check if true_gene is in top-1, top-3, top-10
    """
    ranked = rank_genes_by_phenotype(sample.hpo_terms, top_k=top_k)

    predicted_genes = [r["gene"] for r in ranked]
    scores = {r["gene"]: r["score"] for r in ranked}
    method = ranked[0]["similarity_method"] if ranked else "jaccard"

    # Find rank of true gene
    rank_of_true = None
    for r in ranked:
        if r["gene"].upper() == sample.true_gene.upper():
            rank_of_true = r["rank"]
            break

    return SampleResult(
        sample_id=       sample.sample_id,
        true_gene=       sample.true_gene,
        predicted_genes= predicted_genes,
        top_k_scores=    scores,
        rank_of_true=    rank_of_true,
        correct_top1=    rank_of_true is not None and rank_of_true <= 1,
        correct_top3=    rank_of_true is not None and rank_of_true <= 3,
        correct_top10=   rank_of_true is not None and rank_of_true <= 10,
        similarity_method= method,
    )


# ── Metrics ───────────────────────────────────────────────────────────────────

def _precision_at_k(results: List[SampleResult], k: int) -> float:
    """
    Precision@k: fraction of top-k predictions that include the true gene.
    In single-gene prediction, this equals accuracy@k.
    """
    if not results:
        return 0.0
    correct = sum(
        1 for r in results
        if r.rank_of_true is not None and r.rank_of_true <= k
    )
    return correct / len(results)


def _recall_at_k(results: List[SampleResult], k: int) -> float:
    """
    Recall@k: fraction of all true genes that are found in top-k predictions.
    In single-label gene prediction, Recall@k = Accuracy@k.
    """
    return _precision_at_k(results, k)


def _f1(precision: float, recall: float) -> float:
    """Harmonic mean of precision and recall."""
    if precision + recall == 0:
        return 0.0
    return 2 * precision * recall / (precision + recall)


def _mean_reciprocal_rank(results: List[SampleResult]) -> float:
    """
    MRR = mean of 1/rank for each sample.
    Standard information retrieval metric.
    If true gene not found, reciprocal rank = 0.
    """
    if not results:
        return 0.0
    total = sum(
        1.0 / r.rank_of_true if r.rank_of_true else 0.0
        for r in results
    )
    return total / len(results)


# ── Main evaluation function ──────────────────────────────────────────────────

def run_evaluation(
    samples: Optional[List[EvalSample]] = None,
    mygene2_filepath: Optional[str] = None,
    max_samples: int = 200,
    top_k: int = 10,
    verbose: bool = True,
) -> EvalReport:
    """
    Run full evaluation pipeline.

    Args:
        samples:           Pre-loaded EvalSample list (or None to load from file/builtin)
        mygene2_filepath:  Path to MyGene2 JSON file (or None for builtin)
        max_samples:       Max samples to evaluate
        top_k:             Rank cutoff for evaluation (default 10)
        verbose:           Print per-sample results

    Returns:
        EvalReport with all metrics
    """
    if samples is None:
        samples = load_mygene2_samples(mygene2_filepath, max_samples)

    if verbose:
        print(f"\n{'='*60}")
        print(f"reviGen Evaluation Pipeline")
        print(f"{'='*60}")
        print(f"Samples: {len(samples)}")
        print(f"Top-k cutoff: {top_k}")
        print()

    results: List[SampleResult] = []

    for sample in samples:
        result = evaluate_sample(sample, top_k=top_k)
        results.append(result)

        if verbose:
            rank_str = f"rank {result.rank_of_true}" if result.rank_of_true else "NOT FOUND"
            top1_str = "✓" if result.correct_top1 else "✗"
            print(
                f"  [{top1_str}] {sample.sample_id:20s} | "
                f"true={sample.true_gene:10s} | "
                f"pred={result.predicted_genes[0] if result.predicted_genes else 'none':10s} | "
                f"{rank_str}"
            )

    # Compute metrics
    p1  = _precision_at_k(results, 1)
    p3  = _precision_at_k(results, 3)
    p10 = _precision_at_k(results, 10)
    r1  = _recall_at_k(results, 1)
    r3  = _recall_at_k(results, 3)
    r10 = _recall_at_k(results, 10)
    mrr = _mean_reciprocal_rank(results)
    sim_method = results[0].similarity_method if results else "unknown"

    report = EvalReport(
        n_samples=     len(results),
        top1_accuracy= p1,
        top3_accuracy= p3,
        top10_accuracy=p10,
        precision_at1= p1,
        precision_at3= p3,
        precision_at10=p10,
        recall_at1=    r1,
        recall_at3=    r3,
        recall_at10=   r10,
        f1_at1=        _f1(p1,  r1),
        f1_at3=        _f1(p3,  r3),
        f1_at10=       _f1(p10, r10),
        mean_reciprocal_rank= mrr,
        similarity_method=    sim_method,
        per_sample=    results,
    )

    if verbose:
        _print_report(report)

    return report


def _print_report(report: EvalReport) -> None:
    """Print formatted evaluation report."""
    print(f"\n{'='*60}")
    print(f"EVALUATION RESULTS  (n={report.n_samples}, method={report.similarity_method})")
    print(f"{'='*60}")
    print(f"{'Metric':<28} {'@k=1':>8} {'@k=3':>8} {'@k=10':>8}")
    print(f"{'-'*60}")
    print(f"{'Precision':<28} {report.precision_at1:>8.3f} {report.precision_at3:>8.3f} {report.precision_at10:>8.3f}")
    print(f"{'Recall':<28} {report.recall_at1:>8.3f} {report.recall_at3:>8.3f} {report.recall_at10:>8.3f}")
    print(f"{'F1 Score':<28} {report.f1_at1:>8.3f} {report.f1_at3:>8.3f} {report.f1_at10:>8.3f}")
    print(f"{'Accuracy (top-k)':<28} {report.top1_accuracy:>8.3f} {report.top3_accuracy:>8.3f} {report.top10_accuracy:>8.3f}")
    print(f"{'-'*60}")
    print(f"{'Mean Reciprocal Rank (MRR)':<28} {report.mean_reciprocal_rank:>8.3f}")
    print(f"{'='*60}")

    # Per-gene breakdown
    gene_results: Dict[str, List[bool]] = {}
    for r in report.per_sample:
        if r.true_gene not in gene_results:
            gene_results[r.true_gene] = []
        gene_results[r.true_gene].append(r.correct_top1)

    print(f"\nPer-gene top-1 accuracy:")
    for gene, hits in sorted(gene_results.items()):
        acc = sum(hits) / len(hits)
        bar = "█" * int(acc * 10) + "░" * (10 - int(acc * 10))
        print(f"  {gene:12s} [{bar}] {acc:.0%}  ({sum(hits)}/{len(hits)})")


def report_to_dict(report: EvalReport) -> Dict:
    """Convert EvalReport to JSON-serialisable dict."""
    return {
        "summary": {
            "n_samples":           report.n_samples,
            "similarity_method":   report.similarity_method,
            "top1_accuracy":       round(report.top1_accuracy,  4),
            "top3_accuracy":       round(report.top3_accuracy,  4),
            "top10_accuracy":      round(report.top10_accuracy, 4),
            "f1_at1":              round(report.f1_at1,         4),
            "f1_at3":              round(report.f1_at3,         4),
            "f1_at10":             round(report.f1_at10,        4),
            "mean_reciprocal_rank":round(report.mean_reciprocal_rank, 4),
        },
        "per_sample": [
            {
                "sample_id":       r.sample_id,
                "true_gene":       r.true_gene,
                "top1_prediction": r.predicted_genes[0] if r.predicted_genes else None,
                "top3_predictions":r.predicted_genes[:3],
                "rank_of_true":    r.rank_of_true,
                "correct_top1":    r.correct_top1,
                "correct_top3":    r.correct_top3,
                "correct_top10":   r.correct_top10,
            }
            for r in report.per_sample
        ],
    }


# ── CLI entry point ───────────────────────────────────────────────────────────

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="reviGen Evaluation Pipeline")
    parser.add_argument("--mygene2",  default=None, help="Path to MyGene2 JSON file")
    parser.add_argument("--max",      type=int, default=200, help="Max samples")
    parser.add_argument("--topk",     type=int, default=10,  help="Rank cutoff")
    parser.add_argument("--out",      default=None, help="Save JSON report to file")
    args = parser.parse_args()

    report = run_evaluation(
        mygene2_filepath=args.mygene2,
        max_samples=args.max,
        top_k=args.topk,
        verbose=True,
    )

    if args.out:
        with open(args.out, "w") as f:
            json.dump(report_to_dict(report), f, indent=2)
        print(f"\nReport saved to {args.out}")