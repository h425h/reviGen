"""
reviGen API — FastAPI Backend — FINAL.
======================================

Endpoints:
  GET  /api/health                — health check + active signals list
  POST /api/reanalysis-score      — core reanalysis scoring (7 signals)
  POST /api/extract-hpo           — NLP: free text → HPO terms (Claude-based)
  POST /api/extract-and-score     — NLP + reanalysis in one call
  POST /api/rank-genes            — phenotype → ranked gene list
  GET  /api/evaluate              — run evaluation pipeline, return F1/metrics
  POST /api/evaluate/single       — evaluate one sample against ground truth
"""

import os
import json
from contextlib import asynccontextmanager
from typing import List, Optional, Dict

from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
import anthropic

from .data_loaders import load_all_datasets, hpo_terms
from .aggregator import compute_reanalysis_score
from .nlp_extractor import extract_hpo_with_metadata
from .gene_ranker import rank_genes_by_phenotype
from .eval_pipeline import (
    run_evaluation, load_mygene2_samples,
    EvalSample, report_to_dict,
)


# ── Pydantic Models ────────────────────────────────────────────────────────────

class VariantInfo(BaseModel):
    chrom: str
    pos:   int
    ref:   str
    alt:   str
    gene:  Optional[str] = None
    hgvs:  Optional[str] = None


class PriorTest(BaseModel):
    test_type:     str
    test_date:     str
    result:        str
    vus_genes:     List[str] = []
    analysis_type: Optional[str] = None  # 'singleton'|'duo'|'trio'


class ClinicalContext(BaseModel):
    patient_sex:               Optional[str]  = None
    consanguineous:            Optional[bool] = None
    cnv_calling_performed:     Optional[bool] = None
    splice_analysis_performed: Optional[bool] = None
    mito_analysis_performed:   Optional[bool] = None
    repeat_expansion_checked:  Optional[bool] = None
    variants:                  Optional[List[VariantInfo]] = None


class ReanalysisRequest(BaseModel):
    original_hpo_terms: List[str]
    current_hpo_terms:  List[str]
    prior_test:         PriorTest
    clinical_context:   Optional[ClinicalContext] = None


class ExtractHpoRequest(BaseModel):
    """Free-text clinical note → HPO terms."""
    clinical_text: str
    context:       Optional[str] = None  # 'original' | 'current' | None


class ExtractAndScoreRequest(BaseModel):
    """NLP extraction + reanalysis scoring in one call."""
    original_clinical_text: str          # Notes from time of original test
    current_clinical_text:  str          # Current clinical notes
    prior_test:             PriorTest
    clinical_context:       Optional[ClinicalContext] = None


class RankGenesRequest(BaseModel):
    """Rank candidate genes from HPO profile."""
    hpo_terms: List[str]
    top_k:     int = 10


class EvalSingleRequest(BaseModel):
    """Evaluate model on one sample with known ground truth."""
    hpo_terms:  List[str]
    true_gene:  str
    sample_id:  Optional[str] = "user_input"


class HealthResponse(BaseModel):
    status:          str
    datasets_loaded: bool
    signals_active:  List[str]
    endpoints:       List[str]


# ── Lifespan ───────────────────────────────────────────────────────────────────

datasets_loaded = False


@asynccontextmanager
async def lifespan(app: FastAPI):
    global datasets_loaded
    data_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")

    print("=" * 60)
    print("reviGen API — FINAL — Starting…")
    print("=" * 60)

    datasets_loaded = load_all_datasets(data_dir)

    print(f"OMIM key:        {'SET' if os.environ.get('OMIM_API_KEY') else 'not set (fallback active)'}")
    print(f"Anthropic key:   {'SET' if os.environ.get('ANTHROPIC_API_KEY') else 'not set (NLP fallback active)'}")
    print(f"Datasets loaded: {datasets_loaded}")
    print("=" * 60)

    yield
    print("reviGen API shutting down.")


# ── App ────────────────────────────────────────────────────────────────────────

app = FastAPI(
    title="reviGen API",
    description="Rare disease genetic reanalysis trigger — 7-signal engine with NLP and evaluation",
    version="3.0.0",
    lifespan=lifespan,
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


# ── 1. Health ──────────────────────────────────────────────────────────────────

@app.get("/api/health", response_model=HealthResponse)
async def health_check():
    return HealthResponse(
        status="healthy",
        datasets_loaded=datasets_loaded,
        signals_active=[
            "OMIM gene-disease surveillance (weight 0.29)",
            "ClinVar VUS reclassification (weight 0.14)",
            "Phenotypic drift — raw (weight 0.11)",
            "Disease profile match — Resnik IC-weighted (weight 0.13)",
            "Inheritance pattern flags (weight 0.12)",
            "Analysis method gaps (weight 0.10)",
            "AlphaMissense variant scoring (weight 0.05, requires variant coords)",
            "Time since test (weight 0.03)",
            "Entropy modifier — multi-VUS uncertainty (up to +0.08)",
        ],
        endpoints=[
            "POST /api/reanalysis-score",
            "POST /api/extract-hpo",
            "POST /api/extract-and-score",
            "POST /api/rank-genes",
            "GET  /api/evaluate",
            "POST /api/evaluate/single",
        ],
    )


# ── 2. Core reanalysis score ───────────────────────────────────────────────────

@app.post("/api/reanalysis-score")
async def compute_score(request: ReanalysisRequest):
    """
    Core endpoint. Takes structured HPO terms + prior test metadata.
    Returns 7-signal reanalysis score with Claude clinical recommendation.
    """
    try:
        ctx = request.clinical_context or ClinicalContext()
        variant_dicts = [v.dict() for v in ctx.variants] if ctx.variants else None

        result = compute_reanalysis_score(
            vus_genes=             request.prior_test.vus_genes,
            test_date=             request.prior_test.test_date,
            original_hpo_terms=    request.original_hpo_terms,
            current_hpo_terms=     request.current_hpo_terms,
            test_type=             request.prior_test.test_type,
            patient_sex=           ctx.patient_sex,
            analysis_type=         request.prior_test.analysis_type,
            consanguineous=        ctx.consanguineous,
            cnv_calling_performed= ctx.cnv_calling_performed,
            splice_analysis_performed= ctx.splice_analysis_performed,
            mito_analysis_performed=   ctx.mito_analysis_performed,
            repeat_expansion_checked=  ctx.repeat_expansion_checked,
            variants=              variant_dicts,
            omim_api_key=          os.environ.get("OMIM_API_KEY", ""),
        )

        result["clinical_recommendation"] = await _get_claude_recommendation(request, result)
        return result

    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


# ── 3. NLP: Free text → HPO terms ─────────────────────────────────────────────

@app.post("/api/extract-hpo")
async def extract_hpo(request: ExtractHpoRequest):
    """
    NLP endpoint. Converts free-text clinical description to HPO terms.

    Uses Claude API for extraction (falls back to keyword matching if
    no API key is set). Output is a list of HP:xxxxxxx strings ready
    to feed into /api/reanalysis-score.

    Example input:
      "7-year-old with hypotonia since birth, absent speech, hand
       stereotypies, and seizures from age 2. MRI shows microcephaly."

    Example output:
      ["HP:0001252", "HP:0001344", "HP:0000733", "HP:0001250", "HP:0000252"]
    """
    if not request.clinical_text.strip():
        raise HTTPException(status_code=400, detail="clinical_text is required")

    try:
        result = extract_hpo_with_metadata(
            clinical_text=request.clinical_text,
            anthropic_api_key=os.environ.get("ANTHROPIC_API_KEY", ""),
        )
        return {
            "hpo_terms":  result["hpo_terms"],
            "term_count": result["term_count"],
            "method":     result["method"],
            "message":    (
                f"Extracted {result['term_count']} HPO terms using {result['method']} method. "
                f"{'(Claude API)' if result['method'] == 'claude' else '(keyword fallback — set ANTHROPIC_API_KEY for better accuracy)'}"
            ),
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


# ── 4. Extract + Score (NLP pipeline in one call) ─────────────────────────────

@app.post("/api/extract-and-score")
async def extract_and_score(request: ExtractAndScoreRequest):
    """
    Full pipeline in one call:
      free text (original notes) → HPO terms
      free text (current notes)  → HPO terms
      → 7-signal reanalysis score

    This is the demo endpoint: paste two clinical notes, get a score.
    """
    try:
        api_key = os.environ.get("ANTHROPIC_API_KEY", "")

        # Extract HPO terms from both sets of notes
        original_result = extract_hpo_with_metadata(request.original_clinical_text, api_key)
        current_result  = extract_hpo_with_metadata(request.current_clinical_text,  api_key)

        original_hpo = original_result["hpo_terms"]
        current_hpo  = current_result["hpo_terms"]

        if not original_hpo and not current_hpo:
            raise HTTPException(
                status_code=422,
                detail="Could not extract HPO terms from either clinical note. "
                       "Try adding more specific phenotype descriptions."
            )

        # Build synthetic ReanalysisRequest and call core scorer
        ctx = request.clinical_context or ClinicalContext()
        variant_dicts = [v.dict() for v in ctx.variants] if ctx.variants else None

        result = compute_reanalysis_score(
            vus_genes=             request.prior_test.vus_genes,
            test_date=             request.prior_test.test_date,
            original_hpo_terms=    original_hpo,
            current_hpo_terms=     current_hpo,
            test_type=             request.prior_test.test_type,
            patient_sex=           ctx.patient_sex,
            analysis_type=         request.prior_test.analysis_type,
            consanguineous=        ctx.consanguineous,
            cnv_calling_performed= ctx.cnv_calling_performed,
            variants=              variant_dicts,
            omim_api_key=          os.environ.get("OMIM_API_KEY", ""),
        )

        # Build a synthetic request for Claude recommendation
        synthetic_request = ReanalysisRequest(
            original_hpo_terms=original_hpo,
            current_hpo_terms= current_hpo,
            prior_test=        request.prior_test,
            clinical_context=  request.clinical_context,
        )
        result["clinical_recommendation"] = await _get_claude_recommendation(
            synthetic_request, result
        )

        # Add NLP extraction metadata
        result["nlp_extraction"] = {
            "original_hpo_terms": original_hpo,
            "current_hpo_terms":  current_hpo,
            "original_method":    original_result["method"],
            "current_method":     current_result["method"],
            "original_term_count":original_result["term_count"],
            "current_term_count": current_result["term_count"],
        }

        return result

    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


# ── 5. Gene ranker ─────────────────────────────────────────────────────────────

@app.post("/api/rank-genes")
async def rank_genes(request: RankGenesRequest):
    """
    Phenotype → ranked gene list.
    Uses disease-intermediate approach: HPO → HPOA disease matching → OMIM gene lookup.
    Returns ranked candidates with similarity scores.
    """
    if not request.hpo_terms:
        raise HTTPException(status_code=400, detail="hpo_terms is required")

    try:
        ranked = rank_genes_by_phenotype(request.hpo_terms, top_k=request.top_k)
        return {
            "ranked_genes":       ranked,
            "gene_count":         len(ranked),
            "similarity_method":  ranked[0]["similarity_method"] if ranked else "none",
            "input_hpo_count":    len(request.hpo_terms),
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


# ── 6. Evaluation: full dataset ────────────────────────────────────────────────

@app.get("/api/evaluate")
async def evaluate(
    mygene2_path: Optional[str] = Query(None, description="Path to MyGene2 JSON file on server"),
    max_samples:  int = Query(50, description="Max samples to evaluate"),
    top_k:        int = Query(10, description="Rank cutoff for evaluation"),
):
    """
    Run evaluation pipeline and return F1/Precision/Recall metrics.

    Uses built-in test cases if no MyGene2 file is provided.
    Built-in cases are 10 literature-validated cases covering the most
    common genes in MyGene2 (FOXG1, CDKL5, MECP2, STXBP1, etc.)

    Returns:
      - F1 score at k=1, k=3, k=10
      - Precision and Recall at each k
      - Mean Reciprocal Rank (MRR)
      - Per-sample results with rank of true gene
    """
    try:
        samples = load_mygene2_samples(mygene2_path, max_samples)
        report  = run_evaluation(samples=samples, top_k=top_k, verbose=False)
        return report_to_dict(report)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


# ── 7. Evaluation: single sample ──────────────────────────────────────────────

@app.post("/api/evaluate/single")
async def evaluate_single(request: EvalSingleRequest):
    """
    Evaluate the model on one sample with a known ground truth gene.

    This is the judge-specified evaluation:
      - Input: HPO terms + true gene
      - Output: was the model correct? At what rank?
      - Returns precision/recall/F1 for this single case

    Example:
      POST /api/evaluate/single
      {"hpo_terms": ["HP:0001263","HP:0002186","HP:0000252"], "true_gene": "FOXG1"}
    """
    try:
        sample = EvalSample(
            sample_id= request.sample_id or "user_input",
            hpo_terms= request.hpo_terms,
            true_gene= request.true_gene.upper().strip(),
            source=    "user_input",
        )

        ranked = rank_genes_by_phenotype(request.hpo_terms, top_k=20)
        predicted_genes = [r["gene"] for r in ranked]

        rank_of_true = None
        for r in ranked:
            if r["gene"].upper() == sample.true_gene.upper():
                rank_of_true = r["rank"]
                break

        correct_top1  = rank_of_true is not None and rank_of_true <= 1
        correct_top3  = rank_of_true is not None and rank_of_true <= 3
        correct_top10 = rank_of_true is not None and rank_of_true <= 10

        # For a single sample: precision@k = 1 if correct else 0
        # F1 = 2*P*R/(P+R); with P=R for single-label, F1 = P = R = accuracy
        def f1(correct: bool) -> float:
            return 1.0 if correct else 0.0

        return {
            "sample_id":       request.sample_id,
            "true_gene":       sample.true_gene,
            "rank_of_true_gene": rank_of_true,
            "top_predictions": ranked[:5],

            "scores": {
                "correct_top1":  correct_top1,
                "correct_top3":  correct_top3,
                "correct_top10": correct_top10,
                "f1_at_1":       f1(correct_top1),
                "f1_at_3":       f1(correct_top3),
                "f1_at_10":      f1(correct_top10),
                "precision_at_1": 1.0 if correct_top1  else 0.0,
                "precision_at_3": 1.0 if correct_top3  else 0.0,
                "precision_at_10":1.0 if correct_top10 else 0.0,
                "recall_at_1":    1.0 if correct_top1  else 0.0,
                "recall_at_3":    1.0 if correct_top3  else 0.0,
                "recall_at_10":   1.0 if correct_top10 else 0.0,
                "reciprocal_rank": round(1.0 / rank_of_true, 4) if rank_of_true else 0.0,
            },

            "verdict": (
                f"✓ CORRECT — {sample.true_gene} predicted at rank {rank_of_true}"
                if rank_of_true else
                f"✗ NOT FOUND — {sample.true_gene} was not in the top {len(ranked)} predictions"
            ),

            "similarity_method": ranked[0]["similarity_method"] if ranked else "none",
        }

    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


# ── Claude recommendation (shared by reanalysis endpoints) ─────────────────────

async def _get_claude_recommendation(request: ReanalysisRequest, result: Dict) -> Dict:
    api_key = os.environ.get("ANTHROPIC_API_KEY", "")

    sb           = result.get("signal_breakdown", {})
    years        = sb.get("time_since_test", {}).get("details", {}).get("years_since_test", 0)
    new_symptoms = sb.get("phenotypic_drift", {}).get("details", {}).get("new_symptoms", [])
    omim_signals = sb.get("omim_gene_disease", {}).get("signals", [])
    inh_flags    = sb.get("inheritance_pattern", {}).get("flags", [])
    gap_flags    = sb.get("analysis_method_gaps", {}).get("gaps", [])
    am_signals   = sb.get("alphamissense", {}).get("signals", [])
    vus_signals  = sb.get("vus_reclassification", {}).get("signals", [])

    if not api_key:
        return _fallback_recommendation(result, years, new_symptoms, omim_signals, inh_flags, gap_flags)

    try:
        client = anthropic.Anthropic(api_key=api_key)

        symptom_names = [s.get("name", s.get("hpo_id", "?")) for s in new_symptoms]
        omim_info = [f"{s['gene']}: {s['disease']} — {'NEW after test' if s.get('is_new_after_test') else 'existing'}"
                     for s in omim_signals]
        vus_info  = [f"{s['gene']}: {s['new_classification']} ({s['review_status']})" for s in vus_signals]
        inh_info  = [f"{f['flag']}: {f['reasoning'][:120]}" for f in inh_flags]
        gap_info  = [f"{g['gap']}: {g['reasoning'][:120]}" for g in gap_flags]
        am_info   = [f"{s.get('gene','?')} {s['variant_id']}: {s['am_class']} ({s['am_pathogenicity']:.3f})" for s in am_signals]

        prompt = f"""Rare disease reanalysis signals for one patient:

Test: {request.prior_test.test_type} on {request.prior_test.test_date} ({years:.1f} years ago)
VUS Genes: {', '.join(request.prior_test.vus_genes) or 'None'}
Score: {result['reanalysis_score']}/100 ({result['confidence']}) — {result.get('urgency','')}

OMIM: {chr(10).join(['- '+s for s in omim_info]) or 'None'}
ClinVar: {chr(10).join(['- '+s for s in vus_info]) or 'None'}
New symptoms: {chr(10).join(['- '+s for s in symptom_names]) or 'None'}
Inheritance flags: {chr(10).join(['- '+s for s in inh_info]) or 'None'}
Analysis gaps: {chr(10).join(['- '+s for s in gap_info]) or 'None'}
AlphaMissense: {chr(10).join(['- '+s for s in am_info]) or 'Not provided'}

Return ONLY valid JSON:
{{
  "top_reasons": ["reason1","reason2","reason3"],
  "checklist": ["action1","action2","action3","action4"],
  "narrative": "2-3 sentences with specific gene/signal references",
  "recommended_reanalysis_type": "trio|CNV_reanalysis|variant_reinterpretation|expanded_panel|genome_upgrade|no_action"
}}"""

        msg = client.messages.create(
            model="claude-sonnet-4-20250514",
            max_tokens=800,
            system="You are a clinical genomics expert. Be specific and actionable. Return only valid JSON.",
            messages=[{"role": "user", "content": prompt}],
        )

        text = msg.content[0].text
        if "```json" in text:
            text = text.split("```json")[1].split("```")[0]
        elif "```" in text:
            text = text.split("```")[1].split("```")[0]

        return json.loads(text.strip())

    except Exception as e:
        print(f"Claude API error: {e}")
        return _fallback_recommendation(result, years, new_symptoms, omim_signals, inh_flags, gap_flags)


def _fallback_recommendation(result, years, new_symptoms, omim_signals, inh_flags, gap_flags) -> Dict:
    reasons, checklist, rtype = [], [], "variant_reinterpretation"

    if omim_signals:
        new = [s for s in omim_signals if s.get("is_new_after_test")]
        if new:
            reasons.append(f"{new[0]['gene']}→{new[0]['disease']} link established after test date")
            checklist.append(f"Re-evaluate {new[0]['gene']} VUS against {new[0]['disease']} criteria")

    for flag in (inh_flags or [])[:2]:
        if flag["flag"] == "AR_SINGLE_HIT":
            reasons.append(f"AR gene {flag['gene']}: second hit likely missed"); rtype = "CNV_reanalysis"
            checklist.append(f"CNV calling for {flag['gene']}")
        elif flag["flag"] == "AD_DE_NOVO_MISSED":
            reasons.append(f"De novo in {flag['gene']} unconfirmable from singleton"); rtype = "trio"
            checklist.append("Order trio sequencing")
        elif flag["flag"] == "XLINKED_FEMALE":
            reasons.append(f"X-linked {flag['gene']} variant in female — check XLD")
            checklist.append(f"Review {flag['gene']} variant in X-linked dominant context")

    for gap in (gap_flags or [])[:2]:
        if gap["gap"] == "CNV_CALLING_NOT_PERFORMED":
            reasons.append("CNV calling not done (10.8% of diagnoses require it)")
            checklist.append("Retrospective CNV calling on existing BAM")
            if rtype == "variant_reinterpretation": rtype = "CNV_reanalysis"
        elif gap["gap"] == "GENE_PANEL_SCOPE_LIMITED":
            reasons.append("Panel scope excludes newly characterised genes"); rtype = "expanded_panel"
            checklist.append("Upgrade to whole exome sequencing")

    if len(reasons) < 3:
        reasons.append(f"{years:.1f} years of knowledge accumulation since test")
        checklist.append("Review VUS in updated ClinVar and OMIM")
    if len(checklist) < 4:
        checklist += ["Discuss with rare disease specialist", "Update HPO phenotype profile"]

    s = result.get("reanalysis_score", 0)
    return {
        "top_reasons": reasons[:3],
        "checklist":   checklist[:4],
        "narrative": (
            f"Score {s}/100 — "
            f"{'urgent' if s>=80 else 'recommended' if s>=65 else 'consider'} reanalysis. "
            f"{reasons[0] if reasons else ''}. Recommended: {rtype.replace('_',' ')}."
        ),
        "recommended_reanalysis_type": rtype,
    }