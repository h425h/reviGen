"""
reviGen Portal API — main_portal.py
=====================================

Full portal backend supporting both Patient and Doctor portals.

Endpoints:
  Health
    GET  /api/health

  Patient data
    POST /api/patient/save-symptoms       — save patient symptom form
    POST /api/patient/save-genomic        — save genomic data + trigger reanalysis
    GET  /api/patient/{patient_id}/summary — full patient record

  NLP + Reanalysis
    POST /api/extract-hpo                 — free text → HPO terms
    POST /api/reanalysis-score            — structured HPO → 7-signal score
    POST /api/extract-and-score           — free text → HPO → 7-signal score (full NLP pipeline)

  Referrals
    POST /api/referrals/create            — doctor creates referral, pushed to patient
    GET  /api/referrals/{patient_id}      — list referrals for patient
    GET  /api/referrals/{referral_id}/download — download referral PDF (stub)

  Chat
    POST /api/chat/send                   — send message (patient or doctor)
    GET  /api/chat/{patient_id}           — get full conversation thread

  Diagnosis
    POST /api/diagnosis/publish           — doctor publishes diagnosis to patient portal
    GET  /api/diagnosis/{patient_id}      — get published diagnosis

  Gene ranking (internal, used by reanalysis)
    POST /api/rank-genes                  — HPO terms → ranked candidate genes

Storage: In-memory dict (demo). Swap for SQLite / PostgreSQL in production.
"""

import os
import json
import uuid
from contextlib import asynccontextmanager
from datetime import datetime, timezone
from typing import List, Optional, Dict, Any

import anthropic
from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field

# ── Internal signal modules (same as existing pipeline) ──────────────────────
# These imports assume you run from the project root and signals/ is a package.
# If running standalone, the aggregator will gracefully handle missing modules.

try:
    # Relative imports — works when running as `uvicorn backend.main_portal:app`
    from .aggregator import compute_reanalysis_score
    from .signals.nlp_extractor import extract_hpo_from_notes
    from .signals.gene_ranker import rank_genes_by_phenotype
    SIGNALS_AVAILABLE = True
except ImportError:
    try:
        # Fallback: absolute imports — works when running directly from backend/
        from aggregator import compute_reanalysis_score
        from signals.nlp_extractor import extract_hpo_from_notes
        from signals.gene_ranker import rank_genes_by_phenotype
        SIGNALS_AVAILABLE = True
    except ImportError:
        SIGNALS_AVAILABLE = False
        print("[WARN] Signal modules not found — running in demo/stub mode.")


# ═══════════════════════════════════════════════════════════════════════════════
#  IN-MEMORY STORE  (replace with real DB in production)
# ═══════════════════════════════════════════════════════════════════════════════

class Store:
    """Simple in-memory store for demo. Thread-safe enough for single-process hackathon use."""

    def __init__(self):
        self.patients: Dict[str, Dict]     = {}
        self.messages: Dict[str, List]     = {}   # patient_id → list of messages
        self.referrals: Dict[str, List]    = {}   # patient_id → list of referrals
        self.diagnoses: Dict[str, Dict]    = {}   # patient_id → diagnosis dict
        self.reanalysis: Dict[str, Dict]   = {}   # patient_id → latest reanalysis result

        # Seed demo patients so the portals work immediately on first load
        self._seed_demo_data()

    def _seed_demo_data(self):
        """Pre-populate three demo patients matching the UI drop-down."""

        # ── Demo Patient 1: FOXG1 ─────────────────────────────────────────────
        self.patients["P047"] = {
            "patient_id":  "P047",
            "display":     "Patient #0047 — Female, age 7",
            "demographics": {
                "date_of_birth":    "2019-01-15",
                "sex":              "Female",
                "country_of_origin":"United Arab Emirates",
                "consanguinity":    "No",
            },
            "symptoms": {
                "past_symptoms":    "Feeding difficulties as an infant, hypotonia at birth, motor delay noticed at 6 months.",
                "current_symptoms": ["Developmental delay","Absent / limited speech","Hypotonia (low muscle tone)","Autistic features","Seizures","Regression"],
                "current_free":     "7-year-old with global developmental delay, absent speech, apraxia, autistic behaviour, choreiform movements, and dysarthria. Progressive course since 18 months.",
                "recent_changes":   "Started having absence seizures since March. New regression in communication.",
                "recent_onset":     "2026-03-01",
                "frequency":        "Daily",
                "family_history":   "No known family history of neurological conditions. Parents non-consanguineous.",
                "prior_diagnoses":  "Unspecified developmental delay, Epilepsy NOS",
            },
            "genomic": {
                "test_type":    "Gene Panel",
                "test_date":    "2019-06-01",
                "result":       "Variant of Uncertain Significance (VUS)",
                "vus_genes":    ["FOXG1"],
                "cnv_called":   "No",
                "analysis_type":"Singleton",
                "lab":          "GeneDx",
                "notes":        "HP:0001263, HP:0001252, HP:0011968 — FOXG1 c.460C>T p.Arg154Cys classified VUS",
            },
            "hpo_original": ["HP:0001252","HP:0001270","HP:0011968"],
            "hpo_current":  ["HP:0001263","HP:0002186","HP:0000252","HP:0001252","HP:0002072","HP:0000729","HP:0001344","HP:0001270"],
        }

        # ── Demo Patient 2: STXBP1 ───────────────────────────────────────────
        self.patients["P031"] = {
            "patient_id":  "P031",
            "display":     "Patient #0031 — Male, age 12",
            "demographics": {
                "date_of_birth":    "2013-08-22",
                "sex":              "Male",
                "country_of_origin":"Egypt",
                "consanguinity":    "Yes — first cousins",
            },
            "symptoms": {
                "past_symptoms":    "Neonatal seizures, hypotonia from birth.",
                "current_symptoms": ["Seizures","Hypotonia (low muscle tone)","Developmental delay","Absent / limited speech"],
                "current_free":     "12-year-old with refractory epilepsy, intellectual disability, absent speech, and hypotonia. Seizures are daily despite multiple AEDs.",
                "recent_changes":   "Seizure frequency increased. Now having tonic-clonic clusters.",
                "recent_onset":     "2025-11-01",
                "frequency":        "Daily",
                "family_history":   "Parents are first cousins. Older sibling deceased — cause unknown.",
                "prior_diagnoses":  "Ohtahara syndrome NOS, infantile spasms",
            },
            "genomic": {
                "test_type":    "Whole Exome Sequencing (WES)",
                "test_date":    "2018-03-15",
                "result":       "Variant of Uncertain Significance (VUS)",
                "vus_genes":    ["STXBP1","KCNQ2"],
                "cnv_called":   "No",
                "analysis_type":"Singleton",
                "lab":          "Centogene",
                "notes":        "STXBP1 p.Val84Asp VUS, KCNQ2 p.Arg144Cys VUS — both of uncertain significance",
            },
            "hpo_original": ["HP:0001250","HP:0001252","HP:0001344"],
            "hpo_current":  ["HP:0001250","HP:0001252","HP:0001344","HP:0001263","HP:0001249","HP:0002376"],
        }

        # ── Demo Patient 3: PTHR1/Nina case ──────────────────────────────────
        self.patients["P058"] = {
            "patient_id":  "P058",
            "display":     "Patient #0058 — Female, age 4",
            "demographics": {
                "date_of_birth":    "2021-05-10",
                "sex":              "Female",
                "country_of_origin":"Canada",
                "consanguinity":    "No",
            },
            "symptoms": {
                "past_symptoms":    "Short stature noted at 1 year. Rickets diagnosed at 18 months.",
                "current_symptoms": ["Short stature","Skeletal abnormalities","Feeding difficulties"],
                "current_free":     "4-year-old with severe rickets, hypercalcemia, nephrocalcinosis, and failure to thrive. Parathyroid hormone levels markedly elevated.",
                "recent_changes":   "Hypercalcemia worsened. Admitted for IV hydration.",
                "recent_onset":     "2025-12-01",
                "frequency":        "Constant",
                "family_history":   "No family history. De novo suspected.",
                "prior_diagnoses":  "Vitamin D-resistant rickets",
            },
            "genomic": {
                "test_type":    "Whole Exome Sequencing (WES)",
                "test_date":    "2023-02-01",
                "result":       "Negative / No variant found",
                "vus_genes":    [],
                "cnv_called":   "Yes",
                "analysis_type":"Trio",
                "lab":          "Ambry Genetics",
                "notes":        "No pathogenic or likely pathogenic variants. PTHR1 not included in analysis panel at time of testing.",
            },
            "hpo_original": ["HP:0001510","HP:0002748","HP:0003165"],
            "hpo_current":  ["HP:0001510","HP:0002748","HP:0003165","HP:0003072","HP:0000121","HP:0001508"],
        }

        # Seed referrals for P047
        self.referrals["P047"] = [
            {
                "referral_id":    "REF-001",
                "patient_id":     "P047",
                "specialist_name":"Dr. Sarah Okonkwo",
                "specialty":      "Paediatric Neurology",
                "hospital":       "Boston Children's Hospital",
                "urgency":        "Urgent (within 1 week)",
                "reason":         "Referral for evaluation of progressive neurodevelopmental regression, seizure management, and FOXG1 gene variant review. Please assess for FOXG1 syndrome phenotype.",
                "issued_by":      "Dr. Marcus Reinholt",
                "issued_at":      "2026-02-14T09:00:00Z",
                "status":         "issued",
            },
            {
                "referral_id":    "REF-002",
                "patient_id":     "P047",
                "specialist_name":"Dr. Marcus Reinholt",
                "specialty":      "Clinical Genetics",
                "hospital":       "MGH Genetics",
                "urgency":        "Urgent (within 1 week)",
                "reason":         "Referral for trio whole exome reanalysis in context of new phenotypic features since 2019 gene panel. Prior CNV analysis not performed — recommend retrospective calling on archived BAM.",
                "issued_by":      "Dr. Marcus Reinholt",
                "issued_at":      "2026-02-20T11:30:00Z",
                "status":         "issued",
            },
        ]

        # Seed chat for P047
        self.messages["P047"] = [
            {
                "message_id": "MSG-001",
                "patient_id": "P047",
                "sender":     "doctor",
                "sender_name":"Dr. Reinholt",
                "content":    "Hello! I've reviewed your updated symptom profile. A few questions — have the hand stereotypies been present since birth, or did they develop after the regression at 18 months?",
                "sent_at":    "2026-02-18T10:42:00Z",
            },
            {
                "message_id": "MSG-002",
                "patient_id": "P047",
                "sender":     "patient",
                "sender_name":"Patient / Family",
                "content":    "They started around 20 months — after the regression. Before that she had normal purposeful hand movements.",
                "sent_at":    "2026-02-18T14:15:00Z",
            },
            {
                "message_id": "MSG-003",
                "patient_id": "P047",
                "sender":     "doctor",
                "sender_name":"Dr. Reinholt",
                "content":    "Thank you, that's helpful. That timeline is consistent with what we're seeing in the reanalysis. I've sent referrals to Dr. Okonkwo and the genetics team. Please check the Referrals tab.",
                "sent_at":    "2026-02-19T09:05:00Z",
            },
        ]

        # Seed published diagnosis for P047
        self.diagnoses["P047"] = {
            "patient_id":          "P047",
            "suspected_diagnosis": "FOXG1 Syndrome (OMIM:613454)",
            "status":              "Suspected — awaiting confirmation",
            "next_step":           "Trio whole exome reanalysis",
            "reanalysis_score":    84,
            "similarity_score":    0.89,
            "similarity_method":   "Resnik IC-weighted",
            "signals_active":      ["OMIM HIGH ⚡", "X-linked female flag", "CNV not called", "Phenotypic drift +5 new symptoms"],
            "clinical_summary":    (
                "Based on the clinical profile, phenotypic trajectory (regression at 18 months, stereotypies, "
                "absent speech, seizures), and the FOXG1 VUS identified in the 2019 gene panel, FOXG1 syndrome "
                "is the leading candidate. The Resnik similarity score between the patient's HPO profile and the "
                "canonical FOXG1 syndrome profile is 0.89 — very high phenotypic specificity. "
                "Recommending urgent trio reanalysis with CNV calling."
            ),
            "doctor_public_notes": (
                "Referrals to paediatric neurology and clinical genetics have been issued. "
                "Please bring all prior test reports to upcoming appointments. "
                "We are here to support you through this process."
            ),
            "issued_by":   "Dr. Marcus Reinholt",
            "issued_at":   "2026-02-24T14:00:00Z",
            "published":   True,
        }

db = Store()


# ═══════════════════════════════════════════════════════════════════════════════
#  PYDANTIC MODELS
# ═══════════════════════════════════════════════════════════════════════════════

class DemographicsIn(BaseModel):
    date_of_birth:     Optional[str] = None
    sex:               Optional[str] = None
    country_of_origin: Optional[str] = None
    consanguinity:     Optional[str] = None

class SymptomsIn(BaseModel):
    past_symptoms:    Optional[str]       = None
    current_symptoms: Optional[List[str]] = []
    current_free:     Optional[str]       = None
    recent_changes:   Optional[str]       = None
    recent_onset:     Optional[str]       = None
    frequency:        Optional[str]       = None
    family_history:   Optional[str]       = None
    prior_diagnoses:  Optional[str]       = None

class SaveSymptomsRequest(BaseModel):
    patient_id:   str
    demographics: Optional[DemographicsIn] = None
    symptoms:     SymptomsIn

class GenomicIn(BaseModel):
    test_type:     Optional[str]       = None
    test_date:     Optional[str]       = None
    result:        Optional[str]       = None
    vus_genes:     Optional[List[str]] = []
    cnv_called:    Optional[str]       = None
    analysis_type: Optional[str]       = None
    lab:           Optional[str]       = None
    notes:         Optional[str]       = None

class SaveGenomicRequest(BaseModel):
    patient_id: str
    genomic:    GenomicIn

class ExtractHpoRequest(BaseModel):
    clinical_text: str
    context:       Optional[str] = None  # 'original' | 'current'

class ExtractAndScoreRequest(BaseModel):
    patient_id:              Optional[str] = None
    original_clinical_text:  str
    current_clinical_text:   str
    test_type:               Optional[str]       = "Gene Panel"
    test_date:               Optional[str]       = "2019-01-01"
    vus_genes:               Optional[List[str]] = []
    analysis_type:           Optional[str]       = "singleton"
    patient_sex:             Optional[str]       = None
    cnv_calling_performed:   Optional[bool]      = None
    consanguineous:          Optional[bool]      = None

class ReanalysisRequest(BaseModel):
    patient_id:          Optional[str]       = None
    original_hpo_terms:  List[str]
    current_hpo_terms:   List[str]
    test_type:           str                 = "Gene Panel"
    test_date:           str                 = "2019-01-01"
    vus_genes:           List[str]           = []
    analysis_type:       Optional[str]       = "singleton"
    patient_sex:         Optional[str]       = None
    cnv_calling_performed: Optional[bool]    = None
    consanguineous:      Optional[bool]      = None

class RankGenesRequest(BaseModel):
    hpo_terms: List[str]
    top_k:     int = 10

class CreateReferralRequest(BaseModel):
    patient_id:       str
    specialist_name:  str
    specialty:        str
    hospital:         str
    urgency:          str
    reason:           str
    issued_by:        str = "Dr. (reviGen)"

class SendMessageRequest(BaseModel):
    patient_id:   str
    sender:       str        # 'patient' | 'doctor'
    sender_name:  str
    content:      str

class PublishDiagnosisRequest(BaseModel):
    patient_id:            str
    suspected_diagnosis:   str
    status:                str
    next_step:             str
    clinical_summary:      str
    doctor_public_notes:   Optional[str] = None
    issued_by:             str = "Dr. (reviGen)"


# ═══════════════════════════════════════════════════════════════════════════════
#  APP SETUP
# ═══════════════════════════════════════════════════════════════════════════════

@asynccontextmanager
async def lifespan(app: FastAPI):
    print("=" * 60)
    print("reviGen Portal API — Starting")
    print(f"  ANTHROPIC_API_KEY : {'SET' if os.environ.get('ANTHROPIC_API_KEY') else 'not set (NLP fallback active)'}")
    print(f"  OMIM_API_KEY      : {'SET' if os.environ.get('OMIM_API_KEY') else 'not set (OMIM fallback active)'}")
    print(f"  Signals available : {SIGNALS_AVAILABLE}")
    print(f"  Demo patients     : {list(db.patients.keys())}")
    print("=" * 60)
    yield
    print("reviGen Portal API — Shutdown")


app = FastAPI(
    title="reviGen Portal API",
    description="Full patient + doctor portal backend for reviGen rare disease reanalysis tool",
    version="4.0.0",
    lifespan=lifespan,
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


# ═══════════════════════════════════════════════════════════════════════════════
#  HEALTH
# ═══════════════════════════════════════════════════════════════════════════════

@app.get("/api/health")
async def health():
    return {
        "status":           "healthy",
        "signals_available": SIGNALS_AVAILABLE,
        "demo_patients":    list(db.patients.keys()),
        "anthropic_key":    bool(os.environ.get("ANTHROPIC_API_KEY")),
        "omim_key":         bool(os.environ.get("OMIM_API_KEY")),
    }


# ═══════════════════════════════════════════════════════════════════════════════
#  PATIENT LIST (for doctor portal drop-down)
# ═══════════════════════════════════════════════════════════════════════════════

@app.get("/api/patients")
async def list_patients():
    return {
        "patients": [
            {"patient_id": pid, "display": p.get("display", pid)}
            for pid, p in db.patients.items()
        ]
    }


# ═══════════════════════════════════════════════════════════════════════════════
#  PATIENT DATA
# ═══════════════════════════════════════════════════════════════════════════════

@app.get("/api/patient/{patient_id}/summary")
async def get_patient_summary(patient_id: str):
    """Full patient record — used to populate doctor's side panel."""
    p = db.patients.get(patient_id)
    if not p:
        raise HTTPException(404, f"Patient {patient_id} not found")
    return {
        **p,
        "has_referrals": bool(db.referrals.get(patient_id)),
        "has_diagnosis":  bool(db.diagnoses.get(patient_id)),
        "message_count":  len(db.messages.get(patient_id, [])),
        "reanalysis":     db.reanalysis.get(patient_id),
    }


@app.post("/api/patient/save-symptoms")
async def save_symptoms(req: SaveSymptomsRequest):
    """Save patient symptom form. Creates patient record if not exists."""
    if req.patient_id not in db.patients:
        db.patients[req.patient_id] = {
            "patient_id": req.patient_id,
            "display":    f"Patient {req.patient_id}",
        }
    p = db.patients[req.patient_id]
    if req.demographics:
        p["demographics"] = req.demographics.dict()
    p["symptoms"] = req.symptoms.dict()
    return {"ok": True, "patient_id": req.patient_id}


@app.post("/api/patient/save-genomic")
async def save_genomic(req: SaveGenomicRequest):
    """
    Save patient genomic data.
    If prior VUS genes and test date are present, automatically triggers
    the reanalysis pipeline using existing HPO terms (if available).
    """
    if req.patient_id not in db.patients:
        db.patients[req.patient_id] = {"patient_id": req.patient_id}

    p = db.patients[req.patient_id]
    p["genomic"] = req.genomic.dict()

    result_summary = {"ok": True, "patient_id": req.patient_id, "reanalysis_triggered": False}

    # Auto-trigger reanalysis if we have HPO terms + VUS genes + test date
    if (req.genomic.vus_genes and req.genomic.test_date
            and p.get("hpo_original") and p.get("hpo_current")):
        try:
            score_result = _run_reanalysis(
                vus_genes=req.genomic.vus_genes,
                test_date=req.genomic.test_date,
                original_hpo=p["hpo_original"],
                current_hpo=p["hpo_current"],
                test_type=req.genomic.test_type or "Gene Panel",
                analysis_type=req.genomic.analysis_type,
                cnv_calling_performed=(req.genomic.cnv_called == "Yes"),
                patient_sex=p.get("demographics", {}).get("sex"),
            )
            db.reanalysis[req.patient_id] = score_result
            result_summary["reanalysis_triggered"] = True
            result_summary["score"] = score_result.get("reanalysis_score")
        except Exception as e:
            result_summary["reanalysis_error"] = str(e)

    return result_summary


# ═══════════════════════════════════════════════════════════════════════════════
#  NLP + REANALYSIS
# ═══════════════════════════════════════════════════════════════════════════════

@app.post("/api/extract-hpo")
async def extract_hpo(req: ExtractHpoRequest):
    """
    Free-text clinical note → HPO terms.
    Uses Claude API (falls back to keyword matching if no API key).
    """
    if not req.clinical_text.strip():
        raise HTTPException(400, "clinical_text is required")
    result = _extract_hpo(req.clinical_text)
    return {
        "hpo_terms":  result["hpo_terms"],
        "term_count": result["term_count"],
        "method":     result["method"],
    }


@app.post("/api/extract-and-score")
async def extract_and_score(req: ExtractAndScoreRequest):
    """
    Full NLP pipeline in one call:
      original free text → HPO terms
      current free text  → HPO terms
      → 7-signal reanalysis score + Claude clinical recommendation

    This is the primary endpoint for the 'Clinical Notes (NLP)' tab.
    """
    api_key = os.environ.get("ANTHROPIC_API_KEY", "")

    orig_result = _extract_hpo(req.original_clinical_text)
    curr_result = _extract_hpo(req.current_clinical_text)

    original_hpo = orig_result["hpo_terms"]
    current_hpo  = curr_result["hpo_terms"]

    score_result = _run_reanalysis(
        vus_genes=req.vus_genes or [],
        test_date=req.test_date or "2019-01-01",
        original_hpo=original_hpo,
        current_hpo=current_hpo,
        test_type=req.test_type or "Gene Panel",
        analysis_type=req.analysis_type,
        cnv_calling_performed=req.cnv_calling_performed,
        patient_sex=req.patient_sex,
        consanguineous=req.consanguineous,
    )

    # Get Claude clinical recommendation
    score_result["clinical_recommendation"] = await _claude_recommendation(
        vus_genes=req.vus_genes or [],
        test_date=req.test_date or "2019-01-01",
        score_result=score_result,
    )

    # Store if patient_id given
    if req.patient_id and req.patient_id in db.patients:
        db.patients[req.patient_id]["hpo_original"] = original_hpo
        db.patients[req.patient_id]["hpo_current"]  = current_hpo
        db.reanalysis[req.patient_id] = score_result

    score_result["nlp_extraction"] = {
        "original_hpo_terms":  original_hpo,
        "current_hpo_terms":   current_hpo,
        "original_method":     orig_result["method"],
        "current_method":      curr_result["method"],
        "original_term_count": orig_result["term_count"],
        "current_term_count":  curr_result["term_count"],
    }

    return score_result


@app.post("/api/reanalysis-score")
async def reanalysis_score(req: ReanalysisRequest):
    """Structured HPO → 7-signal reanalysis score. Used by the structured input tab."""
    score_result = _run_reanalysis(
        vus_genes=req.vus_genes,
        test_date=req.test_date,
        original_hpo=req.original_hpo_terms,
        current_hpo=req.current_hpo_terms,
        test_type=req.test_type,
        analysis_type=req.analysis_type,
        cnv_calling_performed=req.cnv_calling_performed,
        patient_sex=req.patient_sex,
        consanguineous=req.consanguineous,
    )
    score_result["clinical_recommendation"] = await _claude_recommendation(
        vus_genes=req.vus_genes,
        test_date=req.test_date,
        score_result=score_result,
    )

    if req.patient_id and req.patient_id in db.patients:
        db.patients[req.patient_id]["hpo_original"] = req.original_hpo_terms
        db.patients[req.patient_id]["hpo_current"]  = req.current_hpo_terms
        db.reanalysis[req.patient_id] = score_result

    return score_result


@app.post("/api/rank-genes")
async def rank_genes(req: RankGenesRequest):
    """HPO terms → ranked candidate genes (disease-intermediate method)."""
    if not req.hpo_terms:
        raise HTTPException(400, "hpo_terms is required")
    if not SIGNALS_AVAILABLE:
        return {"ranked_genes": _stub_gene_ranks(req.hpo_terms), "gene_count": 5, "similarity_method": "stub"}
    ranked = rank_genes_by_phenotype(req.hpo_terms, top_k=req.top_k)
    return {
        "ranked_genes":      ranked,
        "gene_count":        len(ranked),
        "similarity_method": ranked[0]["similarity_method"] if ranked else "none",
        "input_hpo_count":   len(req.hpo_terms),
    }


# ═══════════════════════════════════════════════════════════════════════════════
#  REFERRALS
# ═══════════════════════════════════════════════════════════════════════════════

@app.get("/api/referrals/{patient_id}")
async def get_referrals(patient_id: str):
    refs = db.referrals.get(patient_id, [])
    return {"patient_id": patient_id, "referrals": refs, "count": len(refs)}


@app.post("/api/referrals/create")
async def create_referral(req: CreateReferralRequest):
    """Doctor creates a referral — it is immediately visible in patient portal."""
    if req.patient_id not in db.patients:
        raise HTTPException(404, f"Patient {req.patient_id} not found")

    referral = {
        "referral_id":    f"REF-{uuid.uuid4().hex[:6].upper()}",
        "patient_id":     req.patient_id,
        "specialist_name":req.specialist_name,
        "specialty":      req.specialty,
        "hospital":       req.hospital,
        "urgency":        req.urgency,
        "reason":         req.reason,
        "issued_by":      req.issued_by,
        "issued_at":      datetime.now(timezone.utc).isoformat(),
        "status":         "issued",
    }

    if req.patient_id not in db.referrals:
        db.referrals[req.patient_id] = []
    db.referrals[req.patient_id].append(referral)

    return {"ok": True, "referral": referral}


# ═══════════════════════════════════════════════════════════════════════════════
#  CHAT
# ═══════════════════════════════════════════════════════════════════════════════

@app.get("/api/chat/{patient_id}")
async def get_chat(patient_id: str):
    msgs = db.messages.get(patient_id, [])
    return {"patient_id": patient_id, "messages": msgs, "count": len(msgs)}


@app.post("/api/chat/send")
async def send_message(req: SendMessageRequest):
    """Send a message from patient or doctor. Both portals call the same endpoint."""
    if req.patient_id not in db.patients:
        raise HTTPException(404, f"Patient {req.patient_id} not found")
    if req.sender not in ("patient", "doctor"):
        raise HTTPException(400, "sender must be 'patient' or 'doctor'")
    if not req.content.strip():
        raise HTTPException(400, "content is required")

    msg = {
        "message_id":  f"MSG-{uuid.uuid4().hex[:6].upper()}",
        "patient_id":  req.patient_id,
        "sender":      req.sender,
        "sender_name": req.sender_name,
        "content":     req.content.strip(),
        "sent_at":     datetime.now(timezone.utc).isoformat(),
    }

    if req.patient_id not in db.messages:
        db.messages[req.patient_id] = []
    db.messages[req.patient_id].append(msg)

    return {"ok": True, "message": msg}


# ═══════════════════════════════════════════════════════════════════════════════
#  DIAGNOSIS
# ═══════════════════════════════════════════════════════════════════════════════

@app.get("/api/diagnosis/{patient_id}")
async def get_diagnosis(patient_id: str):
    diag = db.diagnoses.get(patient_id)
    if not diag:
        return {"patient_id": patient_id, "published": False, "diagnosis": None}
    return {"patient_id": patient_id, "published": diag.get("published", False), "diagnosis": diag}


@app.post("/api/diagnosis/publish")
async def publish_diagnosis(req: PublishDiagnosisRequest):
    """Doctor publishes diagnosis — immediately visible in patient's Diagnosis tab."""
    if req.patient_id not in db.patients:
        raise HTTPException(404, f"Patient {req.patient_id} not found")

    # Pull in latest reanalysis score if available
    reanalysis = db.reanalysis.get(req.patient_id, {})

    diag = {
        "patient_id":           req.patient_id,
        "suspected_diagnosis":  req.suspected_diagnosis,
        "status":               req.status,
        "next_step":            req.next_step,
        "clinical_summary":     req.clinical_summary,
        "doctor_public_notes":  req.doctor_public_notes or "",
        "reanalysis_score":     reanalysis.get("reanalysis_score", 0),
        "similarity_score":     reanalysis.get("signal_breakdown", {})
                                    .get("disease_match", {})
                                    .get("score", 0),
        "signals_active":       _extract_active_signals(reanalysis),
        "issued_by":            req.issued_by,
        "issued_at":            datetime.now(timezone.utc).isoformat(),
        "published":            True,
    }

    db.diagnoses[req.patient_id] = diag
    return {"ok": True, "diagnosis": diag}


# ═══════════════════════════════════════════════════════════════════════════════
#  INTERNAL HELPERS
# ═══════════════════════════════════════════════════════════════════════════════

def _extract_hpo(text: str) -> Dict:
    """
    Wrap NLP extractor with graceful fallback.
    Normalises the richer extract_hpo_from_notes output into the simple
    {hpo_terms, term_count, method} shape the rest of main_portal.py expects.
    """
    if SIGNALS_AVAILABLE:
        try:
            raw = extract_hpo_from_notes(
                clinical_text=text,
                api_key=os.environ.get("ANTHROPIC_API_KEY", ""),
            )
            # extract_hpo_from_notes returns a rich dict — normalise it
            hpo_terms = raw.get("hpo_terms", [])
            method    = raw.get("extraction_method", "claude_api")
            # Map internal method names to friendly labels
            method_label = "claude" if "claude" in method else "keyword"
            return {
                "hpo_terms":  hpo_terms,
                "term_count": len(hpo_terms),
                "method":     method_label,
                # Pass through rich fields for potential future use
                "hpo_with_names":  raw.get("hpo_with_names", []),
                "temporal_hpo":    raw.get("temporal_hpo", []),
                "excluded_terms":  raw.get("excluded_terms", []),
                "unmapped_terms":  raw.get("unmapped_terms", []),
            }
        except Exception as e:
            print(f"[NLP] Error: {e}")

    # Keyword fallback (inline, no signals dependency)
    return _keyword_extract(text)


def _run_reanalysis(
    vus_genes, test_date, original_hpo, current_hpo,
    test_type="Gene Panel", analysis_type=None,
    cnv_calling_performed=None, patient_sex=None, consanguineous=None,
) -> Dict:
    """Run 7-signal reanalysis pipeline or return a demo stub."""
    if SIGNALS_AVAILABLE:
        try:
            return compute_reanalysis_score(
                vus_genes=vus_genes,
                test_date=test_date,
                original_hpo_terms=original_hpo,
                current_hpo_terms=current_hpo,
                test_type=test_type,
                patient_sex=patient_sex,
                analysis_type=analysis_type,
                consanguineous=consanguineous,
                cnv_calling_performed=cnv_calling_performed,
                omim_api_key=os.environ.get("OMIM_API_KEY", ""),
            )
        except Exception as e:
            print(f"[Reanalysis] Error: {e}")

    return _demo_score_stub(vus_genes, test_date, original_hpo, current_hpo, cnv_calling_performed)


async def _claude_recommendation(vus_genes: List[str], test_date: str, score_result: Dict) -> Dict:
    """Generate Claude clinical recommendation. Falls back to rule-based if no API key."""
    api_key = os.environ.get("ANTHROPIC_API_KEY", "")
    sb      = score_result.get("signal_breakdown", {})
    score   = score_result.get("reanalysis_score", 0)
    urgency = score_result.get("urgency", "")

    if not api_key:
        return _rule_based_recommendation(score_result)

    try:
        client = anthropic.Anthropic(api_key=api_key)

        omim_signals = sb.get("omim_gene_disease", {}).get("signals", [])
        inh_flags    = sb.get("inheritance_pattern", {}).get("flags", [])
        gap_flags    = sb.get("analysis_method_gaps", {}).get("gaps", [])
        drift_detail = sb.get("phenotypic_drift", {}).get("details", {})
        new_symptoms = drift_detail.get("new_symptoms", drift_detail.get("new_hpo_terms", []))
        entropy      = sb.get("entropy_modifier", {}).get("entropy_value", 0)

        prompt = f"""Rare disease reanalysis signals:
Test: {test_date} | Genes: {', '.join(vus_genes) or 'None'} | Score: {score}/100 | {urgency}

OMIM new signals: {json.dumps([{'gene':s.get('gene'),'disease':s.get('disease'),'new':s.get('is_new_after_test')} for s in omim_signals[:3]])}
Inheritance flags: {json.dumps([f.get('flag') for f in inh_flags[:3]])}
Analysis gaps: {json.dumps([g.get('gap') for g in gap_flags[:3]])}
New HPO terms since test: {len(new_symptoms)}
Entropy (multi-VUS uncertainty): {entropy:.2f}

Return ONLY valid JSON:
{{
  "top_reasons": ["reason1","reason2","reason3"],
  "checklist": ["action1","action2","action3","action4"],
  "narrative": "2-3 sentences citing specific gene names and signals",
  "recommended_reanalysis_type": "trio|CNV_reanalysis|variant_reinterpretation|expanded_panel|genome_upgrade|no_action"
}}"""

        msg = client.messages.create(
            model="claude-sonnet-4-20250514",
            max_tokens=600,
            system="You are a clinical genomics expert. Be specific and actionable. Return only valid JSON, no markdown.",
            messages=[{"role": "user", "content": prompt}],
        )
        text = msg.content[0].text.strip()
        if "```" in text:
            text = text.split("```")[1].lstrip("json").strip()
            text = text.split("```")[0].strip()
        return json.loads(text)
    except Exception as e:
        print(f"[Claude] Recommendation error: {e}")
        return _rule_based_recommendation(score_result)


def _rule_based_recommendation(score_result: Dict) -> Dict:
    """Rule-based fallback when Claude API is unavailable."""
    sb      = score_result.get("signal_breakdown", {})
    score   = score_result.get("reanalysis_score", 0)
    reasons, checklist = [], []

    # OMIM signal
    omim = sb.get("omim_gene_disease", {}).get("signals", [])
    new  = [s for s in omim if s.get("is_new_after_test")]
    if new:
        g = new[0]
        reasons.append(f"{g['gene']} → {g['disease']} association established after original test")
        checklist.append(f"Re-evaluate {g['gene']} VUS against updated {g['disease']} criteria")

    # Inheritance flags
    for flag in sb.get("inheritance_pattern", {}).get("flags", [])[:2]:
        f = flag.get("flag", "")
        if f == "AR_SINGLE_HIT":
            reasons.append("AR gene — likely second hit (CNV/deep intronic) missed")
            checklist.append("Perform CNV calling on archived BAM")
        elif f == "AD_DE_NOVO_MISSED":
            reasons.append("De novo variant unconfirmable from singleton analysis")
            checklist.append("Order trio sequencing")
        elif f == "XLINKED_FEMALE":
            reasons.append("X-linked dominant variant in female — carrier classification likely incorrect")
            checklist.append("Review X-linked dominant penetrance in females")

    # Analysis gaps
    for gap in sb.get("analysis_method_gaps", {}).get("gaps", [])[:2]:
        g = gap.get("gap", "")
        if g == "CNV_CALLING_NOT_PERFORMED":
            reasons.append("CNV analysis not performed — diagnostic in ~11% of cases")
            checklist.append("Retrospective CNV calling on archived BAM")
        elif g == "GENE_PANEL_SCOPE_LIMITED":
            reasons.append("Panel scope excludes recently characterised genes")
            checklist.append("Upgrade to whole exome sequencing")

    # Entropy
    entropy = sb.get("entropy_modifier", {}).get("entropy_value", 0)
    if entropy > 0.6:
        reasons.append(f"High diagnostic uncertainty across {sb.get('entropy_modifier',{}).get('n_vus_genes',2)} VUS genes (entropy {entropy:.2f})")

    reasons   = (reasons + ["Knowledge base substantially updated since test"])[:3]
    checklist = (checklist + ["Review VUS in updated ClinVar and OMIM", "Consult rare disease specialist"])[:4]

    rtype = "trio" if score >= 80 else "CNV_reanalysis" if score >= 65 else "variant_reinterpretation" if score >= 40 else "no_action"

    return {
        "top_reasons":                 reasons,
        "checklist":                   checklist,
        "narrative":                   f"Score {score}/100. {reasons[0]}. Recommended: {rtype.replace('_', ' ')}.",
        "recommended_reanalysis_type": rtype,
    }


def _extract_active_signals(reanalysis: Dict) -> List[str]:
    """Extract human-readable signal labels from a reanalysis result."""
    signals = []
    sb = reanalysis.get("signal_breakdown", {})

    if sb.get("omim_gene_disease", {}).get("score", 0) > 0.3:
        signals.append("OMIM HIGH ⚡")
    for flag in sb.get("inheritance_pattern", {}).get("flags", [])[:1]:
        signals.append(f"Inheritance: {flag.get('flag','')}")
    for gap in sb.get("analysis_method_gaps", {}).get("gaps", [])[:1]:
        signals.append(f"Gap: {gap.get('gap','').replace('_',' ')}")
    drift = sb.get("phenotypic_drift", {}).get("details", {})
    new_syms = drift.get("new_symptoms", drift.get("new_hpo_terms", []))
    if new_syms:
        signals.append(f"Drift: {len(new_syms)} new symptoms")
    if sb.get("entropy_modifier", {}).get("entropy_value", 0) > 0.5:
        signals.append("High entropy — multi-VUS uncertainty")

    return signals or ["7-signal analysis complete"]


# ── Keyword HPO extraction fallback ───────────────────────────────────────────

HPO_KEYWORDS = {
    "hypotonia": "HP:0001252", "floppy": "HP:0001252", "low muscle tone": "HP:0001252",
    "hypertonia": "HP:0001276", "spasticity": "HP:0001257",
    "seizure": "HP:0001250", "epilepsy": "HP:0001250", "convulsion": "HP:0001250",
    "developmental delay": "HP:0001263", "global developmental delay": "HP:0001263",
    "intellectual disability": "HP:0001249", "intellectual disability": "HP:0001249",
    "absent speech": "HP:0001344", "no speech": "HP:0001344", "nonverbal": "HP:0001344",
    "speech delay": "HP:0000750", "language delay": "HP:0000750",
    "apraxia": "HP:0002186", "dysarthria": "HP:0001260",
    "ataxia": "HP:0001251", "unsteady gait": "HP:0001251",
    "dystonia": "HP:0001332", "stereotypies": "HP:0002063",
    "hand wringing": "HP:0002063", "hand stereotypies": "HP:0002063",
    "choreiform": "HP:0002072", "chorea": "HP:0002072",
    "autism": "HP:0000729", "autistic": "HP:0000729",
    "regression": "HP:0002376", "regress": "HP:0002376",
    "microcephaly": "HP:0000252", "small head": "HP:0000252",
    "macrocephaly": "HP:0000256", "large head": "HP:0000256",
    "feeding difficulties": "HP:0011968", "poor feeding": "HP:0011968",
    "motor delay": "HP:0001270", "delayed motor": "HP:0001270",
    "short stature": "HP:0004322", "growth retardation": "HP:0001510",
    "failure to thrive": "HP:0001508",
    "rickets": "HP:0002748", "hypercalcemia": "HP:0003072",
    "nephrocalcinosis": "HP:0000121",
    "hearing loss": "HP:0000365", "visual impairment": "HP:0000505",
    "skeletal": "HP:0000924", "joint": "HP:0001376",
    "neonatal": "HP:0003623", "infantile": "HP:0003593",
}

def _keyword_extract(text: str) -> Dict:
    text_lower = text.lower()
    found = {}
    for kw, hpo_id in HPO_KEYWORDS.items():
        if kw in text_lower:
            found[hpo_id] = True
    terms = list(found.keys())
    return {"hpo_terms": terms, "term_count": len(terms), "method": "keyword"}


# ── Demo score stub (when signals/ not installed) ─────────────────────────────

def _demo_score_stub(vus_genes, test_date, original_hpo, current_hpo, cnv_done) -> Dict:
    """Return a plausible demo score matching v3 aggregator output when signals unavailable."""
    new_syms = [h for h in (current_hpo or []) if h not in set(original_hpo or [])]
    score = 45
    if vus_genes:         score += 20
    if not cnv_done:      score += 12
    if len(new_syms) > 2: score += 10
    score = min(score, 97)

    confidence = "URGENT" if score >= 80 else "HIGH" if score >= 65 else "MEDIUM"
    urgency    = ("Immediate escalation — reanalysis within days" if score >= 80
                  else "Reanalysis recommended — schedule within weeks" if score >= 65
                  else "Consider reanalysis — review at next clinic visit")

    return {
        "reanalysis_score": score,
        "raw_score":        round(score / 100, 3),
        "base_score":       round(score / 100, 3),
        "entropy_boost":    0.0,
        "confidence":       confidence,
        "urgency":          urgency,
        "overall_signal_strength": "HIGH" if score >= 80 else "MEDIUM",
        "signal_breakdown": {
            "omim_gene_disease": {
                "score": 0.8 if vus_genes else 0.1,
                "signals": [{"gene": g, "disease": "associated disorder", "is_new_after_test": True, "signal_strength": "HIGH"} for g in (vus_genes or [])],
                "signal_strength": "HIGH" if vus_genes else "LOW",
            },
            "vus_reclassification": {"score": 0.1, "signals": [], "signal_strength": "LOW"},
            "phenotypic_drift": {
                "score": 0.6 if new_syms else 0.1,
                "details": {"new_hpo_terms": new_syms, "new_symptoms": [{"hpo_id": h} for h in new_syms]},
                "signal_strength": "HIGH" if len(new_syms) > 3 else "MEDIUM",
            },
            "disease_match": {
                "score": 0.89, "signal_strength": "HIGH",
                "best_matching_disease": vus_genes[0] + " syndrome" if vus_genes else None,
                "similarity_method": "jaccard",
            },
            "inheritance_pattern": {"score": 0.0, "flags": [], "flag_count": 0, "signal_strength": "LOW"},
            "analysis_method_gaps": {
                "score": 0.7 if not cnv_done else 0.0,
                "gaps": [] if cnv_done else [{"gap": "CNV_CALLING_NOT_PERFORMED", "reasoning": "CNV calling was not performed on original data"}],
                "gap_count": 0 if cnv_done else 1,
                "signal_strength": "HIGH" if not cnv_done else "LOW",
            },
            "alphamissense": {"score": 0.0, "signals": [], "signal_strength": "LOW"},
            "time_since_test": {"score": 0.3, "details": {"years_since_test": 6}, "signal_strength": "MEDIUM"},
            "entropy_modifier": {"entropy_value": 0.0, "entropy_boost": 0.0, "n_vus_genes": len(vus_genes or []), "signal_strength": "LOW"},
        },
    }


def _stub_gene_ranks(hpo_terms: List[str]) -> List[Dict]:
    """Return top gene candidates for demo mode."""
    return [
        {"rank": 1, "gene": "FOXG1",  "disease_id": "OMIM:613454", "disease_name": "FOXG1 syndrome",              "score": 0.89, "similarity_method": "stub"},
        {"rank": 2, "gene": "CDKL5",  "disease_id": "OMIM:300672", "disease_name": "CDKL5 deficiency disorder",   "score": 0.72, "similarity_method": "stub"},
        {"rank": 3, "gene": "MECP2",  "disease_id": "OMIM:312750", "disease_name": "Rett syndrome",               "score": 0.65, "similarity_method": "stub"},
        {"rank": 4, "gene": "STXBP1", "disease_id": "OMIM:308350", "disease_name": "Developmental epilepsy DEE4", "score": 0.58, "similarity_method": "stub"},
        {"rank": 5, "gene": "SCN2A",  "disease_id": "OMIM:613721", "disease_name": "DEE11",                       "score": 0.51, "similarity_method": "stub"},
    ]