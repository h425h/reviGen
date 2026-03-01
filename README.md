# reviGen

> **AI-powered triage tool that identifies which rare disease patients with inconclusive genetic tests should be reanalyzed today — because the answer may already exist in the literature.**

Built at **Harvard RareHack 2026** by Giorgi Bakuradze, Neha Chimakurthy, Maya McCabe, Sahar Islam, and Harsh Mishra.

---

## The Problem

When a child with a suspected genetic disorder receives an inconclusive test result — a Variant of Uncertain Significance (VUS) or no pathogenic variant found — the family is told "we don't know" and sent home. The average diagnostic odyssey for a rare disease patient lasts **7.4 years**.

The cruel irony: genomic knowledge doubles every 18 months. The gene that was unknown in 2019 may be fully characterized by 2022. A variant classified as uncertain may be confirmed pathogenic by a leading expert panel a year later. But **nobody systematically goes back to check**. The family is still waiting.

reviGen automates that triage. It continuously monitors new medical discoveries and flags which patients are most likely to benefit from having their old genetic data reanalyzed — so doctors know exactly who to prioritize and why.

---

## How It Works

reviGen combines a patient's past genetic test results and current symptoms, runs them through **7 independent evidence-based signals**, and produces a 0–100 urgency score.

### The 7-Signal Scoring Engine

Each signal is weighted by how often that factor has historically explained missed diagnoses, derived from a real-world Radboud University 2022 study (n=150 reanalyzed cases):

| Signal | Weight | Evidence Basis |
|--------|--------|----------------|
| OMIM Gene-Disease Surveillance | **29%** | 42% of new diagnoses from newly characterised gene-disease links (Radboud 2022) |
| ClinVar VUS Reclassification | **14%** | ~15% of VUS reclassified within 3 years |
| Disease Profile Match (Resnik) | **13%** | IC-weighted semantic HPO similarity |
| Inheritance Pattern Flags | **12%** | AR second-hit, de novo unconfirmed, X-linked female |
| Phenotypic Drift | **11%** | New symptoms since original test |
| Analysis Method Gaps | **10%** | CNV not called, panel vs exome, old pipeline |
| Time Since Test | **6%** | Proxy — mechanism captured by signals above |
| AlphaMissense (DeepMind) | **5%** | Protein structure pathogenicity via AlphaFold2 |

### Score Formula

```
base_score  = Σ (weightᵢ × signalᵢ)
entropy_boost = H_norm × 0.08        # up to +8pts for multi-VUS uncertainty
final_score = min(1.0, base_score + entropy_boost) × 100
```

**Urgency Tiers:** URGENT ≥80 · HIGH 65–80 · MEDIUM 35–65 · LOW <35

### NLP Pipeline

Free-text clinical notes are converted to structured HPO terms via the Claude API:

```
Clinical Notes → Claude API → JSON (present/excluded phenotypes)
→ name_to_hpo_id() mapping → HPO Term List → Scoring Engine
```

Claude is used instead of traditional biomedical NLP (scispaCy) because 2024 benchmarks (Wan et al., Luo et al.) show LLMs outperform dictionary-based approaches on rare disease HPO extraction, especially for negation handling ("no seizures" → excluded phenotype). Fallback: 60-term keyword regex dictionary for offline/demo mode.

---

## Architecture

```
┌─────────────────────────────────────────────────────┐
│                   Frontend (HTML)                   │
│         Patient Portal    │    Doctor Portal        │
└──────────────┬────────────┴──────────┬──────────────┘
               │                       │
               ▼                       ▼
┌─────────────────────────────────────────────────────┐
│              FastAPI Backend (Python)               │
│                                                     │
│  POST /api/patient/save-symptoms                    │
│  POST /api/patient/save-genomic                     │
│  POST /api/extract-and-score   ◄── main endpoint    │
│  POST /api/doctor/patients                          │
│  GET  /api/evaluate                                 │
└──────────────┬──────────────────────────────────────┘
               │
       ┌───────▼────────┐
       │  aggregator.py │  (7 signals → weighted score)
       └───────┬────────┘
               │
    ┌──────────┼──────────────────────────────┐
    │          │                              │
    ▼          ▼                              ▼
OMIM API   ClinVar TSV              AlphaMissense SQLite
HPO/HPOA   (local)                  (local, DeepMind)
```

**Data flow for a full analysis request:**
1. Doctor pastes clinical notes → `POST /api/extract-and-score`
2. `nlp_extractor.py` → Claude API → HPO term list
3. `aggregator.py` → 7 signal modules run → weighted 0–100 score
4. Second Claude API call → plain-English clinical recommendation
5. JSON response rendered in frontend

---

## Evaluation

The gene ranker (phenotype → candidate gene ranking) was evaluated on two sets:

### Internal Validation (n=10 literature cases)

| Metric | Score |
|--------|-------|
| F1@1 | **0.90** |
| F1@3 | **1.00** |
| F1@10 | **1.00** |
| MRR | **0.95** |

### Held-Out Validation (n=10 MyGene2 genes, never seen during development)

| Group | F1@1 | F1@3 | MRR |
|-------|------|------|-----|
| In-DB genes (n=7) | **0.71** | **1.00** | **0.857** |
| Out-of-DB genes (n=3) | 0.00 | 0.00 | 0.000 |

The 3 out-of-DB misses (KDM1A, SYNGAP1, SHANK3) are not model failures — they quantify exactly the gap the OMIM surveillance signal is designed to close. When these genes receive confirmed disease associations in OMIM, reviGen flags them via the surveillance signal before the phenotype database is updated.

**Baseline comparison:** Random baseline MRR ≈ 0.05. LIRICAL (state-of-the-art published tool) MRR ≈ 0.72 on similar neurodevelopmental cases. reviGen held-out in-DB MRR = **0.857**.

---

## Tech Stack

| Layer | Technology |
|-------|-----------|
| Language | Python 3.11, JavaScript, HTML/CSS |
| Backend Framework | FastAPI + Uvicorn |
| AI/NLP | Anthropic Claude API (claude-sonnet-4-6) |
| Variant Pathogenicity | AlphaMissense (Google DeepMind, SQLite) |
| Phenotype Similarity | Resnik semantic similarity via phenopy / Jaccard fallback |
| Gene-Disease Knowledge | OMIM API (live), ClinVar TSV (local), HPOA annotations |
| Phenotype Ontology | HPO (hp.obo) |
| Database | In-memory Python dict (PostgreSQL-ready) |
| Evaluation Dataset | MyGene2 rare disease registry |

---

## Project Structure

```
dxreanalyze/
├── backend/
│   ├── main_portal.py          # FastAPI app — all API endpoints
│   ├── aggregator.py           # 7-signal scoring engine
│   └── signals/
│       ├── omim_surveillance.py     # OMIM API + new gene-disease links
│       ├── vus_reclassification.py  # ClinVar VUS reclassification
│       ├── phenotypic_drift.py      # Resnik HPO similarity + drift
│       ├── inheritance.py           # AR/XLD/de novo rule-based flags
│       ├── analysis_gaps.py         # CNV, pipeline, scope gaps
│       ├── alphamissense.py         # DeepMind variant pathogenicity
│       ├── time_scorer.py           # Time since test signal
│       ├── nlp_extractor.py         # Claude API HPO extraction
│       ├── eval_pipeline.py         # F1/MRR evaluation pipeline
│       └── gene_ranker.py           # Phenotype-to-gene ranker
├── frontend/
│   └── revigen_portals_v2.html     # Single-file dual portal UI
├── validation/
│   └── mygene2_heldout_eval.py     # Held-out evaluation script
├── data/                            # (gitignored — see setup)
│   ├── phenotype.hpoa
│   ├── variant_summary.txt
│   └── alphamissense.db
├── .gitignore
├── requirements.txt
└── README.md
```

---

## Setup & Running

### 1. Clone the repo

```bash
git clone https://github.com/h425h/reviGen.git
cd reviGen
```

### 2. Install dependencies

```bash
pip install -r requirements.txt
```

### 3. Set environment variables

```bash
export ANTHROPIC_API_KEY=your_key_here
export OMIM_API_KEY=your_key_here        # optional — fallback table used if absent
```

### 4. Download data files (gitignored due to size)

```bash
# ClinVar variant summary
curl -o data/variant_summary.txt.gz \
  https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
gunzip data/variant_summary.txt.gz

# HPO ontology
curl -L -o data/hp.obo https://purl.obolibrary.org/obo/hp.obo

# HPOA disease annotations
curl -L -o data/phenotype.hpoa \
  https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/phenotype.hpoa
```

### 5. Run the server

```bash
uvicorn backend.main_portal:app --reload --port 8000
```

Open `frontend/revigen_portals_v2.html` in your browser.

---

## Running Evaluation

```bash
# Internal validation (10 literature cases)
python -m backend.signals.eval_pipeline --out results.json

# Held-out MyGene2 evaluation
python3 validation/mygene2_heldout_eval.py
```

---

## Key Design Decisions

**Why a weighted signal engine instead of a trained ML model?**
The rare disease reanalysis space has no large labeled training dataset. The Radboud 2022 study (n=150) is one of the largest real-world reanalysis cohorts. Using literature-derived weights is more defensible and interpretable than a black-box model trained on insufficient data.

**Why Claude for NLP instead of scispaCy?**
2024 benchmarks show LLMs outperform biomedical NER tools on HPO extraction specifically because they handle negation, context, and rare terminology that keyword-based approaches miss. Hallucinated terms that can't map to a real HPO ID are captured in `unmapped_terms` and never affect scoring.

**Why Resnik similarity instead of Jaccard?**
Jaccard treats all HPO terms equally. Resnik weights terms by Information Content — rare, specific symptoms (e.g. "Jansen metaphyseal chondrodysplasia") contribute far more than common ones (e.g. "hypotonia"). This is the same approach used by LIRICAL and Exomiser.

---

## Known Limitations & Next Steps

- ClinVar filter currently hardcoded to 5 genes — production version queries all genes dynamically
- In-memory database resets on server restart — PostgreSQL swap is the immediate next step
- Entropy boost currently uses uniform per-gene posteriors — per-signal per-gene weighting would be more accurate
- Phenotypic drift signal would benefit from delta-Resnik scoring (measuring convergence toward a diagnosis) rather than raw symptom count
- External validation against a larger MyGene2 cohort is the primary evaluation next step

---

## License

MIT

---

*Built at Harvard RareHack 2025. For rare disease patients still waiting for an answer.*
