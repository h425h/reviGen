import React, { useState } from 'react';

const API_BASE_URL = 'http://localhost:8000';

// ── Helpers ───────────────────────────────────────────────────────────────────
const scoreConfig = (score) => {
  if (score >= 80) return { ring: '#dc2626', bg: '#fef2f2', label: 'URGENT',         text: '#991b1b' };
  if (score >= 65) return { ring: '#ef4444', bg: '#fff1f1', label: 'HIGH PRIORITY',   text: '#b91c1c' };
  if (score >= 35) return { ring: '#f59e0b', bg: '#fffbeb', label: 'MEDIUM PRIORITY', text: '#92400e' };
  return               { ring: '#22c55e', bg: '#f0fdf4', label: 'LOW PRIORITY',     text: '#166534' };
};

const strengthStyle = (s) => ({
  HIGH:   { bg: '#fee2e2', color: '#991b1b', border: '#fca5a5' },
  MEDIUM: { bg: '#fef3c7', color: '#92400e', border: '#fcd34d' },
  LOW:    { bg: '#dcfce7', color: '#166534', border: '#86efac' },
}[s] || { bg: '#f3f4f6', color: '#374151', border: '#d1d5db' });

const inp = {
  width: '100%', boxSizing: 'border-box',
  border: '1px solid #d1d5db', borderRadius: 8,
  padding: '9px 12px', fontSize: 13, color: '#111827',
  background: '#fafafa', outline: 'none', fontFamily: 'inherit',
};
const lbl = {
  display: 'block', fontSize: 11, fontWeight: 700,
  color: '#6b7280', marginBottom: 5, letterSpacing: '0.07em', textTransform: 'uppercase',
};
const secLbl = {
  fontSize: 11, fontWeight: 700, color: '#6b7280',
  letterSpacing: '0.07em', textTransform: 'uppercase', marginBottom: 10,
};

// ── Sub-components ────────────────────────────────────────────────────────────
function Badge({ strength }) {
  const s = strengthStyle(strength);
  return (
    <span style={{
      background: s.bg, color: s.color, border: `1px solid ${s.border}`,
      borderRadius: 4, fontSize: 10, fontWeight: 700,
      padding: '2px 8px', letterSpacing: '0.08em', textTransform: 'uppercase',
    }}>{strength}</span>
  );
}

function SignalCard({ label, score, strength, detail, isNew }) {
  return (
    <div style={{
      background: '#fff', border: `1px solid ${isNew ? '#bfdbfe' : '#e5e7eb'}`,
      borderRadius: 10, padding: '12px 14px', position: 'relative',
      display: 'flex', flexDirection: 'column', gap: 5,
    }}>
      {isNew && (
        <span style={{
          position: 'absolute', top: -8, right: 10,
          background: '#2563eb', color: '#fff',
          fontSize: 9, fontWeight: 700, padding: '1px 7px', borderRadius: 20,
        }}>NEW</span>
      )}
      <div style={{ fontSize: 10, color: '#6b7280', fontWeight: 600, letterSpacing: '0.06em', textTransform: 'uppercase' }}>{label}</div>
      <div style={{ fontSize: 24, fontWeight: 800, color: '#111827', fontFamily: 'Georgia, serif' }}>
        {Math.round(score * 100)}<span style={{ fontSize: 12, color: '#9ca3af' }}>%</span>
      </div>
      <Badge strength={strength} />
      {detail && <div style={{ fontSize: 11, color: '#9ca3af', marginTop: 2 }}>{detail}</div>}
    </div>
  );
}

function ScoreGauge({ score }) {
  const cfg = scoreConfig(score);
  const r = 54, c = 2 * Math.PI * r, d = (score / 100) * c;
  return (
    <div style={{ display: 'flex', flexDirection: 'column', alignItems: 'center', gap: 8 }}>
      <svg width={140} height={140} viewBox="0 0 140 140">
        <circle cx={70} cy={70} r={r} fill="none" stroke="#f3f4f6" strokeWidth={10} />
        <circle cx={70} cy={70} r={r} fill="none" stroke={cfg.ring} strokeWidth={10}
          strokeDasharray={`${d} ${c - d}`} strokeLinecap="round"
          transform="rotate(-90 70 70)" style={{ transition: 'stroke-dasharray 0.8s ease' }} />
        <text x={70} y={64} textAnchor="middle" fontSize={30} fontWeight={800} fill="#111827" fontFamily="Georgia, serif">{Math.round(score)}</text>
        <text x={70} y={82} textAnchor="middle" fontSize={11} fill="#9ca3af">/100</text>
      </svg>
      <div style={{
        background: cfg.bg, color: cfg.text, border: `1px solid ${cfg.ring}30`,
        borderRadius: 20, padding: '4px 16px', fontSize: 11, fontWeight: 800, letterSpacing: '0.1em',
      }}>{cfg.label}</div>
    </div>
  );
}

function CollapsibleFlags({ flags, label, color }) {
  const [open, setOpen] = useState(false);
  if (!flags?.length) return null;
  return (
    <div style={{ marginTop: 6 }}>
      <button onClick={() => setOpen(!open)} style={{
        background: 'none', border: 'none', cursor: 'pointer',
        color: color || '#dc2626', fontSize: 11, fontWeight: 700, padding: 0,
      }}>
        {open ? '▾' : '▸'} {flags.length} {label}
      </button>
      {open && (
        <div style={{ marginTop: 6, display: 'flex', flexDirection: 'column', gap: 5 }}>
          {flags.map((f, i) => (
            <div key={i} style={{
              background: '#fafafa', border: '1px solid #e5e7eb',
              borderRadius: 7, padding: '8px 12px', fontSize: 11,
            }}>
              <span style={{ fontWeight: 700, color: '#111827' }}>{f.flag || f.gap || 'Signal'}</span>
              {f.gene && !['ALL','MULTIPLE'].includes(f.gene) && (
                <span style={{ marginLeft: 6, color: '#2563eb', fontWeight: 600 }}>({f.gene})</span>
              )}
              <div style={{ marginTop: 3, color: '#6b7280', lineHeight: 1.4 }}>
                {(f.reasoning || '').substring(0, 160)}{(f.reasoning || '').length > 160 ? '…' : ''}
              </div>
            </div>
          ))}
        </div>
      )}
    </div>
  );
}

// ── Evaluation Panel ──────────────────────────────────────────────────────────
function EvalPanel() {
  const [loading, setLoading] = useState(false);
  const [report,  setReport]  = useState(null);
  const [error,   setError]   = useState(null);

  // Single-case eval
  const [singleHpo,  setSingleHpo]  = useState('');
  const [singleGene, setSingleGene] = useState('');
  const [singleResult, setSingleResult] = useState(null);

  const runFullEval = async () => {
    setLoading(true); setError(null); setReport(null);
    try {
      const res = await fetch(`${API_BASE_URL}/api/evaluate?max_samples=10&top_k=10`);
      if (!res.ok) throw new Error(`API error: ${res.status}`);
      setReport(await res.json());
    } catch (e) { setError(e.message); }
    finally { setLoading(false); }
  };

  const runSingleEval = async () => {
    if (!singleHpo.trim() || !singleGene.trim()) return;
    setLoading(true); setError(null); setSingleResult(null);
    try {
      const hpoTerms = singleHpo.split(/[,\n]/).map(t => t.trim()).filter(Boolean);
      const res = await fetch(`${API_BASE_URL}/api/evaluate/single`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ hpo_terms: hpoTerms, true_gene: singleGene.trim().toUpperCase() }),
      });
      if (!res.ok) throw new Error(`API error: ${res.status}`);
      setSingleResult(await res.json());
    } catch (e) { setError(e.message); }
    finally { setLoading(false); }
  };

  const btnStyle = (active) => ({
    padding: '9px 18px', borderRadius: 7, border: 'none', cursor: 'pointer',
    fontWeight: 700, fontSize: 13,
    background: active ? '#2563eb' : '#e5e7eb',
    color: active ? '#fff' : '#6b7280',
  });

  return (
    <div style={{ display: 'flex', flexDirection: 'column', gap: 20 }}>

      {/* Single case */}
      <div style={{ background: '#f8fafc', border: '1px solid #e5e7eb', borderRadius: 10, padding: 18 }}>
        <div style={secLbl}>Single-Case Evaluation (Judge Demo)</div>
        <p style={{ fontSize: 12, color: '#6b7280', margin: '0 0 12px' }}>
          Enter HPO terms + known gene → see if the model predicts correctly. Returns F1 / Precision / Recall.
        </p>
        <div style={{ display: 'flex', flexDirection: 'column', gap: 10 }}>
          <div>
            <label style={lbl}>HPO Terms (comma-separated)</label>
            <input value={singleHpo} onChange={e => setSingleHpo(e.target.value)}
              placeholder="HP:0001263, HP:0002186, HP:0000252, HP:0001252, HP:0001344"
              style={inp} />
          </div>
          <div>
            <label style={lbl}>True Gene (ground truth)</label>
            <input value={singleGene} onChange={e => setSingleGene(e.target.value)}
              placeholder="FOXG1" style={{ ...inp, width: 200 }} />
          </div>
          <div style={{ display: 'flex', gap: 8, flexWrap: 'wrap' }}>
            <button onClick={runSingleEval} disabled={loading} style={btnStyle(true)}>
              {loading ? 'Running…' : 'Evaluate Single Case'}
            </button>
            {/* Quick-fill buttons */}
            {[
              { label: 'FOXG1', hpo: 'HP:0001263,HP:0002186,HP:0000252,HP:0001252,HP:0001332,HP:0001344,HP:0001270', gene: 'FOXG1' },
              { label: 'CDKL5', hpo: 'HP:0001250,HP:0001263,HP:0001252,HP:0000252,HP:0001344,HP:0001270', gene: 'CDKL5' },
              { label: 'MECP2', hpo: 'HP:0001263,HP:0002186,HP:0000252,HP:0001332,HP:0001344,HP:0001250', gene: 'MECP2' },
            ].map(d => (
              <button key={d.label} onClick={() => { setSingleHpo(d.hpo); setSingleGene(d.gene); }}
                style={{ ...btnStyle(false), fontSize: 11 }}>
                {d.label} demo
              </button>
            ))}
          </div>
        </div>

        {singleResult && (
          <div style={{ marginTop: 14 }}>
            <div style={{
              padding: '12px 16px', borderRadius: 8,
              background: singleResult.scores.correct_top1 ? '#f0fdf4' : '#fef2f2',
              border: `1px solid ${singleResult.scores.correct_top1 ? '#86efac' : '#fca5a5'}`,
              fontSize: 13, fontWeight: 700,
              color: singleResult.scores.correct_top1 ? '#166534' : '#991b1b',
            }}>
              {singleResult.verdict}
            </div>
            <div style={{ display: 'grid', gridTemplateColumns: 'repeat(3,1fr)', gap: 10, marginTop: 10 }}>
              {[
                { k: 1,  f1: singleResult.scores.f1_at_1,  p: singleResult.scores.precision_at_1,  r: singleResult.scores.recall_at_1 },
                { k: 3,  f1: singleResult.scores.f1_at_3,  p: singleResult.scores.precision_at_3,  r: singleResult.scores.recall_at_3 },
                { k: 10, f1: singleResult.scores.f1_at_10, p: singleResult.scores.precision_at_10, r: singleResult.scores.recall_at_10 },
              ].map(({ k, f1, p, r }) => (
                <div key={k} style={{ background: '#fff', border: '1px solid #e5e7eb', borderRadius: 8, padding: '10px 14px' }}>
                  <div style={{ fontSize: 10, color: '#6b7280', fontWeight: 700, marginBottom: 4 }}>@k={k}</div>
                  <div style={{ fontSize: 20, fontWeight: 800, color: f1 > 0 ? '#16a34a' : '#dc2626' }}>
                    F1: {f1.toFixed(2)}
                  </div>
                  <div style={{ fontSize: 11, color: '#6b7280' }}>P: {p.toFixed(2)} · R: {r.toFixed(2)}</div>
                </div>
              ))}
            </div>
            <div style={{ marginTop: 10 }}>
              <div style={secLbl}>Top-5 Predictions</div>
              {singleResult.top_predictions.map((p, i) => (
                <div key={i} style={{
                  display: 'flex', justifyContent: 'space-between', alignItems: 'center',
                  padding: '6px 10px', borderRadius: 6, marginBottom: 4,
                  background: p.gene === singleResult.true_gene ? '#f0fdf4' : '#fafafa',
                  border: `1px solid ${p.gene === singleResult.true_gene ? '#86efac' : '#e5e7eb'}`,
                  fontSize: 12,
                }}>
                  <span style={{ fontWeight: 700 }}>#{p.rank} {p.gene}</span>
                  <span style={{ color: '#6b7280' }}>{p.disease_name}</span>
                  <span style={{ fontWeight: 600, color: '#2563eb' }}>{(p.score * 100).toFixed(1)}%</span>
                </div>
              ))}
            </div>
          </div>
        )}
      </div>

      {/* Full eval */}
      <div style={{ background: '#f8fafc', border: '1px solid #e5e7eb', borderRadius: 10, padding: 18 }}>
        <div style={secLbl}>Full Evaluation (10 Built-in Test Cases)</div>
        <p style={{ fontSize: 12, color: '#6b7280', margin: '0 0 12px' }}>
          Runs against 10 literature-validated cases. Computes aggregate F1, Precision, Recall, MRR.
        </p>
        <button onClick={runFullEval} disabled={loading} style={btnStyle(true)}>
          {loading ? 'Running…' : 'Run Full Evaluation'}
        </button>

        {error && (
          <div style={{ marginTop: 12, background: '#fef2f2', border: '1px solid #fca5a5', borderRadius: 8, padding: 12, fontSize: 12, color: '#991b1b' }}>
            {error}
          </div>
        )}

        {report && (
          <div style={{ marginTop: 14 }}>
            <div style={{ display: 'grid', gridTemplateColumns: 'repeat(3,1fr)', gap: 10, marginBottom: 14 }}>
              {[
                { k: 1,  f1: report.summary.f1_at1,  p: report.summary.precision_at1,  r: report.summary.recall_at1  },
                { k: 3,  f1: report.summary.f1_at3,  p: report.summary.precision_at3,  r: report.summary.recall_at3  },
                { k: 10, f1: report.summary.f1_at10, p: report.summary.precision_at10, r: report.summary.recall_at10 },
              ].map(({ k, f1, p, r }) => (
                <div key={k} style={{ background: '#fff', border: '1px solid #e5e7eb', borderRadius: 8, padding: '12px 14px' }}>
                  <div style={{ fontSize: 10, color: '#6b7280', fontWeight: 700, marginBottom: 4 }}>F1 @k={k}</div>
                  <div style={{ fontSize: 28, fontWeight: 800, color: '#111827', fontFamily: 'Georgia,serif' }}>{f1.toFixed(2)}</div>
                  <div style={{ fontSize: 11, color: '#6b7280' }}>P={p.toFixed(2)} · R={r.toFixed(2)}</div>
                </div>
              ))}
            </div>
            <div style={{ fontSize: 12, color: '#6b7280', marginBottom: 10 }}>
              MRR: <strong>{report.summary.mean_reciprocal_rank?.toFixed(3)}</strong> ·
              Method: <strong>{report.summary.similarity_method}</strong> ·
              n={report.summary.n_samples}
            </div>
            <div style={secLbl}>Per-sample results</div>
            {report.per_sample.map((r, i) => (
              <div key={i} style={{
                display: 'flex', justifyContent: 'space-between', alignItems: 'center',
                padding: '6px 10px', borderRadius: 6, marginBottom: 4,
                background: r.correct_top1 ? '#f0fdf4' : r.correct_top3 ? '#fffbeb' : '#fef2f2',
                border: `1px solid ${r.correct_top1 ? '#86efac' : r.correct_top3 ? '#fcd34d' : '#fca5a5'}`,
                fontSize: 11,
              }}>
                <span style={{ fontWeight: 700, minWidth: 120 }}>{r.sample_id}</span>
                <span>true: <strong>{r.true_gene}</strong></span>
                <span>pred: <strong>{r.top1_prediction || '—'}</strong></span>
                <span style={{ color: r.correct_top1 ? '#16a34a' : '#dc2626', fontWeight: 700 }}>
                  {r.correct_top1 ? '✓ top-1' : r.correct_top3 ? '~ top-3' : r.correct_top10 ? '~ top-10' : '✗'}
                </span>
              </div>
            ))}
          </div>
        )}
      </div>
    </div>
  );
}

// ── Main App ──────────────────────────────────────────────────────────────────
export default function App() {
  const [activeTab, setActiveTab] = useState('structured'); // 'structured' | 'nlp' | 'eval'

  // Structured input state
  const [originalHpo,   setOriginalHpo]   = useState('');
  const [currentHpo,    setCurrentHpo]    = useState('');
  const [testType,      setTestType]      = useState('exome');
  const [testDate,      setTestDate]      = useState('');
  const [testResult,    setTestResult]    = useState('negative');
  const [vusGenes,      setVusGenes]      = useState('');
  const [analysisType,  setAnalysisType]  = useState('');
  const [patientSex,    setPatientSex]    = useState('');
  const [consanguineous,setConsanguineous]= useState('');
  const [cnvCalling,    setCnvCalling]    = useState('');
  const [showAdvanced,  setShowAdvanced]  = useState(false);

  // NLP input state
  const [originalText,  setOriginalText]  = useState('');
  const [currentText,   setCurrentText]   = useState('');

  // Shared
  const [loading, setLoading] = useState(false);
  const [results, setResults] = useState(null);
  const [error,   setError]   = useState(null);
  const [checked, setChecked] = useState({});

  const parseTerms = (t) => t.split(/[,\n]/).map(x => x.trim()).filter(Boolean);

  const handleStructuredAnalyze = async () => {
    setLoading(true); setError(null); setResults(null); setChecked({});
    try {
      const res = await fetch(`${API_BASE_URL}/api/reanalysis-score`, {
        method: 'POST', headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          original_hpo_terms: parseTerms(originalHpo),
          current_hpo_terms:  parseTerms(currentHpo),
          prior_test: {
            test_type: testType, test_date: testDate, result: testResult,
            vus_genes: parseTerms(vusGenes).map(g => g.toUpperCase()),
            analysis_type: analysisType || null,
          },
          clinical_context: {
            patient_sex: patientSex || null,
            consanguineous: consanguineous === 'true' ? true : consanguineous === 'false' ? false : null,
            cnv_calling_performed: cnvCalling === 'true' ? true : cnvCalling === 'false' ? false : null,
          },
        }),
      });
      if (!res.ok) throw new Error(`API error: ${res.status}`);
      setResults(await res.json());
    } catch (e) { setError(e.message); }
    finally { setLoading(false); }
  };

  const handleNlpAnalyze = async () => {
    if (!originalText.trim() || !currentText.trim() || !testDate) {
      setError('Original notes, current notes, and test date are required.'); return;
    }
    setLoading(true); setError(null); setResults(null); setChecked({});
    try {
      const res = await fetch(`${API_BASE_URL}/api/extract-and-score`, {
        method: 'POST', headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          original_clinical_text: originalText,
          current_clinical_text:  currentText,
          prior_test: {
            test_type: testType, test_date: testDate, result: testResult,
            vus_genes: parseTerms(vusGenes).map(g => g.toUpperCase()),
            analysis_type: analysisType || null,
          },
          clinical_context: {
            patient_sex: patientSex || null,
            cnv_calling_performed: cnvCalling === 'true' ? true : cnvCalling === 'false' ? false : null,
          },
        }),
      });
      if (!res.ok) throw new Error(`API error: ${res.status}`);
      setResults(await res.json());
    } catch (e) { setError(e.message); }
    finally { setLoading(false); }
  };

  // Demo fills
  const demos = [
    {
      label: 'FOXG1 — gene panel 2019',
      fill: () => {
        setOriginalHpo('HP:0001252, HP:0001270, HP:0011968');
        setCurrentHpo('HP:0001252, HP:0001270, HP:0011968, HP:0001263, HP:0002186, HP:0000729, HP:0002072, HP:0001260');
        setOriginalText('Infant with hypotonia, motor delay, and feeding difficulties at 6 months.');
        setCurrentText('7-year-old with global developmental delay, absent speech, apraxia, autistic behaviour, choreiform movements, and dysarthria. Progressive course since 18 months.');
        setTestType('gene panel'); setTestDate('2019-06-01'); setTestResult('negative');
        setVusGenes('FOXG1'); setAnalysisType('singleton'); setPatientSex('F'); setCnvCalling('false');
      },
    },
    {
      label: 'DHTKD1 (AR) — 2021 singleton',
      fill: () => {
        setOriginalHpo('HP:0001252, HP:0000252, HP:0001263');
        setCurrentHpo('HP:0001252, HP:0000252, HP:0001263, HP:0003355, HP:0011968');
        setOriginalText('Young child with hypotonia, microcephaly, and global developmental delay.');
        setCurrentText('School-age child with hypotonia, microcephaly, global developmental delay, aminoaciduria, and feeding difficulties.');
        setTestType('exome'); setTestDate('2021-03-15'); setTestResult('vus');
        setVusGenes('DHTKD1'); setAnalysisType('singleton'); setPatientSex('M'); setCnvCalling('false');
      },
    },
    {
      label: "Nina Nazar — PTHR1 2015",
      fill: () => {
        setOriginalHpo('HP:0003521, HP:0002750, HP:0000843');
        setCurrentHpo('HP:0003521, HP:0002750, HP:0000843, HP:0000518, HP:0002748, HP:0003075');
        setOriginalText('Adult with disproportionate short stature, delayed bone maturation, and hyperparathyroidism.');
        setCurrentText('Adult with short stature, skeletal dysplasia, hyperparathyroidism, cataracts, rickets, and hypophosphataemia.');
        setTestType('exome'); setTestDate('2015-01-01'); setTestResult('negative');
        setVusGenes('PTHR1'); setAnalysisType('singleton'); setPatientSex('F'); setCnvCalling('false');
      },
    },
  ];

  const sb = results?.signal_breakdown;

  const tabStyle = (t) => ({
    padding: '8px 20px', borderRadius: '8px 8px 0 0',
    border: '1px solid #e5e7eb', borderBottom: activeTab === t ? '1px solid #fff' : '1px solid #e5e7eb',
    background: activeTab === t ? '#fff' : '#f8fafc',
    cursor: 'pointer', fontWeight: 700, fontSize: 13,
    color: activeTab === t ? '#2563eb' : '#6b7280',
    marginRight: 4,
  });

  return (
    <div style={{ minHeight: '100vh', background: '#f8fafc', fontFamily: "'DM Sans','Segoe UI',sans-serif" }}>

      {/* Header */}
      <div style={{
        background: '#fff', borderBottom: '1px solid #e5e7eb', padding: '0 32px',
        display: 'flex', alignItems: 'center', justifyContent: 'space-between', height: 64,
      }}>
        <div style={{ display: 'flex', alignItems: 'baseline', gap: 10 }}>
          <span style={{ fontSize: 22, fontWeight: 900, color: '#0f172a', fontFamily: 'Georgia,serif', letterSpacing: '-0.5px' }}>
            Revi<span style={{ color: '#2563eb' }}>gen</span>
          </span>
          <span style={{ fontSize: 12, color: '#94a3b8', fontWeight: 500 }}>Reviewing & Reviving Rare Diseases</span>
        </div>
        <div style={{ display: 'flex', alignItems: 'center', gap: 8 }}>
          <div style={{ width: 8, height: 8, borderRadius: '50%', background: '#22c55e' }} />
          <span style={{ fontSize: 12, color: '#6b7280' }}>HackRare 2026</span>
        </div>
      </div>

      <div style={{ maxWidth: 1200, margin: '0 auto', padding: '32px 24px' }}>

        {/* Hero */}
        <div style={{ marginBottom: 24 }}>
          <h1 style={{ fontSize: 26, fontWeight: 800, color: '#0f172a', margin: 0, fontFamily: 'Georgia,serif' }}>
            Genetic Reanalysis Trigger
          </h1>
          <p style={{ fontSize: 14, color: '#64748b', marginTop: 6, maxWidth: 620 }}>
            7-signal engine · OMIM surveillance · Resnik IC-weighted phenotype matching ·
            Shannon entropy for multi-VUS uncertainty · NLP extraction · F1-evaluated
          </p>
        </div>

        {/* Tabs */}
        <div style={{ marginBottom: -1, position: 'relative', zIndex: 1 }}>
          <button style={tabStyle('structured')} onClick={() => setActiveTab('structured')}>
            📋 Structured Input
          </button>
          <button style={tabStyle('nlp')} onClick={() => setActiveTab('nlp')}>
            🧠 Clinical Notes (NLP)
          </button>
          <button style={tabStyle('eval')} onClick={() => setActiveTab('eval')}>
            📊 Evaluation / F1
          </button>
        </div>

        {/* Evaluation tab */}
        {activeTab === 'eval' && (
          <div style={{ background: '#fff', border: '1px solid #e5e7eb', borderRadius: '0 10px 10px 10px', padding: 28 }}>
            <EvalPanel />
          </div>
        )}

        {/* Input + Results (shared between structured and NLP tabs) */}
        {activeTab !== 'eval' && (
          <div style={{ display: 'grid', gridTemplateColumns: '1fr 1.1fr', gap: 24, alignItems: 'start' }}>

            {/* ── LEFT: Input Panel ── */}
            <div style={{ background: '#fff', border: '1px solid #e5e7eb', borderRadius: '0 10px 10px 10px', padding: 28 }}>

              {/* Structured input */}
              {activeTab === 'structured' && (
                <div style={{ display: 'flex', flexDirection: 'column', gap: 14 }}>
                  <div>
                    <label style={lbl}>HPO Terms — at original test</label>
                    <textarea value={originalHpo} onChange={e => setOriginalHpo(e.target.value)}
                      placeholder="HP:0001252, HP:0001270, HP:0011968"
                      rows={3} style={{ ...inp, resize: 'vertical' }} />
                  </div>
                  <div>
                    <label style={lbl}>HPO Terms — current</label>
                    <textarea value={currentHpo} onChange={e => setCurrentHpo(e.target.value)}
                      placeholder="HP:0001252, HP:0001270, HP:0001263, HP:0002186"
                      rows={3} style={{ ...inp, resize: 'vertical' }} />
                  </div>
                </div>
              )}

              {/* NLP input */}
              {activeTab === 'nlp' && (
                <div style={{ display: 'flex', flexDirection: 'column', gap: 14 }}>
                  <div style={{ background: '#eff6ff', border: '1px solid #bfdbfe', borderRadius: 8, padding: '10px 14px', fontSize: 12, color: '#1e3a5f' }}>
                    ✨ <strong>NLP mode:</strong> Paste free-text clinical notes. Claude extracts HPO terms automatically.
                  </div>
                  <div>
                    <label style={lbl}>Clinical Notes — at time of original test</label>
                    <textarea value={originalText} onChange={e => setOriginalText(e.target.value)}
                      placeholder="e.g. Infant with hypotonia, motor delay, and feeding difficulties at 6 months..."
                      rows={4} style={{ ...inp, resize: 'vertical' }} />
                  </div>
                  <div>
                    <label style={lbl}>Clinical Notes — current presentation</label>
                    <textarea value={currentText} onChange={e => setCurrentText(e.target.value)}
                      placeholder="e.g. 7-year-old with absent speech, hand stereotypies, seizures, and progressive regression..."
                      rows={4} style={{ ...inp, resize: 'vertical' }} />
                  </div>
                </div>
              )}

              {/* Shared test metadata */}
              <div style={{ marginTop: 14, display: 'flex', flexDirection: 'column', gap: 12 }}>
                <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: 12 }}>
                  <div>
                    <label style={lbl}>Test Type</label>
                    <select value={testType} onChange={e => setTestType(e.target.value)} style={inp}>
                      <option value="exome">Exome Sequencing</option>
                      <option value="gene panel">Gene Panel</option>
                      <option value="genome">Genome Sequencing</option>
                      <option value="other">Other</option>
                    </select>
                  </div>
                  <div>
                    <label style={lbl}>Test Date</label>
                    <input type="date" value={testDate} onChange={e => setTestDate(e.target.value)} style={inp} />
                  </div>
                </div>
                <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: 12 }}>
                  <div>
                    <label style={lbl}>Result</label>
                    <select value={testResult} onChange={e => setTestResult(e.target.value)} style={inp}>
                      <option value="negative">Negative</option>
                      <option value="vus">VUS Identified</option>
                      <option value="inconclusive">Inconclusive</option>
                    </select>
                  </div>
                  <div>
                    <label style={lbl}>Analysis Type</label>
                    <select value={analysisType} onChange={e => setAnalysisType(e.target.value)} style={inp}>
                      <option value="">Unknown</option>
                      <option value="singleton">Singleton</option>
                      <option value="duo">Duo</option>
                      <option value="trio">Trio</option>
                    </select>
                  </div>
                </div>
                <div>
                  <label style={lbl}>VUS Genes</label>
                  <input value={vusGenes} onChange={e => setVusGenes(e.target.value)}
                    placeholder="FOXG1, CDKL5" style={inp} />
                </div>

                {/* Advanced context */}
                <div style={{ borderTop: '1px solid #f1f5f9', paddingTop: 10 }}>
                  <button onClick={() => setShowAdvanced(!showAdvanced)} style={{
                    background: 'none', border: 'none', cursor: 'pointer',
                    fontSize: 12, color: '#2563eb', fontWeight: 700, padding: 0,
                  }}>
                    {showAdvanced ? '▾' : '▸'} Clinical Context (improves inheritance & gap signals)
                  </button>
                  {showAdvanced && (
                    <div style={{ marginTop: 10, display: 'grid', gridTemplateColumns: '1fr 1fr', gap: 10 }}>
                      <div>
                        <label style={lbl}>Sex</label>
                        <select value={patientSex} onChange={e => setPatientSex(e.target.value)} style={inp}>
                          <option value="">Unknown</option>
                          <option value="M">Male</option>
                          <option value="F">Female</option>
                        </select>
                      </div>
                      <div>
                        <label style={lbl}>CNV Calling?</label>
                        <select value={cnvCalling} onChange={e => setCnvCalling(e.target.value)} style={inp}>
                          <option value="">Unknown</option>
                          <option value="true">Yes</option>
                          <option value="false">No</option>
                        </select>
                      </div>
                    </div>
                  )}
                </div>

                <button
                  onClick={activeTab === 'nlp' ? handleNlpAnalyze : handleStructuredAnalyze}
                  disabled={loading || !testDate}
                  style={{
                    width: '100%', padding: '12px 0', borderRadius: 8, border: 'none',
                    background: loading || !testDate ? '#e5e7eb' : '#2563eb',
                    color: loading || !testDate ? '#9ca3af' : '#fff',
                    fontSize: 14, fontWeight: 700, cursor: loading || !testDate ? 'not-allowed' : 'pointer',
                  }}>
                  {loading ? 'Analyzing…' : activeTab === 'nlp' ? '🧠 Extract HPO + Run Score' : 'Run Reanalysis Score'}
                </button>
              </div>

              {/* Demo cases */}
              <div style={{ marginTop: 16, paddingTop: 14, borderTop: '1px solid #f1f5f9' }}>
                <div style={{ fontSize: 11, color: '#94a3b8', fontWeight: 600, marginBottom: 8, letterSpacing: '0.06em', textTransform: 'uppercase' }}>
                  Demo Cases
                </div>
                {demos.map((d, i) => (
                  <button key={i} onClick={d.fill} style={{
                    width: '100%', padding: '7px 12px', marginBottom: 5,
                    background: '#eff6ff', border: '1px solid #bfdbfe', borderRadius: 6,
                    color: '#1d4ed8', fontSize: 12, fontWeight: 600, cursor: 'pointer', textAlign: 'left',
                  }}>
                    ▶ {d.label}
                  </button>
                ))}
              </div>
            </div>

            {/* ── RIGHT: Results Panel ── */}
            <div style={{ background: '#fff', border: '1px solid #e5e7eb', borderRadius: 10, padding: 28, minHeight: 400 }}>
              <h2 style={{ fontSize: 14, fontWeight: 800, color: '#0f172a', margin: '0 0 20px' }}>
                Analysis Results
              </h2>

              {!loading && !error && !results && (
                <div style={{ display: 'flex', flexDirection: 'column', alignItems: 'center', justifyContent: 'center', height: 320, gap: 12 }}>
                  <div style={{ fontSize: 36 }}>🧬</div>
                  <div style={{ fontSize: 14, color: '#94a3b8', textAlign: 'center', maxWidth: 260 }}>
                    Run analysis to see the 7-signal reanalysis score
                  </div>
                </div>
              )}

              {loading && (
                <div style={{ display: 'flex', flexDirection: 'column', alignItems: 'center', justifyContent: 'center', height: 320, gap: 16 }}>
                  <div style={{ width: 40, height: 40, border: '3px solid #e5e7eb', borderTop: '3px solid #2563eb', borderRadius: '50%', animation: 'spin 0.8s linear infinite' }} />
                  <div style={{ fontSize: 13, color: '#6b7280' }}>
                    {activeTab === 'nlp' ? 'Extracting HPO terms + scoring…' : 'Querying OMIM, ClinVar & HPO databases…'}
                  </div>
                  <style>{`@keyframes spin { to { transform: rotate(360deg); } }`}</style>
                </div>
              )}

              {error && (
                <div style={{ background: '#fef2f2', border: '1px solid #fca5a5', borderRadius: 8, padding: 14, color: '#991b1b', fontSize: 13 }}>
                  {error}
                </div>
              )}

              {results && (
                <div style={{ display: 'flex', flexDirection: 'column', gap: 20 }}>

                  {/* NLP extraction summary */}
                  {results.nlp_extraction && (
                    <div style={{ background: '#eff6ff', border: '1px solid #bfdbfe', borderRadius: 8, padding: '10px 14px', fontSize: 12 }}>
                      <div style={{ fontWeight: 700, color: '#2563eb', marginBottom: 4 }}>✨ NLP Extraction</div>
                      <div style={{ color: '#1e3a5f' }}>
                        Original: <strong>{results.nlp_extraction.original_term_count} HPO terms</strong> ·
                        Current: <strong>{results.nlp_extraction.current_term_count} HPO terms</strong> ·
                        Method: <strong>{results.nlp_extraction.original_method}</strong>
                      </div>
                      <div style={{ marginTop: 4, color: '#6b7280', fontSize: 11 }}>
                        {results.nlp_extraction.current_hpo_terms.join(', ')}
                      </div>
                    </div>
                  )}

                  {/* Score Gauge */}
                  <div style={{ display: 'flex', flexDirection: 'column', alignItems: 'center', gap: 4 }}>
                    <ScoreGauge score={results.reanalysis_score} />
                    {results.urgency && (
                      <div style={{ fontSize: 12, color: '#6b7280', fontStyle: 'italic' }}>{results.urgency}</div>
                    )}
                    {results.entropy_boost > 0 && (
                      <div style={{ fontSize: 11, color: '#d97706', fontWeight: 600 }}>
                        +{(results.entropy_boost * 100).toFixed(1)} pts entropy boost
                        (H={sb?.entropy_modifier?.entropy_value?.toFixed(2)})
                      </div>
                    )}
                  </div>

                  {/* Narrative */}
                  {results.clinical_recommendation?.narrative && (
                    <div style={{ background: '#eff6ff', border: '1px solid #bfdbfe', borderRadius: 10, padding: 14 }}>
                      <div style={{ fontSize: 10, fontWeight: 700, color: '#2563eb', letterSpacing: '0.08em', textTransform: 'uppercase', marginBottom: 6 }}>
                        Clinical Summary
                      </div>
                      <p style={{ margin: 0, fontSize: 13, color: '#1e3a5f', lineHeight: 1.6 }}>
                        {results.clinical_recommendation.narrative}
                      </p>
                      {results.clinical_recommendation.recommended_reanalysis_type && (
                        <div style={{ marginTop: 8, fontSize: 11, fontWeight: 700, color: '#1d4ed8' }}>
                          Recommended: {results.clinical_recommendation.recommended_reanalysis_type.replace(/_/g, ' ').toUpperCase()}
                        </div>
                      )}
                    </div>
                  )}

                  {/* 7-signal grid (2×3 + entropy) */}
                  <div>
                    <div style={secLbl}>Signal Breakdown</div>
                    <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr 1fr', gap: 8 }}>
                      <SignalCard label="OMIM Gene-Disease" isNew
                        score={sb?.omim_gene_disease?.score || 0}
                        strength={sb?.omim_gene_disease?.signal_strength || 'LOW'}
                        detail={`${sb?.omim_gene_disease?.signal_count || 0} gene(s)`} />
                      <SignalCard label="VUS Reclassification"
                        score={sb?.vus_reclassification?.score || 0}
                        strength={sb?.vus_reclassification?.signal_strength || 'LOW'}
                        detail={`${sb?.vus_reclassification?.signal_count || 0} ClinVar signal(s)`} />
                      <SignalCard label="Disease Match (Resnik)"
                        score={sb?.disease_match?.score || 0}
                        strength={sb?.disease_match?.signal_strength || 'LOW'}
                        detail={sb?.disease_match?.similarity_method === 'resnik' ? 'IC-weighted' : 'Jaccard fallback'} />
                      <SignalCard label="Inheritance Flags" isNew
                        score={sb?.inheritance_pattern?.score || 0}
                        strength={sb?.inheritance_pattern?.signal_strength || 'LOW'}
                        detail={`${sb?.inheritance_pattern?.flag_count || 0} flag(s)`} />
                      <SignalCard label="Analysis Gaps" isNew
                        score={sb?.analysis_method_gaps?.score || 0}
                        strength={sb?.analysis_method_gaps?.signal_strength || 'LOW'}
                        detail={`${sb?.analysis_method_gaps?.gap_count || 0} gap(s)`} />
                      <SignalCard label="Phenotypic Drift"
                        score={sb?.phenotypic_drift?.score || 0}
                        strength={sb?.phenotypic_drift?.details?.signal_strength || 'LOW'}
                        detail={`${sb?.phenotypic_drift?.details?.new_symptom_count || 0} new symptoms`} />
                    </div>
                    {/* Entropy modifier */}
                    {sb?.entropy_modifier?.entropy_value > 0 && (
                      <div style={{
                        marginTop: 8, background: '#fffbeb', border: '1px solid #fcd34d',
                        borderRadius: 8, padding: '8px 12px', fontSize: 11, color: '#92400e',
                      }}>
                        <strong>Entropy modifier:</strong> H={sb.entropy_modifier.entropy_value.toFixed(3)} →
                        +{(sb.entropy_modifier.entropy_boost * 100).toFixed(1)} pts
                        ({sb.entropy_modifier.n_vus_genes} VUS genes, uncertain which is causal)
                      </div>
                    )}
                  </div>

                  {/* OMIM signals */}
                  {sb?.omim_gene_disease?.signals?.length > 0 && (
                    <div>
                      <div style={secLbl}>OMIM Gene-Disease Signals</div>
                      {sb.omim_gene_disease.signals.map((sig, i) => (
                        <div key={i} style={{
                          background: sig.is_new_after_test ? '#fefce8' : '#f0fdf4',
                          border: `1px solid ${sig.is_new_after_test ? '#fde047' : '#86efac'}`,
                          borderRadius: 8, padding: '8px 12px', marginBottom: 6, fontSize: 12,
                        }}>
                          <strong>{sig.gene}</strong> → {sig.disease}
                          {sig.is_new_after_test && <span style={{ marginLeft: 8, color: '#d97706', fontWeight: 700, fontSize: 10 }}>⚡ NEW AFTER TEST</span>}
                          <div style={{ fontSize: 11, color: '#6b7280', marginTop: 2 }}>{sig.reasoning}</div>
                        </div>
                      ))}
                    </div>
                  )}

                  {/* Inheritance + gaps (collapsible) */}
                  {sb?.inheritance_pattern?.flags?.length > 0 && (
                    <div>
                      <div style={secLbl}>Inheritance Flags</div>
                      <CollapsibleFlags flags={sb.inheritance_pattern.flags} label="inheritance flags" color="#dc2626" />
                    </div>
                  )}
                  {sb?.analysis_method_gaps?.gaps?.length > 0 && (
                    <div>
                      <div style={secLbl}>Analysis Gaps</div>
                      <CollapsibleFlags flags={sb.analysis_method_gaps.gaps.map(g => ({ ...g, flag: g.gap }))} label="gaps" color="#d97706" />
                    </div>
                  )}

                  {/* New symptoms */}
                  {sb?.phenotypic_drift?.details?.new_symptoms?.length > 0 && (
                    <div>
                      <div style={secLbl}>New Symptoms</div>
                      <div style={{ display: 'flex', flexWrap: 'wrap', gap: 6 }}>
                        {sb.phenotypic_drift.details.new_symptoms.map((s, i) => (
                          <span key={i} style={{
                            background: '#eff6ff', color: '#1d4ed8', border: '1px solid #bfdbfe',
                            borderRadius: 20, fontSize: 11, fontWeight: 600, padding: '4px 10px',
                          }}>{s.name || s.hpo_id}</span>
                        ))}
                      </div>
                    </div>
                  )}

                  {/* Why + Checklist */}
                  {results.clinical_recommendation?.top_reasons && (
                    <div>
                      <div style={secLbl}>Why Reanalysis?</div>
                      {results.clinical_recommendation.top_reasons.map((r, i) => (
                        <div key={i} style={{ display: 'flex', gap: 8, fontSize: 13, color: '#374151', marginBottom: 4 }}>
                          <span style={{ color: '#2563eb', fontWeight: 700 }}>→</span>{r}
                        </div>
                      ))}
                    </div>
                  )}
                  {results.clinical_recommendation?.checklist && (
                    <div>
                      <div style={secLbl}>Reanalysis Checklist</div>
                      {results.clinical_recommendation.checklist.map((item, i) => (
                        <label key={i} style={{ display: 'flex', gap: 10, alignItems: 'flex-start', cursor: 'pointer', marginBottom: 6 }}>
                          <input type="checkbox" checked={checked[i] || false}
                            onChange={() => setChecked(p => ({ ...p, [i]: !p[i] }))}
                            style={{ marginTop: 2, accentColor: '#2563eb' }} />
                          <span style={{ fontSize: 13, color: checked[i] ? '#9ca3af' : '#374151', textDecoration: checked[i] ? 'line-through' : 'none' }}>
                            {item}
                          </span>
                        </label>
                      ))}
                    </div>
                  )}

                  {/* Weight disclosure */}
                  <details style={{ fontSize: 11, color: '#9ca3af' }}>
                    <summary style={{ cursor: 'pointer', fontWeight: 600 }}>Signal weights (literature-backed)</summary>
                    <div style={{ marginTop: 8, display: 'flex', flexDirection: 'column', gap: 2 }}>
                      {[['OMIM gene-disease','29%','Radboud 2022: 42% of yield'],['ClinVar VUS','14%','10% of yield'],
                        ['Disease match (Resnik)','13%','IC-weighted similarity'],['Phenotypic drift','11%','Raw symptom change'],
                        ['Inheritance flags','12%','AR 2nd hit, de novo, XLD female'],['Analysis gaps','10%','19% bioinformatics, NEJM CNV'],
                        ['AlphaMissense','5%','Structure-informed (if coords)'],['Time since test','3%','Mechanism above'],
                        ['Entropy modifier','≤8pts','Multi-VUS uncertainty H=-Σp·log₂p']
                      ].map(([l,w,n]) => (
                        <div key={l} style={{ display: 'flex', gap: 8 }}>
                          <span style={{ fontWeight: 700, minWidth: 175 }}>{l}</span>
                          <span style={{ color: '#2563eb', minWidth: 40 }}>{w}</span>
                          <span>{n}</span>
                        </div>
                      ))}
                    </div>
                  </details>

                </div>
              )}
            </div>
          </div>
        )}

        <div style={{ textAlign: 'center', marginTop: 32, fontSize: 12, color: '#94a3b8' }}>
          Revigen · HackRare 2026 · Built for patients like Nina Nazar · 32 years is too long to wait
        </div>
      </div>
    </div>
  );
}