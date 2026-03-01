"""
AlphaMissense Variant Pathogenicity Signal — NEW SIGNAL.

AlphaMissense (Google DeepMind, 2023) predicts missense variant pathogenicity
using protein structure context (AlphaFold2 embeddings).

Performance: 91.5% accuracy classifying known pathogenic vs benign variants.
Outperforms CADD, REVEL, and prior ML tools on rare disease benchmarks.

Thresholds (from AlphaMissense paper, Table S2):
  am_pathogenicity ≥ 0.564 → likely_pathogenic   (91.5% PPV)
  am_pathogenicity ≤ 0.340 → likely_benign        (90.1% NPV)
  Between             → ambiguous

Data: CC-BY licensed, freely downloadable.
  hg38: https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz (~1.5 GB)
  hg19: https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg19.tsv.gz

Setup (one-time, run before using this module):
  python -c "from signals.alphamissense import setup_alphamissense_db; setup_alphamissense_db()"

The TSV is loaded into a local SQLite database for fast O(1) coordinate lookups.
No external API calls at query time — everything runs offline.

Column names in the TSV:
  #CHROM  POS  REF  ALT  genome  uniprot_id  transcript_id  protein_variant  am_pathogenicity  am_class
"""

import os
import sqlite3
import gzip
import csv
from typing import Optional, Dict, List


# Path to the SQLite database (created from AlphaMissense TSV)
DEFAULT_DB_PATH = os.path.join(
    os.path.dirname(os.path.dirname(__file__)),
    "data", "alphamissense.db"
)

# Pathogenicity thresholds from the AlphaMissense paper
PATHOGENIC_THRESHOLD = 0.564
BENIGN_THRESHOLD     = 0.340


def setup_alphamissense_db(
    tsv_gz_path: str,
    db_path: str = DEFAULT_DB_PATH
) -> bool:
    """
    One-time setup: load AlphaMissense TSV.gz into SQLite.
    Creates indexed table for fast coordinate lookups.
    
    Args:
        tsv_gz_path: Path to AlphaMissense_hg38.tsv.gz
        db_path:     Output SQLite database path
    
    Returns:
        True on success
    """
    print(f"Loading AlphaMissense from {tsv_gz_path}...")
    print("This takes ~10 minutes for the full dataset. Do this once.")

    os.makedirs(os.path.dirname(db_path), exist_ok=True)
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    cursor.execute("""
        CREATE TABLE IF NOT EXISTS alphamissense (
            chrom TEXT,
            pos   INTEGER,
            ref   TEXT,
            alt   TEXT,
            am_pathogenicity REAL,
            am_class TEXT,
            protein_variant TEXT,
            uniprot_id TEXT
        )
    """)
    cursor.execute("DELETE FROM alphamissense")

    with gzip.open(tsv_gz_path, 'rt', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        batch = []
        for i, row in enumerate(reader):
            # Skip comment/header rows (TSV has a #CHROM header)
            if row.get('#CHROM', '').startswith('#'):
                continue
            try:
                batch.append((
                    row.get('#CHROM', row.get('CHROM', '')),
                    int(row['POS']),
                    row['REF'],
                    row['ALT'],
                    float(row['am_pathogenicity']),
                    row['am_class'],
                    row.get('protein_variant', ''),
                    row.get('uniprot_id', ''),
                ))
            except (ValueError, KeyError):
                continue

            if len(batch) >= 100_000:
                cursor.executemany(
                    "INSERT INTO alphamissense VALUES (?,?,?,?,?,?,?,?)",
                    batch
                )
                batch = []
                if i % 1_000_000 == 0:
                    print(f"  {i:,} rows loaded...")
                    conn.commit()

        if batch:
            cursor.executemany(
                "INSERT INTO alphamissense VALUES (?,?,?,?,?,?,?,?)",
                batch
            )

    print("Creating index...")
    cursor.execute(
        "CREATE INDEX IF NOT EXISTS idx_coords ON alphamissense(chrom, pos, ref, alt)"
    )
    conn.commit()
    conn.close()
    print(f"AlphaMissense database ready at {db_path}")
    return True


def lookup_alphamissense(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    db_path: str = DEFAULT_DB_PATH
) -> Optional[Dict]:
    """
    Look up a single variant in the AlphaMissense SQLite database.
    
    Args:
        chrom: Chromosome (e.g., 'chr14' or '14')
        pos:   Position (GRCh38)
        ref:   Reference allele
        alt:   Alternate allele
        db_path: Path to SQLite database
    
    Returns:
        Dict with am_pathogenicity, am_class, score contribution, or None if not found
    """
    if not os.path.exists(db_path):
        return None

    # Normalise chromosome format — try both 'chr14' and '14'
    chrom_variants = [chrom, chrom.replace('chr', ''), f'chr{chrom.replace("chr", "")}']

    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()

        for chrom_v in chrom_variants:
            cursor.execute(
                """SELECT am_pathogenicity, am_class, protein_variant, uniprot_id
                   FROM alphamissense
                   WHERE chrom=? AND pos=? AND ref=? AND alt=?
                   LIMIT 1""",
                (chrom_v, pos, ref, alt)
            )
            row = cursor.fetchone()
            if row:
                conn.close()
                am_score, am_class, protein_var, uniprot = row
                return _format_result(chrom, pos, ref, alt, am_score, am_class, protein_var)

        conn.close()
        return None  # Variant not in AlphaMissense (e.g., non-missense, non-coding)

    except sqlite3.Error:
        return None


def _format_result(
    chrom: str, pos: int, ref: str, alt: str,
    am_score: float, am_class: str, protein_var: str
) -> Dict:
    """Format AlphaMissense result into signal dict."""
    
    # Normalise class label from TSV
    am_class_norm = am_class.lower().replace(' ', '_')
    if 'pathogenic' in am_class_norm:
        am_class_label = 'likely_pathogenic'
    elif 'benign' in am_class_norm:
        am_class_label = 'likely_benign'
    else:
        am_class_label = 'ambiguous'

    # Signal score: pathogenic contribution normalised to [0, 1]
    # Score = (am_score - benign_threshold) / (1 - benign_threshold)
    # Maps [0.340, 1.0] → [0, 1]; values below threshold → negative → clamp to 0
    normalized_score = max(0.0, (am_score - BENIGN_THRESHOLD) / (1 - BENIGN_THRESHOLD))
    normalized_score = min(1.0, normalized_score)

    if am_class_label == 'likely_pathogenic':
        signal_strength = 'HIGH'
    elif am_class_label == 'ambiguous':
        signal_strength = 'MEDIUM'
    else:
        signal_strength = 'LOW'

    return {
        "signal_type":      "alphamissense",
        "variant_id":       f"{chrom}:{pos}:{ref}>{alt}",
        "protein_change":   protein_var,
        "am_pathogenicity": round(am_score, 4),
        "am_class":         am_class_label,
        "score":            round(normalized_score, 3),
        "signal_strength":  signal_strength,
        "source":           "AlphaMissense_hg38",
        "reasoning": (
            f"AlphaMissense scores {chrom}:{pos}:{ref}>{alt} "
            f"({protein_var}) as {am_class_label} "
            f"(pathogenicity={am_score:.3f}). "
            f"Threshold: ≥0.564 = likely pathogenic."
        ),
    }


def check_variants_alphamissense(
    variants: List[Dict],
    db_path: str = DEFAULT_DB_PATH
) -> List[Dict]:
    """
    Check a list of variants against AlphaMissense.
    
    Each variant dict should have: chrom, pos, ref, alt
    (and optionally: gene, classification, hgvs)
    
    Returns list of AlphaMissense signals for variants that:
      - Are found in the database
      - Are classified as likely_pathogenic or ambiguous
    """
    signals = []
    for var in (variants or []):
        chrom = str(var.get('chrom', var.get('chromosome', '')))
        pos   = var.get('pos', var.get('position'))
        ref   = var.get('ref', '')
        alt   = var.get('alt', '')

        if not all([chrom, pos, ref, alt]):
            continue

        result = lookup_alphamissense(chrom, int(pos), ref, alt, db_path)
        if result and result['am_class'] != 'likely_benign':
            # Attach gene info if provided
            if 'gene' in var:
                result['gene'] = var['gene']
            signals.append(result)

    return signals