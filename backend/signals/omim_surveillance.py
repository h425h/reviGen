"""
OMIM New Gene-Disease Association Surveillance — NEW SIGNAL.

This is the single highest-impact signal for rare disease reanalysis.
Radboud 2022 study: 42% of new diagnoses came from newly characterised
gene-disease associations that didn't exist at time of original testing.

OMIM adds 250-300 new gene-disease associations annually. A gene that
returned 'no disease association' in 2018 may now have 40+ published cases.

API docs: https://www.omim.org/api
Free key for research: https://www.omim.org/api#apikey

phenotypeMappingKey values (per OMIM spec):
  1 = nondiseases
  2 = chromosome deletion or duplication syndrome
  3 = molecular basis of disorder known      ← THIS IS THE ONE WE WANT
  4 = mendelian inheritance
"""

import os
import requests
from datetime import datetime
from typing import List, Dict, Optional


OMIM_API_BASE = "https://api.omim.org/api"

# Fallback: hardcoded known associations for demo/offline use
# Maps gene → (omim_disease_id, disease_name, year_characterized)
KNOWN_ASSOCIATIONS_FALLBACK = {
    "FOXG1":  ("OMIM:613454", "FOXG1 syndrome", "2008"),
    "CDKL5":  ("OMIM:300672", "CDKL5 deficiency disorder", "2012"),
    "SETBP1": ("OMIM:616078", "Schinzel-Giedion midface retraction syndrome", "2010"),
    "ANKRD11":("OMIM:148050", "KBG syndrome", "2011"),
    "PTHR1":  ("OMIM:156400", "Jansen metaphyseal chondrodysplasia", "1994"),
    "STXBP1": ("OMIM:612164", "Developmental and epileptic encephalopathy 4", "2008"),
    "KCNQ2":  ("OMIM:613720", "Developmental and epileptic encephalopathy 7", "2012"),
    "SCN2A":  ("OMIM:613721", "Developmental and epileptic encephalopathy 11", "2009"),
    "GRIN2A": ("OMIM:245570", "Epilepsy-aphasia spectrum", "2010"),
    "KCNT1":  ("OMIM:615005", "Epileptic encephalopathy, early infantile, 14", "2012"),
    "ADGRV1": ("OMIM:276900", "Usher syndrome type IIC", "2002"),
    "DYRK1A": ("OMIM:614104", "Intellectual disability, autosomal dominant 7", "2011"),
    "MED13L": ("OMIM:616278", "MED13L syndrome", "2012"),
    "SETD5":  ("OMIM:615761", "Mental retardation, autosomal dominant 23", "2014"),
    "KAT6A":  ("OMIM:616268", "Mental retardation, autosomal dominant 32", "2015"),
}


def check_omim_gene(
    gene_symbol: str,
    patient_test_date: str,
    omim_api_key: Optional[str] = None
) -> Optional[Dict]:
    """
    Check OMIM for new gene-disease associations after the patient's test date.
    
    Tries live OMIM API first; falls back to hardcoded table for demo.
    
    Args:
        gene_symbol:       Gene to check (e.g., 'FOXG1')
        patient_test_date: ISO date string 'YYYY-MM-DD'
        omim_api_key:      OMIM API key (from env OMIM_API_KEY)
    
    Returns:
        Signal dict if association found, None otherwise
    """
    gene_upper = gene_symbol.upper().strip()
    
    try:
        test_dt = datetime.strptime(patient_test_date, "%Y-%m-%d")
    except (ValueError, TypeError):
        return None

    api_key = omim_api_key or os.environ.get("OMIM_API_KEY", "")
    
    if api_key:
        return _check_omim_live(gene_upper, test_dt, api_key)
    else:
        return _check_omim_fallback(gene_upper, test_dt)


def _check_omim_live(gene: str, test_dt: datetime, api_key: str) -> Optional[Dict]:
    """
    Hit the OMIM API for real-time gene-disease association data.
    Endpoint: /api/entry/search?search={gene}&include=geneMap
    """
    try:
        resp = requests.get(
            f"{OMIM_API_BASE}/entry/search",
            params={
                "search": gene,
                "include": "geneMap",
                "format": "json",
                "apiKey": api_key,
            },
            timeout=10,
        )
        resp.raise_for_status()
        data = resp.json()
    except Exception as e:
        # On API failure, fall back to hardcoded table
        return _check_omim_fallback(gene, test_dt)

    entries = (
        data.get("omim", {})
        .get("searchResponse", {})
        .get("entryList", [])
    )

    for entry_wrapper in entries:
        entry = entry_wrapper.get("entry", {})
        gene_map = entry.get("geneMap", {})
        
        # phenotypeMappingKey=3 → molecular basis of disorder confirmed
        phenotype_list = gene_map.get("phenotypeMapList", [])
        
        for ph_wrapper in phenotype_list:
            ph = ph_wrapper.get("phenotypeMap", {})
            
            if ph.get("phenotypeMappingKey") != 3:
                continue  # Skip unconfirmed associations
            
            disease_name = ph.get("phenotype", "Unknown disease")
            date_updated = gene_map.get("dateUpdated", "")
            omim_number  = entry.get("mimNumber", "")
            
            # Did this association appear/get confirmed after the test?
            is_new = False
            if date_updated:
                try:
                    updated_dt = datetime.strptime(date_updated, "%Y-%m-%d")
                    is_new = updated_dt > test_dt
                except ValueError:
                    pass
            
            # Score: new associations get 0.85, existing confirmed get 0.55
            # (existing still matters — gene confirmed pathogenic post-test)
            score = 0.85 if is_new else 0.55
            strength = "HIGH" if is_new else "MEDIUM"
            
            return {
                "signal_type":    "omim_gene_disease",
                "gene":           gene,
                "disease":        disease_name,
                "omim_id":        f"OMIM:{omim_number}" if omim_number else "Unknown",
                "date_updated":   date_updated,
                "is_new_after_test": is_new,
                "score":          score,
                "signal_strength": strength,
                "source":         "OMIM_API",
                "reasoning": (
                    f"OMIM confirms molecular basis for {disease_name} in {gene}. "
                    f"{'Association updated AFTER test date — major reanalysis trigger.' if is_new else 'Association existed at test time but gene-disease link is confirmed.'}"
                ),
            }

    return None  # Gene has no confirmed disease association in OMIM


def _check_omim_fallback(gene: str, test_dt: datetime) -> Optional[Dict]:
    """
    Fallback using hardcoded table when no OMIM API key is available.
    Used for demo/hackathon without live API access.
    """
    if gene not in KNOWN_ASSOCIATIONS_FALLBACK:
        return None
    
    omim_id, disease_name, year_str = KNOWN_ASSOCIATIONS_FALLBACK[gene]
    
    try:
        assoc_year = int(year_str)
        assoc_dt   = datetime(assoc_year, 6, 1)  # Approximate mid-year
        is_new     = assoc_dt > test_dt
    except (ValueError, TypeError):
        is_new = False
    
    score    = 0.85 if is_new else 0.55
    strength = "HIGH" if is_new else "MEDIUM"
    
    return {
        "signal_type":       "omim_gene_disease",
        "gene":              gene,
        "disease":           disease_name,
        "omim_id":           omim_id,
        "date_updated":      year_str,
        "is_new_after_test": is_new,
        "score":             score,
        "signal_strength":   strength,
        "source":            "OMIM_FALLBACK",
        "reasoning": (
            f"Gene {gene} is causative for {disease_name} ({omim_id}). "
            f"{'Disease characterised AFTER test date ({year_str}) — reanalysis strongly indicated.'.format(year_str=year_str) if is_new else f'Confirmed gene-disease association ({year_str}). VUS in this gene warrants expert review.'}"
        ),
    }


def check_omim_all_genes(
    vus_genes: List[str],
    patient_test_date: str,
    omim_api_key: Optional[str] = None
) -> List[Dict]:
    """
    Check all VUS genes for OMIM signals.
    Returns list of signals, one per gene with a confirmed disease.
    """
    signals = []
    for gene in (vus_genes or []):
        result = check_omim_gene(gene, patient_test_date, omim_api_key)
        if result:
            signals.append(result)
    return signals