"""
Data loaders for DxReanalyze.
Loads HPO ontology, HPOA phenotype annotations, and ClinVar variant data.
"""
import os
import pandas as pd
import networkx as nx
from typing import Dict, List, Tuple, Optional

# Module-level globals to store loaded data
hpo_terms: Dict[str, Dict] = {}
hpo_graph: Optional[nx.DiGraph] = None
disease_phenotypes: Dict[str, List[str]] = {}
clinvar_df: Optional[pd.DataFrame] = None

# Target genes for filtering
TARGET_GENES = ['FOXG1', 'CDKL5', 'SETBP1', 'ANKRD11', 'PTHR1']


def parse_hpo_obo(filepath: str) -> Tuple[Dict[str, Dict], nx.DiGraph]:
    """
    Parse hp.obo into dict {hpo_id: {name, parents: []}} and networkx DiGraph.
    """
    global hpo_terms, hpo_graph

    terms = {}
    graph = nx.DiGraph()

    current_term = None
    current_id = None

    print(f"Loading HPO ontology from {filepath}...")

    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()

            if line == '[Term]':
                # Save previous term if exists
                if current_id and current_term:
                    terms[current_id] = current_term
                    graph.add_node(current_id, **current_term)
                current_term = {'name': '', 'parents': []}
                current_id = None
            elif line.startswith('id: HP:'):
                current_id = line.split('id: ')[1].strip()
            elif line.startswith('name: ') and current_term is not None:
                current_term['name'] = line.split('name: ')[1].strip()
            elif line.startswith('is_a: HP:') and current_term is not None:
                parent_id = line.split('is_a: ')[1].split(' !')[0].strip()
                current_term['parents'].append(parent_id)
                if current_id:
                    graph.add_edge(parent_id, current_id)

    # Don't forget the last term
    if current_id and current_term:
        terms[current_id] = current_term
        graph.add_node(current_id, **current_term)

    print(f"Loaded {len(terms)} HPO terms")

    hpo_terms = terms
    hpo_graph = graph

    return terms, graph


def parse_hpoa(filepath: str) -> Dict[str, List[str]]:
    """
    Parse phenotype.hpoa into dict {disease_id: [hpo_terms]}.
    Skips comment lines starting with #.
    """
    global disease_phenotypes

    mappings = {}

    print(f"Loading HPOA from {filepath}...")

    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            # Skip comment lines
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue

            # HPOA format: database_id, disease_name, qualifier, hpo_id, ...
            disease_id = parts[0]
            hpo_id = parts[3] if len(parts) > 3 else None

            if disease_id and hpo_id and hpo_id.startswith('HP:'):
                if disease_id not in mappings:
                    mappings[disease_id] = []
                if hpo_id not in mappings[disease_id]:
                    mappings[disease_id].append(hpo_id)

    print(f"Loaded {len(mappings)} disease-phenotype mappings")

    disease_phenotypes = mappings
    return mappings


def load_clinvar(filepath: str) -> pd.DataFrame:
    """
    Load variant_summary.txt in chunks of 50,000 rows.
    Keep only rows where GeneSymbol is in target genes and
    ClinicalSignificance contains relevant classifications.
    """
    global clinvar_df

    print(f"Loading ClinVar from {filepath}...")
    print("This may take a moment due to file size...")

    # Columns to keep
    keep_cols = ['GeneSymbol', 'ClinicalSignificance', 'LastEvaluated',
                 'ReviewStatus', 'PhenotypeList', 'Type']

    # Read in chunks
    chunk_size = 50000
    chunks = []
    chunk_num = 0

    try:
        for chunk in pd.read_csv(filepath, sep='\t', chunksize=chunk_size,
                                  low_memory=False, on_bad_lines='skip'):
            chunk_num += 1

            if chunk_num % 10 == 0:
                print(f"Processing chunk {chunk_num}...")

            # Filter by gene symbol
            if 'GeneSymbol' in chunk.columns:
                gene_mask = chunk['GeneSymbol'].isin(TARGET_GENES)
                filtered = chunk[gene_mask]

                # Filter by clinical significance
                if 'ClinicalSignificance' in filtered.columns and len(filtered) > 0:
                    sig_mask = filtered['ClinicalSignificance'].str.contains(
                        'Pathogenic|Likely pathogenic|Uncertain significance',
                        case=False, na=False
                    )
                    filtered = filtered[sig_mask]

                # Keep only required columns
                available_cols = [c for c in keep_cols if c in filtered.columns]
                if len(filtered) > 0:
                    chunks.append(filtered[available_cols])

    except Exception as e:
        print(f"Error reading ClinVar file: {e}")
        return pd.DataFrame()

    if chunks:
        result = pd.concat(chunks, ignore_index=True)
    else:
        result = pd.DataFrame(columns=keep_cols)

    print(f"ClinVar filtered to {len(result)} rows")

    clinvar_df = result
    return result


def get_hpo_name(hpo_id: str) -> str:
    """Get human-readable name for an HPO term."""
    if hpo_id in hpo_terms:
        return hpo_terms[hpo_id].get('name', hpo_id)
    return hpo_id


def get_clinvar_data() -> pd.DataFrame:
    """Get the loaded ClinVar DataFrame."""
    global clinvar_df
    if clinvar_df is None:
        return pd.DataFrame()
    return clinvar_df


def load_all_datasets(data_dir: str) -> bool:
    """Load all datasets from the data directory."""
    try:
        hpo_path = os.path.join(data_dir, 'hp.obo')
        hpoa_path = os.path.join(data_dir, 'phenotype.hpoa')
        clinvar_path = os.path.join(data_dir, 'variant_summary.txt')

        # Check files exist
        for path in [hpo_path, hpoa_path, clinvar_path]:
            if not os.path.exists(path):
                print(f"Missing file: {path}")
                return False

        # Load datasets
        parse_hpo_obo(hpo_path)
        parse_hpoa(hpoa_path)
        load_clinvar(clinvar_path)

        return True
    except Exception as e:
        print(f"Error loading datasets: {e}")
        return False
