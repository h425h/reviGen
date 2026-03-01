from backend.signals.nlp_extractor import extract_hpo_from_notes

def extract_hpo_with_metadata(clinical_text, api_key=None, reference_date=None):
    return extract_hpo_from_notes(clinical_text, api_key=api_key, reference_date=reference_date)
