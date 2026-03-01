"""
DxReanalyze Validation Script.
Tests the API with hardcoded clinical cases.
"""
import requests
import json
import sys

API_BASE_URL = "http://localhost:8000"


def test_case(name: str, payload: dict) -> dict:
    """Run a test case against the API."""
    print(f"\n{'='*60}")
    print(f"CASE: {name}")
    print('='*60)

    try:
        response = requests.post(
            f"{API_BASE_URL}/api/reanalysis-score",
            json=payload,
            timeout=30
        )
        response.raise_for_status()
        return response.json()
    except requests.exceptions.ConnectionError:
        print("ERROR: Cannot connect to API. Is the server running?")
        print(f"Make sure the API is running at {API_BASE_URL}")
        return None
    except Exception as e:
        print(f"ERROR: {e}")
        return None


def print_results(result: dict, expected_confidence: str) -> bool:
    """Print test results and check if passed."""
    if not result:
        print("RESULT: FAIL - No response from API")
        return False

    score = result.get('reanalysis_score', 0)
    confidence = result.get('confidence', 'UNKNOWN')

    print(f"\nReanalysis Score: {score}/100")
    print(f"Confidence: {confidence}")

    # Signal breakdown
    print("\nSignal Breakdown:")
    breakdown = result.get('signal_breakdown', {})

    vus = breakdown.get('vus_reclassification', {})
    print(f"  VUS Reclassification: {vus.get('score', 0):.3f} (weight: {vus.get('weight', 0)}) - {vus.get('signal_strength', 'N/A')}")
    if vus.get('signals'):
        for sig in vus['signals']:
            print(f"    - {sig.get('gene', '?')}: {sig.get('old_classification', '?')} → {sig.get('new_classification', '?')}")

    drift = breakdown.get('phenotypic_drift', {})
    drift_details = drift.get('details', {})
    print(f"  Phenotypic Drift: {drift.get('score', 0):.3f} (weight: {drift.get('weight', 0)}) - {drift_details.get('signal_strength', 'N/A')}")
    print(f"    New symptoms: {drift_details.get('new_symptom_count', 0)}")
    if drift_details.get('new_symptoms'):
        for sym in drift_details['new_symptoms'][:5]:  # Show first 5
            print(f"      - {sym.get('hpo_id', '?')}: {sym.get('name', 'Unknown')}")

    time_sig = breakdown.get('time_since_test', {})
    time_details = time_sig.get('details', {})
    print(f"  Time Since Test: {time_sig.get('score', 0):.3f} (weight: {time_sig.get('weight', 0)}) - {time_details.get('signal_strength', 'N/A')}")
    print(f"    Years elapsed: {time_details.get('years_since_test', 0):.1f}")

    # Clinical recommendation
    rec = result.get('clinical_recommendation', {})
    if rec:
        print("\nTop Reasons for Reanalysis:")
        for i, reason in enumerate(rec.get('top_reasons', []), 1):
            print(f"  {i}. {reason}")

        print("\nChecklist:")
        for item in rec.get('checklist', []):
            print(f"  [ ] {item}")

        print(f"\nNarrative: {rec.get('narrative', 'N/A')}")

    # Check pass/fail
    passed = confidence == expected_confidence
    status = "PASS" if passed else "FAIL"
    print(f"\n{'='*60}")
    print(f"RESULT: {status}")
    print(f"  Expected confidence: {expected_confidence}")
    print(f"  Actual confidence: {confidence}")
    print('='*60)

    return passed


def main():
    print("\n" + "="*60)
    print("DxReanalyze Validation Suite")
    print("="*60)

    # Check API health first
    try:
        health = requests.get(f"{API_BASE_URL}/api/health", timeout=5)
        health_data = health.json()
        print(f"\nAPI Status: {health_data.get('status', 'unknown')}")
        print(f"Datasets Loaded: {health_data.get('datasets_loaded', False)}")
    except requests.exceptions.ConnectionError:
        print("\nERROR: Cannot connect to API server.")
        print(f"Please start the API server first:")
        print(f"  cd dxreanalyze")
        print(f"  python -m uvicorn backend.main:app --reload --port 8000")
        sys.exit(1)

    # Test Case 1: FOXG1 Syndrome
    case1 = {
        "original_hpo_terms": ["HP:0001252", "HP:0001270", "HP:0011968"],
        "current_hpo_terms": [
            "HP:0001252", "HP:0001270", "HP:0011968",
            "HP:0001263", "HP:0002186", "HP:0000729",
            "HP:0002072", "HP:0001260"
        ],
        "prior_test": {
            "test_type": "gene panel",
            "test_date": "2019-06-01",
            "result": "negative",
            "vus_genes": ["FOXG1"]
        }
    }

    result1 = test_case("FOXG1 Syndrome (Child diagnosed late)", case1)
    passed1 = print_results(result1, "HIGH")

    # Test Case 2: Nina Nazar - Jansen's Disease
    case2 = {
        "original_hpo_terms": ["HP:0003521", "HP:0002750", "HP:0000843"],
        "current_hpo_terms": [
            "HP:0003521", "HP:0002750", "HP:0000843",
            "HP:0000518", "HP:0002748", "HP:0003075"
        ],
        "prior_test": {
            "test_type": "exome",
            "test_date": "2015-01-01",
            "result": "negative",
            "vus_genes": ["PTHR1"]
        }
    }

    result2 = test_case("Nina Nazar - Jansen's Metaphyseal Chondrodysplasia", case2)
    passed2 = print_results(result2, "HIGH")

    # Summary
    print("\n" + "="*60)
    print("VALIDATION SUMMARY")
    print("="*60)
    print(f"Case 1 (FOXG1): {'PASS' if passed1 else 'FAIL'}")
    print(f"Case 2 (Nina Nazar): {'PASS' if passed2 else 'FAIL'}")
    print(f"\nOverall: {2 if passed1 and passed2 else (1 if passed1 or passed2 else 0)}/2 cases passed")

    if passed1 and passed2:
        print("\nAll validation cases PASSED!")
        sys.exit(0)
    else:
        print("\nSome validation cases FAILED.")
        sys.exit(1)


if __name__ == "__main__":
    main()
