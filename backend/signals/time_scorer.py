"""
Time-based Signal Scorer.
Scores the urgency of reanalysis based on time since the original test.
"""
from datetime import datetime, date
from typing import Dict


def score_time_since_test(test_date: str) -> Dict:
    """
    Score the time elapsed since the original genetic test.

    Args:
        test_date: Date of original test (YYYY-MM-DD format)

    Returns:
        Dictionary with time-based scoring
    """
    # Parse test date
    try:
        test_datetime = datetime.strptime(test_date, '%Y-%m-%d')
    except (ValueError, TypeError):
        # Return default low score for invalid dates
        return {
            'time_score': 0.0,
            'years_since_test': 0.0,
            'signal_strength': 'LOW'
        }

    # Calculate years elapsed
    today = datetime.now()
    days_elapsed = (today - test_datetime).days
    years_elapsed = days_elapsed / 365.0

    # Calculate time score (maxes out at 5 years)
    time_score = min(1.0, years_elapsed / 5.0)

    # Determine signal strength
    if years_elapsed > 3:
        signal_strength = 'HIGH'
    elif years_elapsed > 1:
        signal_strength = 'MEDIUM'
    else:
        signal_strength = 'LOW'

    return {
        'time_score': round(time_score, 3),
        'years_since_test': round(years_elapsed, 2),
        'signal_strength': signal_strength
    }
