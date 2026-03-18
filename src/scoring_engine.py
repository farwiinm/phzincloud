# src/scoring_engine.py
"""
Henderson-Hasselbalch Scoring Engine for pH-ZinCloud

Computes pH-dependent zinc coordination stability scores based on the
protonation states of coordinating residues.

Core concept:
    At a given pH, each coordinating residue has a probability of being
    deprotonated (available to bind zinc). The site stability score is
    the joint probability that ALL residues are simultaneously deprotonated.

Usage:
    from scoring_engine import hh_probability, site_stability_score, titration_curve
"""

from typing import List, Tuple
import math


# ─────────────────────────────────────────────────────
# FUNCTION 1: Henderson-Hasselbalch probability
# ─────────────────────────────────────────────────────

def hh_probability(pH: float, pKa: float) -> float:
    """
    Compute the probability that a residue is deprotonated at a given pH.

    Uses the Henderson-Hasselbalch equation:
        P(deprotonated) = 1 / (1 + 10^(pKa - pH))

    Args:
        pH:  Environmental pH (e.g. 7.4 for blood, 6.5 for tumour, 4.5 for lysosome)
        pKa: The pKa of the titratable residue

    Returns:
        Float between 0.0 and 1.0.
        0.0 = fully protonated (cannot coordinate zinc)
        1.0 = fully deprotonated (available to coordinate zinc)

    Examples:
        hh_probability(6.0, 6.0) → 0.5    (pH == pKa, exactly 50/50)
        hh_probability(8.0, 6.0) → ~0.990  (pH >> pKa, mostly deprotonated)
        hh_probability(4.0, 6.0) → ~0.010  (pH << pKa, mostly protonated)
    """
    return 1.0 / (1.0 + 10 ** (pKa - pH))


# ─────────────────────────────────────────────────────
# FUNCTION 2: Site stability score
# ─────────────────────────────────────────────────────

def site_stability_score(
    site_residues: List[Tuple[str, float, int]],
    pH: float
) -> dict:
    """
    Compute the stability score for a zinc coordination site at a given pH.

    The overall score is the joint probability that ALL coordinating residues
    are simultaneously deprotonated (i.e. all available to coordinate zinc).

    Args:
        site_residues: List of tuples → (residue_name, pKa_value, tier)
                       e.g. [("HIS", 6.0, 3), ("HIS", 6.54, 2), ("GLU", 4.1, 3)]
        pH:            The pH to score at

    Returns:
        Dict with keys:
            overall_score        → joint probability product (0.0 – 1.0)
            per_residue_scores   → list of (residue_name, pKa, prob) per residue
            minimum_residue_score → the weakest link (lowest individual probability)
            weakest_residue      → name of the weakest residue
            confidence_tier      → worst (highest number) tier among all residues
            n_residues           → number of coordinating residues scored
    """
    if not site_residues:
        return {
            "overall_score": 0.0,
            "per_residue_scores": [],
            "minimum_residue_score": 0.0,
            "weakest_residue": None,
            "confidence_tier": None,
            "n_residues": 0
        }

    per_residue = []
    overall_score = 1.0
    worst_tier = 1  # Will be updated to the highest (least confident) tier

    for residue_name, pka_value, tier in site_residues:

        # Skip residues with no pKa (shouldn't happen, but safety net)
        if pka_value is None:
            continue

        prob = hh_probability(pH, pka_value)
        overall_score *= prob

        per_residue.append({
            "residue_name": residue_name,
            "pka":          pka_value,
            "tier":         tier,
            "probability":  round(prob, 4)
        })

        if tier is not None and tier > worst_tier:
            worst_tier = tier

    if not per_residue:
        return {
            "overall_score": 0.0,
            "per_residue_scores": [],
            "minimum_residue_score": 0.0,
            "weakest_residue": None,
            "confidence_tier": worst_tier,
            "n_residues": 0
        }

    # Find the weakest link (residue with lowest deprotonation probability)
    weakest = min(per_residue, key=lambda x: x["probability"])

    return {
        "overall_score":         round(overall_score, 6),
        "per_residue_scores":    per_residue,
        "minimum_residue_score": weakest["probability"],
        "weakest_residue":       f"{weakest['residue_name']} (pKa={weakest['pka']})",
        "confidence_tier":       worst_tier,
        "n_residues":            len(per_residue)
    }


# ─────────────────────────────────────────────────────
# FUNCTION 3: Titration curve
# ─────────────────────────────────────────────────────

def titration_curve(
    site_residues: List[Tuple[str, float, int]],
    pH_range: List[float] = None
) -> List[dict]:
    """
    Compute stability scores across a pH range to produce a titration curve.

    Args:
        site_residues: List of (residue_name, pKa, tier) tuples
        pH_range:      List of pH values to score at.
                       Defaults to pH 3.0 – 9.0 in 0.5 steps if not provided.

    Returns:
        List of dicts, one per pH step, each containing:
            pH, overall_score, minimum_residue_score, confidence_tier
    """
    if pH_range is None:
        # Default: full physiological range in 0.5 steps
        pH_range = [round(x * 0.5, 1) for x in range(6, 19)]  # 3.0 to 9.0

    curve = []
    for pH in pH_range:
        result = site_stability_score(site_residues, pH)
        curve.append({
            "pH":                   pH,
            "overall_score":        result["overall_score"],
            "minimum_residue_score": result["minimum_residue_score"],
            "confidence_tier":      result["confidence_tier"]
        })

    return curve


# ─────────────────────────────────────────────────────
# FUNCTION 4: Score interpretation label
# ─────────────────────────────────────────────────────

def interpret_score(score: float) -> str:
    """
    Convert a numeric stability score into a human-readable label.
    Used in the dashboard and results table.

    Score ranges based on joint probability thresholds:
        >= 0.75  → High stability   (site likely intact)
        >= 0.40  → Moderate         (partial coordination likely)
        >= 0.10  → Low              (significant proton competition)
        <  0.10  → Disrupted        (zinc likely displaced)
    """
    if score >= 0.75:
        return "High stability"
    elif score >= 0.40:
        return "Moderate stability"
    elif score >= 0.10:
        return "Low stability"
    else:
        return "Disrupted"
    
henderson_hasselbalch = hh_probability
