# tests/test_scoring_engine.py
"""
Tests for scoring_engine.py
Run from project root: python tests/test_scoring_engine.py
"""

import sys
sys.path.insert(0, "src")

from scoring_engine import hh_probability, site_stability_score, titration_curve, interpret_score

print("=" * 55)
print("TEST SUITE: scoring_engine.py")
print("=" * 55)

# ── Test 1: hh_probability formula ──────────────────────
print("\n[1] hh_probability — formula correctness")

assert abs(hh_probability(6.0, 6.0) - 0.5) < 0.001, "pH==pKa must give 0.5"
print("  ✓ hh_probability(6.0, 6.0) = 0.5")

assert hh_probability(8.0, 6.0) > 0.98, "pH >> pKa must give near 1.0"
print(f"  ✓ hh_probability(8.0, 6.0) = {hh_probability(8.0, 6.0):.4f} (>0.98)")

assert hh_probability(4.0, 6.0) < 0.02, "pH << pKa must give near 0.0"
print(f"  ✓ hh_probability(4.0, 6.0) = {hh_probability(4.0, 6.0):.4f} (<0.02)")

assert hh_probability(7.4, 6.0) > 0.95, "pH 7.4 vs pKa 6.0 should be >0.95"
print(f"  ✓ hh_probability(7.4, 6.0) = {hh_probability(7.4, 6.0):.4f} (>0.95)")

# Output must always be between 0 and 1
for ph, pka in [(0.0, 14.0), (14.0, 0.0), (7.0, 7.0)]:
    val = hh_probability(ph, pka)
    assert 0.0 <= val <= 1.0, f"Result must be in [0,1], got {val}"
print("  ✓ Output always in range [0.0, 1.0]")

# ── Test 2: site_stability_score ────────────────────────
print("\n[2] site_stability_score — joint probability")

# Single residue — score should equal individual probability
residues_1 = [("HIS", 6.0, 3)]
result = site_stability_score(residues_1, pH=7.4)
expected = hh_probability(7.4, 6.0)
assert abs(result["overall_score"] - expected) < 0.001
print(f"  ✓ Single residue score matches hh_probability: {result['overall_score']:.4f}")

# Three His site — joint probability
residues_3his = [("HIS", 6.0, 3), ("HIS", 6.0, 3), ("HIS", 6.0, 3)]
result_74 = site_stability_score(residues_3his, pH=7.4)
result_40 = site_stability_score(residues_3his, pH=4.0)

assert result_74["overall_score"] > 0.7, "3-His site at pH 7.4 should be stable"
assert result_40["overall_score"] < 0.1, "3-His site at pH 4.0 should be disrupted"
assert result_74["overall_score"] > result_40["overall_score"], "Higher pH = higher score"
print(f"  ✓ 3-His at pH 7.4: score={result_74['overall_score']:.4f} (>0.7)")
print(f"  ✓ 3-His at pH 4.0: score={result_40['overall_score']:.4f} (<0.1)")
print(f"  ✓ Score decreases as pH decreases")

# Confidence tier — should be worst (highest) tier among residues
mixed_tiers = [("HIS", 6.0, 1), ("CYS", 8.3, 3), ("ASP", 3.9, 2)]
result_mixed = site_stability_score(mixed_tiers, pH=7.0)
assert result_mixed["confidence_tier"] == 3, "Confidence tier should be worst (3)"
print(f"  ✓ Confidence tier = worst tier in site: {result_mixed['confidence_tier']}")

# Empty site
result_empty = site_stability_score([], pH=7.4)
assert result_empty["overall_score"] == 0.0
print(f"  ✓ Empty site returns score=0.0")

# ── Test 3: titration_curve ──────────────────────────────
print("\n[3] titration_curve — pH sweep")

residues = [("HIS", 6.0, 3), ("HIS", 6.54, 2)]
curve = titration_curve(residues)

assert len(curve) == 13, f"Default curve should have 13 points (3.0–9.0), got {len(curve)}"
print(f"  ✓ Default curve has {len(curve)} points (pH 3.0 to 9.0)")

# Scores must increase with pH for these residues
scores = [p["overall_score"] for p in curve]
assert scores == sorted(scores), "Scores must increase monotonically with pH"
print("  ✓ Scores increase monotonically with pH")

# Custom pH range
custom_curve = titration_curve(residues, pH_range=[6.5, 7.0, 7.4])
assert len(custom_curve) == 3
print(f"  ✓ Custom pH range works: {[p['pH'] for p in custom_curve]}")

# ── Test 4: interpret_score ──────────────────────────────
print("\n[4] interpret_score — labels")

assert interpret_score(0.90) == "High stability"
assert interpret_score(0.75) == "High stability"
assert interpret_score(0.50) == "Moderate stability"
assert interpret_score(0.40) == "Moderate stability"
assert interpret_score(0.20) == "Low stability"
assert interpret_score(0.09) == "Disrupted"
assert interpret_score(0.00) == "Disrupted"
print("  ✓ All score labels correct")

print("\n" + "=" * 55)
print("ALL TESTS PASSED ✓")
print("=" * 55)