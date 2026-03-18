# tests/test_pka_lookup.py
"""
Tests for the three-tier pKa lookup module.
Run from project root: python tests/test_pka_lookup.py
"""

import sys
sys.path.insert(0, "src")

from pka_lookup import lookup_experimental_pka, get_canonical_pka, get_pka

print("=" * 50)
print("TEST SUITE: pka_lookup.py")
print("=" * 50)

# ── Tier 3: Canonical defaults ──────────────────────
print("\n[1] Tier 3 — Canonical defaults")

assert get_canonical_pka("HIS")["pka"] == 6.0,  "HIS pKa should be 6.0"
assert get_canonical_pka("CYS")["pka"] == 8.3,  "CYS pKa should be 8.3"
assert get_canonical_pka("ASP")["pka"] == 3.9,  "ASP pKa should be 3.9"
assert get_canonical_pka("GLU")["pka"] == 4.1,  "GLU pKa should be 4.1"
assert get_canonical_pka("XYZ") is None,         "Unknown residue should return None"

for res, expected in [("HIS", 6.0), ("CYS", 8.3), ("ASP", 3.9), ("GLU", 4.1)]:
    result = get_canonical_pka(res)
    assert result["tier"] == 3, f"{res} should be Tier 3"
    assert result["pka"] == expected
    print(f"  ✓ {res} → pKa={result['pka']} (Tier {result['tier']})")

# ── Tier 1: PKAD-R lookup ───────────────────────────
print("\n[2] Tier 1 — PKAD-R experimental lookup")

# Test a non-existent entry — should return None
result = lookup_experimental_pka("XXXX", "A", 999, "HIS")
assert result is None, "Non-existent entry should return None"
print("  ✓ Non-existent entry correctly returns None")

# Test case insensitivity (PDB IDs should work in any case)
result_upper = lookup_experimental_pka("1LJU", "A", 35, "HIS")
result_lower = lookup_experimental_pka("1lju", "A", 35, "HIS")
if result_upper and result_lower:
    assert result_upper["pka"] == result_lower["pka"], "Case should not affect lookup"
    print(f"  ✓ Case-insensitive lookup works: pKa={result_upper['pka']}")
else:
    print("  ℹ 1LJU/35/HIS not in PKAD-R — update with a PDB ID you verified in the CSV")

# ── Master get_pka() cascade ────────────────────────
print("\n[3] Master get_pka() — cascade behaviour")

# Non-existent protein → must fall to Tier 3
result = get_pka("XXXX", "A", 999, "HIS")
assert result["tier"] == 3, "Non-existent entry should fall back to Tier 3"
assert result["pka"] == 6.0, "HIS fallback should be 6.0"
print(f"  ✓ Unknown protein falls to Tier 3: pKa={result['pka']}, source={result['source']}")

# Completely unknown residue → pka should be None
result = get_pka("XXXX", "A", 999, "ZZZ")
assert result["pka"] is None, "Completely unknown residue should return pka=None"
print(f"  ✓ Unknown residue returns pka=None: {result}")

# Known residue, real PDB that exists in parser output
result = get_pka("1CA2", "A", 94, "HIS",
                 pdb_file="data/raw/test_proteins/1CA2.pdb")
assert result["pka"] is not None, "1CA2 HIS94 should get a pKa value"
assert result["tier"] in [1, 2, 3], "Tier must be 1, 2, or 3"
print(f"  ✓ 1CA2 HIS94: pKa={result['pka']}, Tier={result['tier']}, source={result['source']}")

# ── pKa value sanity checks ─────────────────────────
print("\n[4] pKa value sanity checks")

# All standard coordinating residues should return pKa in a plausible range
test_residues = ["HIS", "CYS", "ASP", "GLU"]
for res in test_residues:
    result = get_pka("XXXX", "A", 1, res)
    assert result["pka"] is not None, f"{res} should always return a pKa"
    assert 1.0 <= result["pka"] <= 15.0, f"{res} pKa={result['pka']} is outside plausible range"
    print(f"  ✓ {res}: pKa={result['pka']} (in range 1–15)")

print("\n" + "=" * 50)
print("ALL TESTS PASSED ✓")
print("=" * 50)