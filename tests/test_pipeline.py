# tests/test_pipeline.py
"""
Integration tests for pipeline.py
Run from project root: python tests/test_pipeline.py
"""

import sys
import os
import pandas as pd

sys.path.insert(0, "src")

from pipeline import run_single, run_pipeline, DEFAULT_PH_VALUES

print("=" * 55)
print("TEST SUITE: pipeline.py")
print("=" * 55)

# ── Test 1: run_single on 1CA2 ───────────────────────────
print("\n[1] run_single() on 1CA2")

PDB_1CA2 = "data/raw/test_proteins/1CA2.pdb"
rows = run_single(PDB_1CA2)

assert len(rows) > 0, "1CA2 should produce rows (has zinc)"
print(f"  ✓ {len(rows)} rows produced for 1CA2")

# Every row must have required columns
required_cols = [
    "pdb_id", "zinc_site_id", "residue_name", "residue_seq",
    "pka_value", "pka_tier", "site_stability_7_4",
    "is_ph_switch", "stability_label", "confidence_tier"
]
for col in required_cols:
    assert col in rows[0], f"Missing column: {col}"
print(f"  ✓ All required columns present")

# pdb_id must be 1CA2
assert all(r["pdb_id"] == "1CA2" for r in rows)
print(f"  ✓ pdb_id correctly set to 1CA2")

# Scores must be between 0 and 1
for pH in DEFAULT_PH_VALUES:
    col = f"pH_{str(pH).replace('.', '_')}_score"
    assert col in rows[0], f"Missing pH column: {col}"
    for r in rows:
        assert 0.0 <= r[col] <= 1.0, f"Score out of range at {col}: {r[col]}"
print(f"  ✓ All pH score columns present and in range [0, 1]")

# Biological sense: score at pH 7.4 > score at pH 4.0 for 1CA2
score_74 = rows[0]["pH_7_4_score"]
score_40 = rows[0]["pH_4_0_score"]
assert score_74 > score_40, "pH 7.4 score must be higher than pH 4.0 score"
print(f"  ✓ pH 7.4 score ({score_74:.4f}) > pH 4.0 score ({score_40:.4f})")

# ── Test 2: run_single on protein with no zinc ───────────
print("\n[2] run_single() on no-zinc protein (1UBQ)")

PDB_UBIQ = "data/raw/test_proteins/1UBQ.pdb"
if os.path.exists(PDB_UBIQ):
    rows_ubiq = run_single(PDB_UBIQ)
    assert rows_ubiq == [], "No-zinc protein should return empty list"
    print("  ✓ No-zinc protein returns empty list (no crash)")
else:
    print("  ℹ 1UBQ.pdb not found — skipping no-zinc test")

# ── Test 3: run_pipeline produces a valid CSV ────────────
print("\n[3] run_pipeline() — full batch run")

TEST_OUTPUT = "results/test_pipeline_output.csv"
df = run_pipeline(
    pdb_folder = "data/raw/test_proteins",
    output_csv = TEST_OUTPUT,
    pH_values  = [4.0, 7.4]   # Minimal pH range for speed
)

assert not df.empty, "Pipeline should produce a non-empty DataFrame"
assert os.path.exists(TEST_OUTPUT), "CSV file should be created"
print(f"  ✓ Pipeline produced {len(df)} rows")
print(f"  ✓ CSV saved to {TEST_OUTPUT}")

# All PDB IDs should be uppercase
assert all(df["pdb_id"] == df["pdb_id"].str.upper()), "PDB IDs should be uppercase"
print("  ✓ All PDB IDs are uppercase")

# No null pka_values (every residue must have some pKa)
assert df["pka_value"].notna().all(), "No residue should have null pKa"
print("  ✓ No null pKa values")

# Tier must be 1, 2, or 3
valid_tiers = {1, 2, 3}
assert set(df["pka_tier"].unique()).issubset(valid_tiers), "Tier must be 1, 2, or 3"
print(f"  ✓ Tiers found: {sorted(df['pka_tier'].unique())}")

# Clean up test output
os.remove(TEST_OUTPUT)
print(f"  ✓ Test output file cleaned up")

print("\n" + "=" * 55)
print("ALL TESTS PASSED ✓")
print("=" * 55)