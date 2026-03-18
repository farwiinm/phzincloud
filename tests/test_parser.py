"""
test_parser.py
--------------
Basic assertions to verify parse_zinc_sites.py produces
correct output on known protein structures.

Run with:  python tests/test_parser.py

Expected results are taken from published coordination chemistry
literature (Andreini et al., 2011; Auld, 2001).
"""

import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from parse_zinc_sites import (
    parse_zinc_sites,
    classify_site_type,
    summarise_site
)

# ── Path to test proteins ──────────────────────────────────────────────────
TEST_DIR = os.path.join(
    os.path.dirname(__file__), "..", "data", "raw", "test_proteins"
)

def path(filename):
    return os.path.join(TEST_DIR, filename)


# ── Test 1: Carbonic Anhydrase II (1CA2) ──────────────────────────────────
# Known coordination: His94, His96, His119 (3His site)
# Source: Eriksson et al., 1988; Andreini et al., 2011

def test_1CA2_has_zinc():
    results = parse_zinc_sites(path("1CA2.pdb"))
    assert len(results) > 0, "1CA2 should contain at least one zinc site"
    print("  PASS: 1CA2 has zinc sites")

def test_1CA2_has_three_histidines():
    results = parse_zinc_sites(path("1CA2.pdb"))
    his_residues = [r for r in results if r["residue_name"] == "HIS"]
    assert len(his_residues) >= 3, (
        f"1CA2 should have at least 3 coordinating His, found {len(his_residues)}"
    )
    print(f"  PASS: 1CA2 has {len(his_residues)} coordinating His residue(s)")

def test_1CA2_site_type():
    results = parse_zinc_sites(path("1CA2.pdb"))
    site_ids = dict.fromkeys(r["zinc_site_id"] for r in results)
    for site_id in site_ids:
        site_records = [r for r in results if r["zinc_site_id"] == site_id]
        site_type = site_records[0]["site_type"]
        # Carbonic anhydrase has a 3His site (water is 4th ligand, not counted)
        assert "3His" in site_type or "His" in site_type, (
            f"1CA2 site type should contain His, got: {site_type}"
        )
    print(f"  PASS: 1CA2 site type contains His")

def test_1CA2_distances_reasonable():
    results = parse_zinc_sites(path("1CA2.pdb"))
    for r in results:
        assert 1.5 <= r["distance_to_zinc"] <= 5.0, (
            f"Distance out of range for {r['residue_name']}"
            f"{r['residue_seq_num']}: {r['distance_to_zinc']} Å"
        )
    print("  PASS: All distances in 1CA2 are within 1.5–5.0 Å range")

def test_1CA2_coordinating_atoms():
    results = parse_zinc_sites(path("1CA2.pdb"))
    his_results = [r for r in results if r["residue_name"] == "HIS"]
    for r in his_results:
        assert r["coord_atom"] in ["NE2", "ND1"], (
            f"His coordinating atom should be NE2 or ND1, got: {r['coord_atom']}"
        )
    print("  PASS: 1CA2 His residues coordinate via NE2 or ND1")


# ── Test 2: No-zinc protein ────────────────────────────────────────────────
# The parser must return an empty list, not crash.
# Download any small protein without zinc — e.g. 1UBQ (Ubiquitin)
# If you don't have it, this test uses a simple path check

def test_no_zinc_returns_empty_list():
    # We simulate a no-zinc result by calling with a known no-zinc file
    # If 1UBQ.pdb is not present, we skip gracefully
    no_zinc_path = path("1UBQ.pdb")
    if not os.path.exists(no_zinc_path):
        print("  SKIP: 1UBQ.pdb not found — download it to enable this test")
        return
    results = parse_zinc_sites(no_zinc_path)
    assert results == [], (
        f"No-zinc protein should return empty list, got {len(results)} rows"
    )
    print("  PASS: No-zinc protein returns empty list")


# ── Test 3: classify_site_type ────────────────────────────────────────────

def test_classify_site_type_3his1glu():
    result = classify_site_type(["HIS", "HIS", "HIS", "GLU"])
    assert result == "3His1Glu", f"Expected '3His1Glu', got '{result}'"
    print(f"  PASS: classify_site_type(['HIS','HIS','HIS','GLU']) = '{result}'")

def test_classify_site_type_2his2cys():
    result = classify_site_type(["HIS", "CYS", "HIS", "CYS"])
    assert result == "2His2Cys", f"Expected '2His2Cys', got '{result}'"
    print(f"  PASS: classify_site_type mixed = '{result}'")

def test_classify_site_type_empty():
    result = classify_site_type([])
    assert result == "Unknown", f"Expected 'Unknown' for empty list, got '{result}'"
    print(f"  PASS: classify_site_type([]) = '{result}'")


# ── Test 4: summarise_site ────────────────────────────────────────────────

def test_summarise_site_basic():
    fake_data = [
        {
            "zinc_site_id": "1CA2_ZN_263_A",
            "site_type": "3His",
            "total_ligands": 3,
            "residue_name": "HIS",
            "residue_seq_num": 94,
            "coord_atom": "NE2",
            "distance_to_zinc": 2.05
        },
        {
            "zinc_site_id": "1CA2_ZN_263_A",
            "site_type": "3His",
            "total_ligands": 3,
            "residue_name": "HIS",
            "residue_seq_num": 96,
            "coord_atom": "NE2",
            "distance_to_zinc": 2.10
        },
    ]
    summary = summarise_site(fake_data)
    assert summary["zinc_site_id"] == "1CA2_ZN_263_A"
    assert summary["avg_distance"] == round((2.05 + 2.10) / 2, 3)
    print(f"  PASS: summarise_site produces correct avg_distance")


# ── Run all tests ─────────────────────────────────────────────────────────

def run_all():
    tests = [
        test_1CA2_has_zinc,
        test_1CA2_has_three_histidines,
        test_1CA2_site_type,
        test_1CA2_distances_reasonable,
        test_1CA2_coordinating_atoms,
        test_no_zinc_returns_empty_list,
        test_classify_site_type_3his1glu,
        test_classify_site_type_2his2cys,
        test_classify_site_type_empty,
        test_summarise_site_basic,
    ]

    print(f"\n{'═' * 60}")
    print(f"  pH-ZinCloud Parser Tests")
    print(f"{'═' * 60}\n")

    passed = 0
    failed = 0

    for test_fn in tests:
        try:
            test_fn()
            passed += 1
        except AssertionError as e:
            print(f"  FAIL: {test_fn.__name__} — {e}")
            failed += 1
        except Exception as e:
            print(f"  ERROR: {test_fn.__name__} — {e}")
            failed += 1

    print(f"\n{'─' * 60}")
    print(f"  Results: {passed} passed, {failed} failed")
    print(f"{'─' * 60}\n")

    if failed > 0:
        sys.exit(1)


if __name__ == "__main__":
    run_all()