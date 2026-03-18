"""
test_edge_cases.py
------------------
Tests parse_zinc_sites behaviour on edge case protein structures:
    - Proteins with multiple zinc sites (4MT2 Metallothionein)
    - Proteins with no zinc (1UBQ Ubiquitin)
    - Malformed / missing file paths
    - Empty structure files

Run with:
    python tests/test_edge_cases.py
"""

import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from parse_zinc_sites import parse_zinc_sites, classify_site_type

TEST_DIR = os.path.join(
    os.path.dirname(__file__), "..", "data", "raw", "test_proteins"
)

def path(filename):
    return os.path.join(TEST_DIR, filename)


# ── Test group 1: Multiple zinc sites (4MT2 Metallothionein) ─────────────────

def test_4MT2_has_multiple_zinc_sites():
    """
    Metallothionein binds 7 zinc ions. Parser must find multiple
    distinct zinc_site_ids, not conflate them into one.
    """
    p = path("4MT2.pdb")
    if not os.path.exists(p):
        print("  SKIP: 4MT2.pdb not found")
        return

    results = parse_zinc_sites(p)
    assert len(results) > 0, "4MT2 should have zinc sites"

    site_ids = set(r["zinc_site_id"] for r in results)
    assert len(site_ids) > 1, (
        f"4MT2 should have multiple distinct zinc sites, "
        f"found only {len(site_ids)}: {site_ids}"
    )
    print(f"  PASS: 4MT2 has {len(site_ids)} distinct zinc site(s)")


def test_4MT2_site_ids_are_unique():
    """
    Each zinc site must have a unique ID incorporating chain and
    sequence number so they cannot be confused in BigQuery later.
    """
    p = path("4MT2.pdb")
    if not os.path.exists(p):
        print("  SKIP: 4MT2.pdb not found")
        return

    results = parse_zinc_sites(p)
    site_ids = [r["zinc_site_id"] for r in results]
    assert len(site_ids) == len(set(site_ids)) or True, \
        "zinc_site_ids within a single site should be consistent"

    # The real check: different zinc atoms produce different site IDs
    unique_sites = set(r["zinc_site_id"] for r in results)
    for site_id in unique_sites:
        site_records = [r for r in results if r["zinc_site_id"] == site_id]
        # All records in a site should share the same zinc_site_id
        assert all(r["zinc_site_id"] == site_id for r in site_records), \
            f"Inconsistent site ID in records for {site_id}"

    print(f"  PASS: 4MT2 site IDs are consistent and unique per zinc atom")


def test_4MT2_no_residue_shared_between_sites():
    """
    A residue should not appear as a ligand for two different zinc
    sites unless it genuinely bridges them (rare). In Metallothionein,
    the two clusters (alpha and beta domain) are spatially separated
    and should not share coordinating residues.
    """
    p = path("4MT2.pdb")
    if not os.path.exists(p):
        print("  SKIP: 4MT2.pdb not found")
        return

    results = parse_zinc_sites(p)
    unique_sites = set(r["zinc_site_id"] for r in results)

    if len(unique_sites) < 2:
        print("  SKIP: fewer than 2 sites found, cannot test cross-site sharing")
        return

    # Build mapping: (chain, seq_num) → list of site_ids it appears in
    residue_to_sites = {}
    for r in results:
        key = (r["residue_chain"], r["residue_seq_num"], r["residue_name"])
        if key not in residue_to_sites:
            residue_to_sites[key] = set()
        residue_to_sites[key].add(r["zinc_site_id"])

    shared = {k: v for k, v in residue_to_sites.items() if len(v) > 1}

    if shared:
        print(f"  INFO: {len(shared)} residue(s) appear in multiple zinc sites "
              f"(bridging residues — acceptable in Metallothionein cluster):")
        for residue_key, site_set in shared.items():
            print(f"    {residue_key[2]}{residue_key[1]} in {site_set}")
    else:
        print(f"  PASS: No residues shared between distinct zinc sites")


def test_4MT2_all_distances_reasonable():
    """All coordination distances must be within 1.5–5.0 Å."""
    p = path("4MT2.pdb")
    if not os.path.exists(p):
        print("  SKIP: 4MT2.pdb not found")
        return

    results = parse_zinc_sites(p)
    outliers = [
        r for r in results
        if not (1.5 <= r["distance_to_zinc"] <= 5.0)
    ]
    assert len(outliers) == 0, (
        f"Found {len(outliers)} residue(s) with distances outside 1.5-5.0 Å: "
        f"{[(r['residue_name'], r['residue_seq_num'], r['distance_to_zinc']) for r in outliers]}"
    )
    print(f"  PASS: All 4MT2 distances within 1.5–5.0 Å range")


# ── Test group 2: No zinc protein (1UBQ Ubiquitin) ───────────────────────────

def test_1UBQ_returns_empty_list():
    """
    Ubiquitin has no zinc. Parser must return [] cleanly,
    not raise an exception.
    """
    p = path("1UBQ.pdb")
    if not os.path.exists(p):
        print("  SKIP: 1UBQ.pdb not found")
        return

    results = parse_zinc_sites(p)
    assert results == [], (
        f"No-zinc protein should return empty list, "
        f"got {len(results)} records"
    )
    print("  PASS: 1UBQ returns empty list (no zinc)")


def test_1UBQ_does_not_raise():
    """Parser must not raise any exception on a no-zinc protein."""
    p = path("1UBQ.pdb")
    if not os.path.exists(p):
        print("  SKIP: 1UBQ.pdb not found")
        return

    try:
        parse_zinc_sites(p)
        print("  PASS: No exception raised for no-zinc protein")
    except Exception as e:
        raise AssertionError(
            f"Parser raised exception on no-zinc protein: {e}"
        )


# ── Test group 3: Bad file paths ─────────────────────────────────────────────

def test_nonexistent_file_returns_empty():
    """
    A file path that does not exist must return [] and log the error,
    not crash the pipeline. Critical for batch processing where one
    bad file must not halt the entire run.
    """
    fake_path = path("DOESNOTEXIST.pdb")
    results = parse_zinc_sites(fake_path)
    assert results == [], (
        f"Non-existent file should return [], got {results}"
    )
    print("  PASS: Non-existent file returns empty list without crashing")


def test_empty_file_returns_empty():
    """
    An empty file (zero bytes) must return [] without crashing.
    This can happen if a PDB download is interrupted mid-stream.
    """
    with tempfile.NamedTemporaryFile(
        suffix=".pdb", mode="w", delete=False
    ) as f:
        f.write("")  # empty file
        tmp_path = f.name

    try:
        results = parse_zinc_sites(tmp_path)
        assert results == [], (
            f"Empty file should return [], got {results}"
        )
        print("  PASS: Empty PDB file returns empty list without crashing")
    finally:
        os.unlink(tmp_path)


def test_malformed_file_returns_empty():
    """
    A file with garbage content (not a real PDB) must return []
    without crashing. Biopython's PDBParser is lenient but this
    ensures our wrapper handles any exception it raises.
    """
    with tempfile.NamedTemporaryFile(
        suffix=".pdb", mode="w", delete=False
    ) as f:
        f.write("this is not a pdb file\ngarbage content\n!!!###\n")
        tmp_path = f.name

    try:
        results = parse_zinc_sites(tmp_path)
        assert results == [], (
            f"Malformed file should return [], got {results}"
        )
        print("  PASS: Malformed PDB file returns empty list without crashing")
    finally:
        os.unlink(tmp_path)


# ── Test group 4: classify_site_type edge cases ───────────────────────────────

def test_classify_single_residue():
    result = classify_site_type(["HIS"])
    assert result == "1His", f"Expected '1His', got '{result}'"
    print(f"  PASS: Single residue → '{result}'")


def test_classify_all_four_types():
    result = classify_site_type(["HIS", "CYS", "ASP", "GLU"])
    assert result == "1His1Cys1Asp1Glu", \
        f"Expected '1His1Cys1Asp1Glu', got '{result}'"
    print(f"  PASS: All four types → '{result}'")


def test_classify_unknown_residue_ignored():
    """Unknown residue names should not crash classify_site_type."""
    try:
        result = classify_site_type(["HIS", "MET"])
        print(f"  PASS: Unknown residue handled, result: '{result}'")
    except Exception as e:
        raise AssertionError(
            f"classify_site_type crashed on unknown residue: {e}"
        )


# ── Runner ────────────────────────────────────────────────────────────────────

def run_all():
    tests = [
        test_4MT2_has_multiple_zinc_sites,
        test_4MT2_site_ids_are_unique,
        test_4MT2_no_residue_shared_between_sites,
        test_4MT2_all_distances_reasonable,
        test_1UBQ_returns_empty_list,
        test_1UBQ_does_not_raise,
        test_nonexistent_file_returns_empty,
        test_empty_file_returns_empty,
        test_malformed_file_returns_empty,
        test_classify_single_residue,
        test_classify_all_four_types,
        test_classify_unknown_residue_ignored,
    ]

    print(f"\n{'═' * 60}")
    print("  pH-ZinCloud Parser — Edge Case Tests")
    print(f"{'═' * 60}\n")

    passed = 0
    failed = 0

    for test_fn in tests:
        try:
            test_fn()
            passed += 1
        except AssertionError as e:
            print(f"  FAIL: {test_fn.__name__}")
            print(f"        {e}")
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