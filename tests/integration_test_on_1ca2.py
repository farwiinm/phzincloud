# src/integration_test_1ca2.py
"""
Full integration test on Carbonic Anhydrase II (1CA2).

Known biology:
  - Zinc coordinated by His94, His96, His119 (and a water molecule)
  - All three histidines have pKa ~6.0 (solution value)
  - At pH 7.4 (blood): site should be highly stable (score near 1.0)
  - At pH 4.0 (lysosome): site should be largely disrupted (score near 0.0)

This test validates that your scoring engine produces biologically
sensible results before you build the full pipeline.
"""

import sys
sys.path.insert(0, "src")

from parse_zinc_sites import parse_zinc_sites
from pka_lookup import get_pka
from scoring_engine import site_stability_score, titration_curve, interpret_score

PDB_FILE = "data/raw/test_proteins/1CA2.pdb"
PDB_ID   = "1CA2"

print("=" * 60)
print("INTEGRATION TEST: 1CA2 (Carbonic Anhydrase II)")
print("=" * 60)

# ── Step 1: Parse zinc sites ─────────────────────────────
print("\n[1] Parsing zinc coordination sites...")
sites = parse_zinc_sites(PDB_FILE, cutoff=5.0)
print(f"    Found {len(sites)} coordinating residues")

for s in sites:
    chain   = s.get("residue_chain") or s.get("chain")
    res_num = s.get("residue_seq_num")
    res_name = s.get("residue_name")
    dist    = s.get("distance_to_zinc") or s.get("distance_from_zinc")
    print(f"    {res_name}{res_num} chain={chain}  distance={float(dist):.2f} Å")

# ── Step 2: pKa lookup for each residue ─────────────────
print("\n[2] Looking up pKa values (three-tier system)...")
site_residues = []  # Will hold (residue_name, pKa, tier) tuples

for s in sites:
    chain    = s.get("residue_chain") or s.get("chain")
    res_num  = s.get("residue_seq_num")
    res_name = s.get("residue_name")

    pka_result = get_pka(
        pdb_id       = PDB_ID,
        chain        = chain,
        residue_num  = res_num,
        residue_name = res_name,
        pdb_file     = PDB_FILE
    )

    site_residues.append((res_name, pka_result["pka"], pka_result["tier"]))

    print(
        f"    {res_name}{res_num}: "
        f"pKa={pka_result['pka']}  "
        f"Tier {pka_result['tier']}  "
        f"({pka_result['source']})"
    )

# ── Step 3: Score at key biological pH values ────────────
print("\n[3] Stability scores at key biological pH values:")
print(f"    {'pH':<8} {'Score':<10} {'Label':<22} {'Weakest Residue'}")
print(f"    {'─'*6:<8} {'─'*8:<10} {'─'*20:<22} {'─'*20}")

key_pH_values = [4.0, 4.5, 5.0, 6.0, 6.5, 7.0, 7.4, 8.0, 9.0]

for pH in key_pH_values:
    result = site_stability_score(site_residues, pH)
    label  = interpret_score(result["overall_score"])
    print(
        f"    {pH:<8.1f} "
        f"{result['overall_score']:<10.4f} "
        f"{label:<22} "
        f"{result['weakest_residue']}"
    )

# ── Step 4: Biological sense check ──────────────────────
print("\n[4] Biological sense checks...")

score_74 = site_stability_score(site_residues, 7.4)["overall_score"]
score_40 = site_stability_score(site_residues, 4.0)["overall_score"]
score_60 = site_stability_score(site_residues, 6.0)["overall_score"]

print(f"    pH 7.4 (blood)     → score = {score_74:.4f}", end="  ")
print("✓ HIGH (as expected)" if score_74 > 0.7 else "✗ UNEXPECTED — check pKa values")

print(f"    pH 4.0 (lysosome)  → score = {score_40:.4f}", end="  ")
print("✓ LOW (as expected)" if score_40 < 0.1 else "✗ UNEXPECTED — check formula")

print(f"    pH 6.0 (pKa point) → score = {score_60:.4f}", end="  ")
print("✓ MID (~0.125 for 3-His site)" if 0.05 < score_60 < 0.3 else f"  (note: {score_60:.4f})")

# ── Step 5: Full titration curve ────────────────────────
print("\n[5] Full titration curve (pH 3.0 → 9.0):")
print(f"    {'pH':<8} {'Score':<10} {'Bar'}")
print(f"    {'─'*6:<8} {'─'*8:<10}")

curve = titration_curve(site_residues)
for point in curve:
    bar = "█" * int(point["overall_score"] * 30)
    print(f"    {point['pH']:<8.1f} {point['overall_score']:<10.4f} {bar}")

# ── Step 6: Per-residue breakdown at pH 7.4 ─────────────
print("\n[6] Per-residue breakdown at pH 7.4:")
result_74 = site_stability_score(site_residues, 7.4)
for r in result_74["per_residue_scores"]:
    print(
        f"    {r['residue_name']:<5} "
        f"pKa={r['pka']:<6} "
        f"P(deprotonated)={r['probability']:.4f}  "
        f"Tier {r['tier']}"
    )
print(f"\n    Joint score (product): {result_74['overall_score']:.6f}")
print(f"    Weakest residue:       {result_74['weakest_residue']}")
print(f"    Confidence tier:       {result_74['confidence_tier']}")

print("\n" + "=" * 60)
print("Integration test complete.")
print("=" * 60)