# src/run_pka_assignment.py
"""
Integration script: runs the zinc site parser on all test proteins,
then calls get_pka() on each coordinating residue.

Produces: results/zinc_sites_pka_enriched.csv

Run from project root:
    python src/run_pka_assignment.py
"""

import sys
import os
import pandas as pd

sys.path.insert(0, "src")

from parse_zinc_sites import parse_zinc_sites
from pka_lookup import get_pka

# ─────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────
PDB_FOLDER  = "data/raw/test_proteins"
OUTPUT_FILE = "results/zinc_sites_pka_enriched.csv"


def process_all_proteins(pdb_folder: str) -> pd.DataFrame:
    all_rows = []

    pdb_files = sorted([f for f in os.listdir(pdb_folder) if f.endswith(".pdb")])

    if not pdb_files:
        print(f"ERROR: No PDB files found in {pdb_folder}")
        return pd.DataFrame()

    print(f"Found {len(pdb_files)} PDB files: {pdb_files}\n")

    for filename in pdb_files:
        pdb_path = os.path.join(pdb_folder, filename)
        pdb_id   = filename.replace(".pdb", "").upper()

        print(f"{'─'*50}")
        print(f"Processing: {pdb_id}")

        # ── Step 1: Parse zinc coordination sites ──
        try:
            sites = parse_zinc_sites(pdb_path, cutoff=5.0)
        except Exception as e:
            print(f"  ERROR — parser failed: {e}")
            continue

        if not sites:
            print(f"  No zinc sites found — skipping")
            continue

        print(f"  {len(sites)} coordinating residue(s) found")

        # ── Step 2: Assign pKa to each residue ──────
        for site in sites:

            # ── Key fix: read the actual field names your parser uses ──
            chain       = site.get("residue_chain") or site.get("chain")
            res_num     = site.get("residue_seq_num")
            res_name    = site.get("residue_name")
            zinc_id     = site.get("zinc_site_id") or site.get("zinc_id", "ZN")
            distance    = site.get("distance_to_zinc") or site.get("distance_from_zinc", 0.0)

            # Safety check — skip if any critical field is missing
            if not all([chain, res_num, res_name]):
                print(f"  WARNING: Incomplete site data — {site}")
                continue

            pka_result = get_pka(
                pdb_id       = pdb_id,
                chain        = chain,
                residue_num  = res_num,
                residue_name = res_name,
                pdb_file     = pdb_path
            )

            row = {
                "pdb_id":            pdb_id,
                "zinc_id":           zinc_id,
                "chain":             chain,
                "residue_name":      res_name,
                "residue_seq_num":   res_num,
                "distance_to_zinc":  round(float(distance), 3),
                "pka_value":         pka_result["pka"],
                "pka_tier":          pka_result["tier"],
                "pka_source":        pka_result["source"],
            }
            all_rows.append(row)

            print(
                f"  {res_name:<4}{res_num:<6} chain={chain}  "
                f"pKa={str(pka_result['pka']):<6}  "
                f"Tier {pka_result['tier']}  "
                f"({pka_result['source']})"
            )

    return pd.DataFrame(all_rows)


def print_tier_coverage(df: pd.DataFrame):
    total = len(df)
    print("\n" + "=" * 55)
    print("TIER COVERAGE SUMMARY")
    print("=" * 55)

    for tier in [1, 2, 3]:
        count = len(df[df["pka_tier"] == tier])
        pct   = (count / total * 100) if total > 0 else 0
        label = {
            1: "Tier 1 — Experimental  (PKAD-R) ",
            2: "Tier 2 — Computational (PROPKA) ",
            3: "Tier 3 — Canonical default       "
        }[tier]
        bar = "█" * int(pct / 5)
        print(f"  {label}: {count:>3}/{total}  ({pct:5.1f}%)  {bar}")

    print(f"\n  Unique proteins processed : {df['pdb_id'].nunique()}")
    print(f"  Proteins with zinc sites  : {df['pdb_id'].unique().tolist()}")

    print("\nPER-PROTEIN TIER BREAKDOWN:")
    summary = df.groupby(["pdb_id", "pka_tier"]).size().unstack(fill_value=0)
    summary.columns = [f"Tier {c}" for c in summary.columns]
    print(summary.to_string())
    print("=" * 55)

    


def main():
    print("=" * 55)
    print("pH-ZinCloud — pKa Assignment Pipeline")
    print("=" * 55)

    os.makedirs("results", exist_ok=True)

    df = process_all_proteins(PDB_FOLDER)

    if df.empty:
        print("\nERROR: No results produced.")
        print("Check that your PDB files contain zinc (ZN) records.")
        return

    # Save enriched CSV
    df.to_csv(OUTPUT_FILE, index=False)
    print(f"\n✓ Enriched results saved → {OUTPUT_FILE}")
    print(f"  Rows: {len(df)}  |  Columns: {list(df.columns)}")

    # Tier coverage summary
    print_tier_coverage(df)

    # Preview
    print("\nFIRST 10 ROWS:")
    print(df.head(10).to_string(index=False))


if __name__ == "__main__":
    main()