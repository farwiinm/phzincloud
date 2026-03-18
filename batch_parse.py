"""
batch_parse.py
--------------
Runs parse_zinc_sites over a folder of PDB files and saves
all results to a single CSV file: results/zinc_sites_raw.csv

Usage:
    python batch_parse.py

Expects PDB files in: data/raw/test_proteins/
Outputs CSV to:       results/zinc_sites_raw.csv
"""

import os
import csv
import sys

# Add src to path so we can import from it
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))

from parse_zinc_sites import parse_zinc_sites, summarise_site

# ── Config ────────────────────────────────────────────────────────────────────
PDB_FOLDER  = "data/raw/test_proteins"
OUTPUT_CSV  = "results/zinc_sites_raw.csv"
CUTOFF      = 5.0

# ── Ensure output directory exists ───────────────────────────────────────────
os.makedirs("results", exist_ok=True)
os.makedirs("logs",    exist_ok=True)

# ── CSV column headers ────────────────────────────────────────────────────────
FIELDNAMES = [
    "pdb_id", "zinc_site_id", "zinc_chain", "zinc_seq_num",
    "site_type", "residue_name", "residue_chain", "residue_seq_num",
    "coord_atom", "distance_to_zinc", "cb_coords", "total_ligands"
]


def run_batch(pdb_folder: str, output_csv: str):
    """
    Process all .pdb files in pdb_folder and write results to output_csv.
    """
    pdb_files = [
        os.path.join(pdb_folder, f)
        for f in os.listdir(pdb_folder)
        if f.lower().endswith(".pdb")
    ]

    if not pdb_files:
        print(f"No PDB files found in {pdb_folder}")
        return

    print(f"Found {len(pdb_files)} PDB files. Starting batch parse...\n")

    all_results   = []
    proteins_done = 0
    zinc_found    = 0
    no_zinc       = 0
    errors        = 0

    for pdb_file in sorted(pdb_files):
        pdb_name = os.path.basename(pdb_file)
        try:
            results = parse_zinc_sites(pdb_file, cutoff=CUTOFF)

            if results:
                all_results.extend(results)
                zinc_found += 1

                # Print a summary for this protein
                # Group by zinc_site_id for clean output
                site_ids = dict.fromkeys(r["zinc_site_id"] for r in results)
                for site_id in site_ids:
                    site_records = [r for r in results
                                    if r["zinc_site_id"] == site_id]
                    summary = summarise_site(site_records)
                    print(f"  ✓ {pdb_name}: {summary['zinc_site_id']} "
                          f"— {summary['site_type']} "
                          f"— {summary['residue_summary']}")
            else:
                no_zinc += 1
                print(f"  — {pdb_name}: No zinc sites found.")

        except Exception as e:
            errors += 1
            print(f"  ✗ {pdb_name}: ERROR — {e}")

        proteins_done += 1

    # Write CSV
    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(all_results)

    # Print final summary
    print(f"\n{'─' * 60}")
    print(f"Batch complete.")
    print(f"  Proteins processed : {proteins_done}")
    print(f"  With zinc          : {zinc_found}")
    print(f"  Without zinc       : {no_zinc}")
    print(f"  Errors             : {errors}")
    print(f"  Total residue rows : {len(all_results)}")
    print(f"  Output saved to    : {output_csv}")
    print(f"{'─' * 60}\n")


if __name__ == "__main__":
    run_batch(PDB_FOLDER, OUTPUT_CSV)