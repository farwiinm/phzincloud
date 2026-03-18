# src/pipeline.py
"""
pH-ZinCloud Master Pipeline
Chains all modules together:
    PDB file → parse zinc sites → pKa lookup → HH scoring → results CSV

This is the core of the artifact. Every downstream component
(Docker, Cloud Run, BigQuery, Streamlit) uses this as its engine.

Usage:
    python src/pipeline.py
    Or import and call run_pipeline() from other scripts.
"""
import sys
import os
import pandas as pd
from typing import List, Dict, Any

# ── Path setup ────────────────────────────────────────────────────────────────
# Works correctly whether called as:
#   python src/pipeline.py          (from project root)
#   python pipeline.py              (from inside src/)
#   CMD ["python", "src/pipeline.py"] (from Docker /app)
_SRC_DIR = os.path.dirname(os.path.abspath(__file__))
if _SRC_DIR not in sys.path:
    sys.path.insert(0, _SRC_DIR)

from parse_zinc_sites import parse_zinc_sites
from pka_lookup import get_pka
from scoring_engine import (
    site_stability_score,
    titration_curve,
    interpret_score
)

# ── Configuration from environment variables ──────────────────────────────────
# These are set by Docker/Cloud Run at runtime.
# Sensible local defaults are provided so the script works without Docker too.
PDB_INPUT_DIR       = os.environ.get("PDB_INPUT_DIR",  "data/raw/test_proteins")
OUTPUT_CSV          = os.environ.get("OUTPUT_CSV",      "results/results.csv")
CUTOFF              = float(os.environ.get("CUTOFF",   "5.0"))
LOG_LEVEL           = os.environ.get("LOG_LEVEL",       "INFO")

# The pH values the pipeline scores every site at
DEFAULT_PH_VALUES   = [4.0, 5.0, 6.0, 7.0, 7.4, 8.0, 9.0]

# pH switch detection thresholds
PH_SWITCH_THRESHOLD = 0.5


# ─────────────────────────────────────────────────────
# Core pipeline function for a single PDB file
# ─────────────────────────────────────────────────────
def run_single(
    pdb_file: str,
    pH_values: List[float] = None
) -> List[Dict[str, Any]]:
    """
    Run the full pipeline on a single PDB file.

    Steps:
        1. Parse zinc coordination sites (Biopython)
        2. Look up pKa for each coordinating residue (3-tier system)
        3. Group residues by zinc site
        4. Score each site at every pH value (Henderson-Hasselbalch)
        5. Flag pH-sensitive sites

    Args:
        pdb_file:  Path to the PDB file
        pH_values: List of pH values to score at (default: DEFAULT_PH_VALUES)

    Returns:
        List of row dicts — one row per residue per protein, with
        per-pH scores and site-level flags attached.
        Returns empty list if no zinc sites found or on error.
    """
    if pH_values is None:
        pH_values = DEFAULT_PH_VALUES

    pdb_id = os.path.basename(pdb_file).replace(".pdb", "").upper()
    rows   = []

    # ── Step 1: Parse ────────────────────────────────────
    try:
        sites = parse_zinc_sites(pdb_file, cutoff=CUTOFF)
    except Exception as e:
        print(f"  [ERROR] Parser failed on {pdb_id}: {e}")
        return []

    if not sites:
        print(f"  [INFO]  No zinc sites found in {pdb_id}")
        return []

    # ── Step 2: Group residues by zinc site ID ───────────
    site_groups: Dict[str, list]  = {}
    site_pka_data: Dict[str, list] = {}

    for s in sites:
        chain    = s.get("residue_chain") or s.get("chain")
        res_num  = s.get("residue_seq_num")
        res_name = s.get("residue_name")
        zinc_id  = s.get("zinc_site_id") or s.get("zinc_id", f"{pdb_id}_ZN")
        distance = float(s.get("distance_to_zinc") or s.get("distance_from_zinc", 0.0))

        # ── Step 3: pKa lookup ───────────────────────────
        pka_result = get_pka(
            pdb_id       = pdb_id,
            chain        = chain,
            residue_num  = res_num,
            residue_name = res_name,
            pdb_file     = pdb_file
        )

        if zinc_id not in site_groups:
            site_groups[zinc_id]   = []
            site_pka_data[zinc_id] = []

        site_groups[zinc_id].append({
            "chain":      chain,
            "res_name":   res_name,
            "res_num":    res_num,
            "distance":   distance,
            "pka_value":  pka_result["pka"],
            "pka_tier":   pka_result["tier"],
            "pka_source": pka_result["source"],
        })
        site_pka_data[zinc_id].append(
            (res_name, pka_result["pka"], pka_result["tier"])
        )

    # ── Step 4: Score each site at every pH ─────────────
    for zinc_id, residues in site_groups.items():
        pka_tuples = site_pka_data[zinc_id]

        ph_scores = {}
        for pH in pH_values:
            result = site_stability_score(pka_tuples, pH)
            ph_scores[pH] = result["overall_score"]

        result_74  = site_stability_score(pka_tuples, 7.4)
        label_74   = interpret_score(result_74["overall_score"])
        conf_tier  = result_74["confidence_tier"]
        weakest    = result_74["weakest_residue"]

        # ── Step 5: pH-switch detection ──────────────────
        score_at_80 = ph_scores.get(8.0, site_stability_score(pka_tuples, 8.0)["overall_score"])
        score_at_70 = ph_scores.get(7.0, site_stability_score(pka_tuples, 7.0)["overall_score"])
        score_at_60 = ph_scores.get(6.0, site_stability_score(pka_tuples, 6.0)["overall_score"])

        drop_80_to_60 = score_at_80 - score_at_60
        is_ph_switch  = bool(
            score_at_80   >  0.5  and
            score_at_60   <  0.5  and
            drop_80_to_60 >  0.4  and
            score_at_70   >  0.3
        )

        # ── Build one row per residue ────────────────────
        for res in residues:
            row = {
                "pdb_id":           pdb_id,
                "zinc_site_id":     zinc_id,
                "residue_name":     res["res_name"],
                "residue_seq":      res["res_num"],
                "chain":            res["chain"],
                "distance_to_zinc": res["distance"],
                "pka_value":        res["pka_value"],
                "pka_tier":         res["pka_tier"],
                "pka_source":       res["pka_source"],
            }
            for pH in pH_values:
                col = f"pH_{str(pH).replace('.', '_')}_score"
                row[col] = round(ph_scores[pH], 6)

            row["site_stability_7_4"] = round(result_74["overall_score"], 6)
            row["stability_label"]    = label_74
            row["weakest_residue"]    = weakest
            row["confidence_tier"]    = conf_tier
            row["is_ph_switch"]       = is_ph_switch
            row["n_coordinating"]     = len(residues)
            rows.append(row)

    return rows


# ─────────────────────────────────────────────────────
# Batch runner
# ─────────────────────────────────────────────────────
def run_pipeline(
    pdb_folder: str,
    output_csv: str,
    pH_values: List[float] = None
) -> pd.DataFrame:
    """
    Run the full pipeline on all PDB files in a folder or GCS bucket.
    Accepts both local paths and GCS URIs (gs://bucket/prefix/).
    """
    if pH_values is None:
        pH_values = DEFAULT_PH_VALUES

    # Import GCS utilities (handles both local and cloud paths)
    from gcs_utils import list_pdb_files, download_pdb, upload_results, is_gcs_path

    # Get list of PDB files — works for both local and GCS paths
    pdb_files = list_pdb_files(pdb_folder)

    if not pdb_files:
        print(f"ERROR: No PDB files found at {pdb_folder}")
        return pd.DataFrame()

    print("=" * 60)
    print("pH-ZinCloud — Full Pipeline Run")
    print("=" * 60)
    print(f"Input        : {pdb_folder}")
    print(f"Output CSV   : {output_csv}")
    print(f"pH values    : {pH_values}")
    print(f"PDB files    : {len(pdb_files)} found")
    print("=" * 60)

    all_rows = []

    # Temp directory for GCS downloads
    import tempfile
    tmp_dir = tempfile.mkdtemp()

    for pdb_path in pdb_files:
        # If GCS path, download to temp dir first
        if is_gcs_path(pdb_path):
            try:
                local_pdb = download_pdb(pdb_path, tmp_dir)
            except Exception as e:
                print(f"  [ERROR] Could not download {pdb_path}: {e}")
                continue
        else:
            local_pdb = pdb_path

        pdb_id = os.path.basename(local_pdb).replace(".pdb", "").upper()
        print(f"\nProcessing: {pdb_id}")

        rows = run_single(local_pdb, pH_values)
        if rows:
            print(f"  → {len(rows)} residue rows produced")
            all_rows.extend(rows)
        else:
            print(f"  → No output (no zinc sites or parser error)")

    if not all_rows:
        print("\nERROR: No results produced across all proteins.")
        return pd.DataFrame()

    df = pd.DataFrame(all_rows)

    # Save locally first
    local_output = output_csv if not is_gcs_path(output_csv) else os.path.join(tmp_dir, "results.csv")
    output_dir = os.path.dirname(os.path.abspath(local_output))
    os.makedirs(output_dir, exist_ok=True)
    df.to_csv(local_output, index=False)

    # Upload to GCS if output path is a GCS URI
    if is_gcs_path(output_csv):
        upload_results(local_output, output_csv)

    print(f"\n{'=' * 60}")
    print(f"PIPELINE COMPLETE")
    print(f"  Total rows      : {len(df)}")
    print(f"  Unique proteins : {df['pdb_id'].nunique()}")
    print(f"  Output saved    : {output_csv}")

    return df

# ─────────────────────────────────────────────────────
# Results summary printer
# ─────────────────────────────────────────────────────
def print_results_summary(df: pd.DataFrame):
    """Print a human-readable summary of pipeline results."""
    if df.empty:
        print("No results to summarise.")
        return

    print("\n" + "=" * 60)
    print("RESULTS SUMMARY")
    print("=" * 60)

    print("\nPER-PROTEIN STABILITY AT pH 7.4:")
    print(f"  {'PDB':<8} {'Site':<25} {'Score@7.4':<12} {'Label':<22} {'pH-Switch':<10} {'Tier'}")
    print(f"  {'─'*6:<8} {'─'*23:<25} {'─'*10:<12} {'─'*20:<22} {'─'*8:<10} {'─'*4}")

    site_cols = ["pdb_id", "zinc_site_id", "site_stability_7_4",
                 "stability_label", "is_ph_switch", "confidence_tier"]
    sites_df = df[site_cols].drop_duplicates(subset=["zinc_site_id"])

    for _, row in sites_df.iterrows():
        switch_flag = "⚠ YES" if row["is_ph_switch"] else "no"
        print(
            f"  {row['pdb_id']:<8} "
            f"{str(row['zinc_site_id']):<25} "
            f"{row['site_stability_7_4']:<12.4f} "
            f"{row['stability_label']:<22} "
            f"{switch_flag:<10} "
            f"Tier {row['confidence_tier']}"
        )

    n_switch = sites_df["is_ph_switch"].sum()
    n_total  = len(sites_df)
    print(f"\n  pH-switch sites detected: {n_switch}/{n_total}")

    print("\nTIER COVERAGE:")
    tier_counts = df["pka_tier"].value_counts().sort_index()
    total = len(df)
    for tier, count in tier_counts.items():
        pct   = count / total * 100
        label = {1: "Experimental (PKAD-R)",
                 2: "Computational (PROPKA)",
                 3: "Canonical default"}.get(tier, "Unknown")
        bar   = "█" * int(pct / 4)
        print(f"  Tier {tier} — {label:<25}: {count:>3}/{total} ({pct:5.1f}%)  {bar}")

    print("\nMEAN SITE SCORE BY pH (averaged across all sites):")
    ph_cols    = [c for c in df.columns if c.startswith("pH_") and c.endswith("_score")]
    site_means = df.groupby("zinc_site_id")[ph_cols].mean()
    overall    = site_means.mean()
    for col in ph_cols:
        ph_val     = col.replace("pH_", "").replace("_score", "").replace("_", ".")
        mean_score = overall[col]
        bar        = "█" * int(mean_score * 20)
        print(f"  pH {ph_val:<5}: {mean_score:.4f}  {bar}")

    print("=" * 60)
    print("\n>>> THESIS NOTE: The pH-switch count and tier coverage above")
    print(">>> go directly into Chapter 6 Results. Record them now.\n")


# ─────────────────────────────────────────────────────
# Entry point — reads paths from environment variables
# ─────────────────────────────────────────────────────
if __name__ == "__main__":
    df = run_pipeline(
        pdb_folder = PDB_INPUT_DIR,
        output_csv = OUTPUT_CSV,
        pH_values  = DEFAULT_PH_VALUES
    )
    if not df.empty:
        print_results_summary(df)
        print("\nFIRST 10 ROWS OF OUTPUT:")
        pd.set_option("display.max_columns", 8)
        pd.set_option("display.width", 120)
        print(df.head(10).to_string(index=False))