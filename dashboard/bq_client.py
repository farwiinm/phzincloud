# dashboard/bq_client.py
"""
BigQuery client for pH-ZinCloud dashboard.
Handles all database queries with graceful fallback to local CSV
when BigQuery credentials are not available (local development mode).
"""
import os
import pandas as pd
import logging

logger = logging.getLogger(__name__)

# ── Configuration ─────────────────────────────────────────────────────────────
PROJECT_ID  = os.environ.get("GCP_PROJECT",    "phzincloud")
DATASET_ID  = os.environ.get("BQ_DATASET",     "phzincloud_results")
TABLE_ID    = os.environ.get("BQ_TABLE",        "residue_scores")
LOCAL_CSV   = os.environ.get("LOCAL_CSV",       "results/results.csv")

FULL_TABLE  = f"{PROJECT_ID}.{DATASET_ID}.{TABLE_ID}"


def _get_client():
    """Return a BigQuery client, or None if credentials unavailable."""
    try:
        from google.cloud import bigquery
        return bigquery.Client(project=PROJECT_ID)
    except Exception as e:
        logger.warning(f"BigQuery client unavailable: {e}")
        return None


def load_local_csv() -> pd.DataFrame:
    """Load results from local CSV as fallback."""
    if os.path.exists(LOCAL_CSV):
        return pd.read_csv(LOCAL_CSV)
    return pd.DataFrame()


def get_all_proteins() -> list:
    """Return list of all PDB IDs available in the dataset."""
    client = _get_client()
    if client:
        try:
            query = f"SELECT DISTINCT pdb_id FROM `{FULL_TABLE}` ORDER BY pdb_id"
            df = client.query(query).to_dataframe()
            return df["pdb_id"].tolist()
        except Exception as e:
            logger.warning(f"BigQuery query failed: {e}")

    # Fallback to local CSV
    df = load_local_csv()
    if not df.empty:
        return sorted(df["pdb_id"].unique().tolist())
    return []


def get_protein_data(pdb_id: str) -> pd.DataFrame:
    """
    Fetch all residue rows for a given PDB ID.
    Returns DataFrame with all columns from residue_scores table.
    """
    client = _get_client()
    if client:
        try:
            query = f"""
                SELECT * FROM `{FULL_TABLE}`
                WHERE pdb_id = @pdb_id
                ORDER BY zinc_site_id, residue_seq
            """
            job_config = client.query.__func__.__doc__  # just to import
            from google.cloud import bigquery
            job_config = bigquery.QueryJobConfig(
                query_parameters=[
                    bigquery.ScalarQueryParameter("pdb_id", "STRING", pdb_id.upper())
                ]
            )
            df = client.query(query, job_config=job_config).to_dataframe()
            if not df.empty:
                return df
        except Exception as e:
            logger.warning(f"BigQuery query failed for {pdb_id}: {e}")

    # Fallback to local CSV
    df = load_local_csv()
    if not df.empty:
        return df[df["pdb_id"] == pdb_id.upper()].copy()
    return pd.DataFrame()


def get_site_scores(pdb_id: str, zinc_site_id: str) -> pd.DataFrame:
    """
    Get per-pH stability scores for a specific zinc site.
    Returns one row with all pH score columns.
    """
    df = get_protein_data(pdb_id)
    if df.empty:
        return pd.DataFrame()
    site_df = df[df["zinc_site_id"] == zinc_site_id]
    return site_df.drop_duplicates(subset=["zinc_site_id"])


def run_pipeline_for_pdb(pdb_id: str) -> pd.DataFrame:
    """
    If a PDB ID is not in the database, run the pipeline on it directly.
    Downloads from RCSB, runs full pipeline, returns results DataFrame.
    """
    import sys
    import tempfile
    import requests

    sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
    from pipeline import run_single

    # Download PDB file from RCSB
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
    except Exception as e:
        logger.error(f"Failed to download {pdb_id} from RCSB: {e}")
        return pd.DataFrame()

    # Save to temp file and run pipeline
    with tempfile.NamedTemporaryFile(
        suffix=".pdb", mode="w", delete=False
    ) as f:
        f.write(response.text)
        tmp_path = f.name

    try:
        rows = run_single(tmp_path)
        if rows:
            return pd.DataFrame(rows)
        return pd.DataFrame()
    finally:
        os.unlink(tmp_path)