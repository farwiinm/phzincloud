# src/gcs_utils.py
"""
Google Cloud Storage utilities for pH-ZinCloud.
Handles downloading PDB files from GCS and uploading results.
Falls back gracefully to local file operations when GCS is not configured,
allowing the same pipeline.py to run both locally and in Cloud Run.
"""
import os
import logging

# Only import GCS library if available — allows local runs without GCP credentials
try:
    from google.cloud import storage
    GCS_AVAILABLE = True
except ImportError:
    GCS_AVAILABLE = False


def is_gcs_path(path: str) -> bool:
    """Return True if the path is a GCS URI (gs://bucket/...)."""
    return path.startswith("gs://")


def parse_gcs_path(gcs_uri: str) -> tuple:
    """
    Parse a GCS URI into bucket name and blob path.
    e.g. 'gs://my-bucket/test_proteins/1CA2.pdb' → ('my-bucket', 'test_proteins/1CA2.pdb')
    """
    assert gcs_uri.startswith("gs://"), f"Not a GCS path: {gcs_uri}"
    without_prefix = gcs_uri[5:]
    bucket_name, _, blob_path = without_prefix.partition("/")
    return bucket_name, blob_path


def list_pdb_files(input_path: str) -> list:
    """
    List all .pdb files at a local path or GCS prefix.

    Args:
        input_path: Local folder path or GCS URI prefix (gs://bucket/folder/)

    Returns:
        List of file paths or GCS URIs pointing to .pdb files.
    """
    if is_gcs_path(input_path):
        if not GCS_AVAILABLE:
            raise RuntimeError("google-cloud-storage not installed but GCS path provided.")
        bucket_name, prefix = parse_gcs_path(input_path)
        client = storage.Client()
        bucket = client.bucket(bucket_name)
        blobs = bucket.list_blobs(prefix=prefix)
        return [
            f"gs://{bucket_name}/{blob.name}"
            for blob in blobs
            if blob.name.endswith(".pdb")
        ]
    else:
        return sorted([
            os.path.join(input_path, f)
            for f in os.listdir(input_path)
            if f.endswith(".pdb")
        ])


def download_pdb(gcs_uri: str, local_dir: str) -> str:
    """
    Download a PDB file from GCS to a local temporary directory.

    Args:
        gcs_uri:   GCS URI of the PDB file (gs://bucket/path/file.pdb)
        local_dir: Local directory to download into

    Returns:
        Local path to the downloaded file.
    """
    if not GCS_AVAILABLE:
        raise RuntimeError("google-cloud-storage not installed.")

    bucket_name, blob_path = parse_gcs_path(gcs_uri)
    filename = os.path.basename(blob_path)
    local_path = os.path.join(local_dir, filename)

    client = storage.Client()
    bucket = client.bucket(bucket_name)
    blob = bucket.blob(blob_path)
    blob.download_to_filename(local_path)

    logging.info(f"Downloaded {gcs_uri} → {local_path}")
    return local_path


def upload_results(local_csv: str, output_path: str):
    """
    Upload results CSV to GCS or keep locally depending on output_path.

    Args:
        local_csv:   Local path to the results CSV file.
        output_path: Destination — either a local path or GCS URI.
    """
    if is_gcs_path(output_path):
        if not GCS_AVAILABLE:
            raise RuntimeError("google-cloud-storage not installed.")
        bucket_name, blob_path = parse_gcs_path(output_path)
        client = storage.Client()
        bucket = client.bucket(bucket_name)
        blob = bucket.blob(blob_path)
        blob.upload_from_filename(local_csv)
        logging.info(f"Uploaded results → {output_path}")
    else:
        # Already at local path — nothing to do
        logging.info(f"Results saved locally at {local_csv}")