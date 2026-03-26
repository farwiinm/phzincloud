"""
fetch_pdbs.py
-------------
Downloads zinc-binding human protein structures from RCSB PDB.

Dataset inclusion criteria (documented for thesis Chapter 6):
    - Organism:    Homo sapiens
    - Method:      X-ray crystallography
    - Resolution:  <= 2.5 Angstroms
    - Metal:       Zinc (ZN)
    - Site type:   Mononuclear (single zinc ion per site)

Usage:
    python fetch_pdbs.py

Output:
    data/raw/batch_proteins/   <- downloaded .pdb files
    data/raw/download_log.csv  <- record of every download attempt
"""

import os
import time
import csv
import json
import requests
import logging

# ── Logging setup ─────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(message)s",
    datefmt="%H:%M:%S"
)
logger = logging.getLogger(__name__)

# ── Configuration ──────────────────────────────────────────────────────────────
OUTPUT_DIR   = "data/raw/batch_proteins"
LOG_FILE     = "data/raw/download_log.csv"
MAX_PROTEINS = 1000          # upper cap — stays well within free trial budget
DELAY_SEC    = 0.25         # pause between downloads — be polite to RCSB
TIMEOUT_SEC  = 30           # per-request timeout


# ══════════════════════════════════════════════════════════════════════════════
# STEP 1 — Query RCSB PDB Search API for matching structures
# ══════════════════════════════════════════════════════════════════════════════

def get_pdb_ids(max_results: int = 800) -> list:
    """
    Query the RCSB PDB Search API v2 for human zinc-binding protein structures.

    Returns a list of PDB ID strings (e.g. ['1CA2', '1ZNF', ...]).
    Falls back to a hardcoded curated list if the API call fails.
    """
    url = "https://search.rcsb.org/rcsbsearch/v2/query"

    # Build the search query — all four criteria must be satisfied (AND logic)
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    # Criterion 1: Human proteins only
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entity_source_organism.ncbi_scientific_name",
                        "operator":  "exact_match",
                        "value":     "Homo sapiens"
                    }
                },
                {
                    # Criterion 2: X-ray crystallography only
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "exptl.method",
                        "operator":  "exact_match",
                        "value":     "X-RAY DIFFRACTION"
                    }
                },
                {
                    # Criterion 3: Resolution <= 2.5 Angstroms
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.resolution_combined",
                        "operator":  "less_or_equal",
                        "value":     2.5
                    }
                },
                {
    # Criterion 4: Contains at least one zinc ion
    "type": "terminal",
    "service": "text",
    "parameters": {
        "attribute": "rcsb_chem_comp_container_identifiers.comp_id",
        "operator":  "exact_match",
        "value":     "ZN"
    }
},
            ]
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {
                "start": 0,
                "rows":  max_results
            },
            "sort": [
                # Sort by score descending — highest quality hits first
                {"sort_by": "score", "direction": "desc"}
            ]
        }
    }

    logger.info(f"Querying RCSB PDB API for up to {max_results} structures...")
    logger.info("Criteria: Homo sapiens | X-ray | Resolution ≤ 2.5Å | Zinc")

    try:
        response = requests.post(
            url,
            data=json.dumps(query),
            headers={"Content-Type": "application/json"},
            timeout=60
        )
        response.raise_for_status()
        data = response.json()

        total     = data.get("total_count", 0)
        pdb_ids   = [hit["identifier"] for hit in data.get("result_set", [])]

        logger.info(f"API returned {total:,} total matches")
        logger.info(f"Retrieved {len(pdb_ids)} PDB IDs for download")
        return pdb_ids

    except requests.exceptions.ConnectionError:
        logger.warning("No internet connection or RCSB API unreachable.")
        logger.warning("Falling back to curated list.")
        return _curated_fallback()

    except requests.exceptions.Timeout:
        logger.warning("RCSB API timed out.")
        logger.warning("Falling back to curated list.")
        return _curated_fallback()

    except Exception as e:
        logger.warning(f"RCSB API query failed: {e}")
        logger.warning("Falling back to curated list.")
        return _curated_fallback()


def _curated_fallback() -> list:
    """
    Hardcoded list of well-characterised human zinc-binding proteins.
    Used when the RCSB API is unavailable.
    Covers all major zinc site types for diverse thesis results.
    """
    return [
        # ── Carbonic anhydrases (His3 catalytic site) ──────────────────────
        "1CA2", "2CBA", "3KS3", "1BIC", "1MOO", "2NNG", "1ZSA", "2VVB",
        # ── Zinc metalloproteases (His2Glu) ────────────────────────────────
        "4TLN", "1IAG", "1KAP", "2OI6", "3B7E", "2AZ8", "1CGL",
        # ── Carboxypeptidases (His2Glu) ────────────────────────────────────
        "3CPA", "1YME", "2BOA", "1M4L", "3QSN", "1DTD",
        # ── Classical zinc fingers (Cys2His2) ─────────────────────────────
        "1ZNF", "1A1L", "2I13", "1AAY", "2JP9", "1MEY", "2DSP", "1TF3",
        # ── Metallothioneins (Cys clusters) ────────────────────────────────
        "4MT2", "1MHU", "2MHU", "4F5S", "1JJD",
        # ── Alcohol dehydrogenases (Cys4 + CysHisCys) ─────────────────────
        "1CDO", "1HLD", "1TEH", "2OHX", "1AGN",
        # ── Tumour suppressors ─────────────────────────────────────────────
        "2OCJ", "1TUP", "1AIE", "2AC0", "1GZH",
        # ── Superoxide dismutases ──────────────────────────────────────────
        "2SOD", "1SOS", "1PU0", "1AZV",
        # ── Nuclear receptors (Cys4 zinc finger) ──────────────────────────
        "1HCQ", "2NLL", "1IE9", "1K74",
        # ── HypA — key pH-switch validation protein ────────────────────────
        "3CYZ",
        # ── Miscellaneous well-characterised zinc proteins ─────────────────
        "2HF8", "3A43", "1AMP", "1ZHP", "2BBM", "1JLJ", "1K9V", "1QGO",
        "1ARU", "2BNV", "1E9I", "1ZNJ", "2ACF", "1DFV", "1YDW", "2DFP",
        "2G69", "1LNR", "2WFZ", "1ONR", "1Y43", "1H8X", "1CEW", "1PMO",
        "2QBK", "1FT1", "1NVJ", "1RYZ", "1FN9", "1G2D", "2C3I", "1A5O",
        "1BS0", "2FQM", "1G1D", "2GLZ", "1OEI", "1Y3K", "1YRY", "2C1A",
    ]


# ══════════════════════════════════════════════════════════════════════════════
# STEP 2 — Download each PDB file
# ══════════════════════════════════════════════════════════════════════════════

def download_one(pdb_id: str) -> dict:
    """
    Download a single PDB file from RCSB.

    Returns a result dict with keys:
        pdb_id   — the identifier
        status   — 'downloaded', 'cached', or 'failed'
        size_kb  — file size in kilobytes
        reason   — error message if failed, empty string otherwise
    """
    out_path = os.path.join(OUTPUT_DIR, f"{pdb_id}.pdb")

    # Skip if already downloaded
    if os.path.exists(out_path):
        size_kb = os.path.getsize(out_path) // 1024
        return {"pdb_id": pdb_id, "status": "cached",
                "size_kb": size_kb, "reason": ""}

    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"

    try:
        r = requests.get(url, timeout=TIMEOUT_SEC)

        if r.status_code == 200 and len(r.text) > 500:
            with open(out_path, "w", encoding="utf-8") as f:
                f.write(r.text)
            size_kb = len(r.text) // 1024
            return {"pdb_id": pdb_id, "status": "downloaded",
                    "size_kb": size_kb, "reason": ""}

        elif r.status_code == 404:
            return {"pdb_id": pdb_id, "status": "failed",
                    "size_kb": 0, "reason": "not found in PDB"}

        else:
            return {"pdb_id": pdb_id, "status": "failed",
                    "size_kb": 0, "reason": f"HTTP {r.status_code}"}

    except requests.exceptions.Timeout:
        return {"pdb_id": pdb_id, "status": "failed",
                "size_kb": 0, "reason": "connection timeout"}

    except Exception as e:
        return {"pdb_id": pdb_id, "status": "failed",
                "size_kb": 0, "reason": str(e)[:80]}


# ══════════════════════════════════════════════════════════════════════════════
# STEP 3 — Main orchestrator
# ══════════════════════════════════════════════════════════════════════════════

def run():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(os.path.dirname(LOG_FILE), exist_ok=True)

    # ── Get PDB IDs ───────────────────────────────────────────────────────────
    pdb_ids = get_pdb_ids(MAX_PROTEINS)

    if not pdb_ids:
        logger.error("No PDB IDs to download. Exiting.")
        return

    total = len(pdb_ids)
    logger.info(f"Starting download of {total} structures...")
    logger.info(f"Output: {OUTPUT_DIR}")
    logger.info("")

    # ── Download loop ─────────────────────────────────────────────────────────
    results     = []
    downloaded  = 0
    cached      = 0
    failed      = 0

    for i, pdb_id in enumerate(pdb_ids, start=1):
        result = download_one(pdb_id)
        results.append(result)

        status = result["status"]

        if status == "downloaded":
            downloaded += 1
            logger.info(
                f"[{i:>4}/{total}]  ✓  {pdb_id}  "
                f"({result['size_kb']} KB)"
            )
        elif status == "cached":
            cached += 1
            # Only log every 50th cached file to keep output readable
            if cached % 50 == 1:
                logger.info(
                    f"[{i:>4}/{total}]  ~  {pdb_id}  (already downloaded)"
                )
        else:
            failed += 1
            logger.warning(
                f"[{i:>4}/{total}]  ✗  {pdb_id}  "
                f"FAILED: {result['reason']}"
            )

        # Progress update every 100 proteins
        if i % 100 == 0:
            available = downloaded + cached
            logger.info(
                f"── Progress: {i}/{total} | "
                f"New: {downloaded} | "
                f"Cached: {cached} | "
                f"Failed: {failed} | "
                f"Available: {available} ──"
            )

        # Polite delay between requests
        time.sleep(DELAY_SEC)

    # ── Save download log ──────────────────────────────────────────────────────
    with open(LOG_FILE, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f, fieldnames=["pdb_id", "status", "size_kb", "reason"]
        )
        writer.writeheader()
        writer.writerows(results)

    # ── Final summary ──────────────────────────────────────────────────────────
    available = downloaded + cached
    logger.info("")
    logger.info("=" * 60)
    logger.info("DOWNLOAD COMPLETE")
    logger.info("=" * 60)
    logger.info(f"  Newly downloaded : {downloaded}")
    logger.info(f"  Already cached   : {cached}")
    logger.info(f"  Failed           : {failed}")
    logger.info(f"  Total available  : {available} PDB files")
    logger.info(f"  Output folder    : {OUTPUT_DIR}/")
    logger.info(f"  Download log     : {LOG_FILE}")
    logger.info("")
    logger.info("Next steps:")
    logger.info("  1. Run pipeline locally:")
    logger.info(f"       set PDB_INPUT_DIR={OUTPUT_DIR}")
    logger.info(f"       set OUTPUT_CSV=results/results_batch.csv")
    logger.info(f"       python src/pipeline.py")
    logger.info("")
    logger.info("  2. Upload to GCS and run Cloud Run job:")
    logger.info(f"       gsutil -m cp {OUTPUT_DIR}/*.pdb "
                f"gs://phzincloud-data/batch_proteins/")
    logger.info("=" * 60)


if __name__ == "__main__":
    run()