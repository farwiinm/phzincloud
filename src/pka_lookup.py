# src/pka_lookup.py
"""
Three-Tier pKa Assignment Module for pH-ZinCloud

Tier 1: Experimental lookup from PKAD-R database
Tier 2: PROPKA computational prediction
Tier 3: Canonical literature defaults (fallback)

Usage:
    from pka_lookup import get_pka
    result = get_pka(pdb_id, chain, residue_num, residue_name, pdb_file)
"""

import pandas as pd
import subprocess
import os
import glob
import logging

logger = logging.getLogger(__name__)

# ─────────────────────────────────────────
# TIER 3: Canonical default pKa values
# Source: Lehninger Biochemistry / standard references
# ─────────────────────────────────────────
CANONICAL_PKA = {
    "HIS": 6.0,
    "CYS": 8.3,
    "ASP": 3.9,
    "GLU": 4.1,
    "LYS": 10.5,
    "TYR": 10.1,
    "ARG": 12.5,
}


# ─────────────────────────────────────────
# Load PKAD-R once at module import
# ─────────────────────────────────────────
def _load_pkad() -> pd.DataFrame:
    """Load and clean the PKAD-R database. Called once at module level."""
    ref_files = glob.glob("data/reference/PKAD-R*.csv")
    if not ref_files:
        logger.warning("PKAD-R file not found in data/reference/. Tier 1 disabled.")
        return pd.DataFrame()

    path = ref_files[0]
    try:
        df = pd.read_csv(path)
    except Exception:
        df = pd.read_csv(path, sep="\t")

    # Keep only 'Main' pKa entries (exclude mutants and alternatives)
    df = df[df["pKa Classification"] == "Main"].copy()

    # Remove rows flagged with warnings
    df = df[df["Warning"].isna()].copy()

    # Convert Expt. pKa to float (it is stored as string)
    df["pKa_value"] = pd.to_numeric(df["Expt. pKa"], errors="coerce")
    df = df.dropna(subset=["pKa_value"])

    # Normalise PDB IDs to uppercase for consistent matching
    df["PDB"] = df["PDB"].str.upper().str.strip()
    df["Chain"] = df["Chain"].str.strip()
    df["ResName"] = df["ResName"].str.upper().str.strip()

    logger.info(f"PKAD-R loaded: {len(df)} usable rows from {path}")
    return df


# Load at module level so it's only read from disk once
_PKAD_DF = _load_pkad()
_PROPKA_CACHE: dict = {}


# ─────────────────────────────────────────
# TIER 1: Experimental lookup from PKAD-R
# ─────────────────────────────────────────
def lookup_experimental_pka(
    pdb_id: str,
    chain: str,
    residue_num: int,
    residue_name: str
) -> dict | None:
    """
    Look up an experimental pKa from PKAD-R.

    Returns a dict with keys: pka, source, tier
    Returns None if not found.
    """
    if _PKAD_DF.empty:
        return None

    mask = (
        (_PKAD_DF["PDB"] == pdb_id.upper().strip()) &
        (_PKAD_DF["Chain"] == chain.strip()) &
        (_PKAD_DF["ResID in PDB"] == int(residue_num)) &
        (_PKAD_DF["ResName"] == residue_name.upper().strip())
    )

    matches = _PKAD_DF[mask]

    if matches.empty:
        return None

    # Take the first match (already filtered to 'Main' only)
    row = matches.iloc[0]
    return {
        "pka": float(row["pKa_value"]),
        "source": "PKAD-R (experimental)",
        "tier": 1,
        "uncertainty": row["Expt. Uncertainty"] if pd.notna(row["Expt. Uncertainty"]) else "N/A"
    }


# ─────────────────────────────────────────
# TIER 2: PROPKA computational prediction
# ─────────────────────────────────────────
def _parse_propka_output(pka_file: str, chain: str, residue_num: int) -> float | None:
    """
    Parse a PROPKA .pka output file and extract pKa for a specific residue.
    Looks in the SUMMARY section for the matching chain + residue number.
    """
    if not os.path.exists(pka_file):
        return None

    with open(pka_file, "r") as f:
        lines = f.readlines()

    in_summary = False
    for line in lines:
        if "SUMMARY OF THIS PREDICTION" in line:
            in_summary = True
            continue

        if in_summary:
            # PROPKA summary lines look like:
            # HIS  94 A      6.54     6.50
            parts = line.split()
            if len(parts) >= 4:
                try:
                    res_num = int(parts[1])
                    res_chain = parts[2]
                    pka_val = float(parts[3])
                    if res_num == int(residue_num) and res_chain == chain:
                        return pka_val
                except (ValueError, IndexError):
                    continue

    return None


def estimate_pka_propka(pdb_file: str, chain: str, residue_num: int) -> dict | None:
    """
    Run PROPKA using the Python library API and return pKa for a specific residue.
    Uses atom.res_num (not atom.res_seq) which is the correct attribute in propka 3.5.0.
    Caches the full MolecularContainer per PDB file to avoid re-running for each residue.
    """
    # Prevent unbounded memory growth during large batch runs
    # Cache only needs current protein — clear when it exceeds 50 entries
    if len(_PROPKA_CACHE) > 50:
        _PROPKA_CACHE.clear()

    if not os.path.exists(pdb_file):
        logger.warning(f"PDB file not found for PROPKA: {pdb_file}")
        return None

    try:
        import propka.run as pk

        # Run PROPKA once per PDB file and cache the result
        if pdb_file not in _PROPKA_CACHE:
            try:
                mol = pk.single(pdb_file, optargs=["--quiet"])
                # Use conformation '1A' (first model) — 'AVR' is the average
                conf_key = '1A' if '1A' in mol.conformations else list(mol.conformations.keys())[0]
                _PROPKA_CACHE[pdb_file] = mol.conformations[conf_key]
                logger.info(f"PROPKA ran on {os.path.basename(pdb_file)}: "
                            f"{len(_PROPKA_CACHE[pdb_file].groups)} groups found")
            except Exception as e:
                logger.warning(f"PROPKA run failed on {pdb_file}: {e}")
                _PROPKA_CACHE[pdb_file] = None

        conf = _PROPKA_CACHE.get(pdb_file)
        if conf is None:
            return None

        # Search groups for the matching chain + residue number
        for group in conf.groups:
            try:
                atom = group.atom
                g_chain  = atom.chain_id
                g_resnum = int(atom.res_num)       # ← correct attribute
                g_pka    = group.pka_value

                if g_chain == chain and g_resnum == int(residue_num):
                    if g_pka is not None and 0.5 <= float(g_pka) <= 14.0:
                        return {
                            "pka":    round(float(g_pka), 2),
                            "source": "PROPKA (computational)",
                            "tier":   2
                        }
            except (AttributeError, ValueError, TypeError):
                continue

    except ImportError:
        logger.warning("PROPKA Python library not installed. Tier 2 disabled.")
    except Exception as e:
        logger.warning(f"PROPKA unexpected error on {pdb_file}: {e}")

    return None


# ─────────────────────────────────────────
# TIER 3: Canonical default fallback
# ─────────────────────────────────────────
def get_canonical_pka(residue_name: str) -> dict | None:
    """Return the canonical textbook pKa for a residue type."""
    pka_val = CANONICAL_PKA.get(residue_name.upper())
    if pka_val is None:
        return None
    return {
        "pka": pka_val,
        "source": "Canonical default (Lehninger)",
        "tier": 3
    }


# ─────────────────────────────────────────
# MAIN ENTRY POINT: Tiered lookup
# ─────────────────────────────────────────
def get_pka(
    pdb_id: str,
    chain: str,
    residue_num: int,
    residue_name: str,
    pdb_file: str = None
) -> dict:
    """
    Get pKa for a residue using the three-tier system.

    Tries Tier 1 (experimental) → Tier 2 (PROPKA) → Tier 3 (canonical).
    Always returns a result — never raises.

    Returns dict with keys: pka, source, tier
    """
    # Tier 1: Experimental
    result = lookup_experimental_pka(pdb_id, chain, residue_num, residue_name)
    if result:
        logger.debug(f"Tier 1 hit: {pdb_id} {chain} {residue_num} {residue_name} → {result['pka']}")
        return result

    # Tier 2: PROPKA (only if pdb_file provided)
    if pdb_file:
        result = estimate_pka_propka(pdb_file, chain, residue_num)
        if result:
            logger.debug(f"Tier 2 hit: {pdb_id} {chain} {residue_num} {residue_name} → {result['pka']}")
            return result

    # Tier 3: Canonical default
    result = get_canonical_pka(residue_name)
    if result:
        logger.debug(f"Tier 3 fallback: {residue_name} → {result['pka']}")
        return result

    # Should never reach here for standard residues
    logger.warning(f"No pKa found for {pdb_id} {chain} {residue_num} {residue_name}")
    return {"pka": None, "source": "Not found", "tier": None}