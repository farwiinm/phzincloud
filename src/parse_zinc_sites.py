# src/parse_zinc_sites.py
"""
parse_zinc_sites.py
-------------------
Version: 1.1.0
Status:  Production-ready (validated against 1CA2, 4TLN, 3CPA, 1CDO, 4MT2, 1ZNF)

Parses PDB structure files to identify zinc ions and their
coordinating residues within a defined distance cutoff.

Validation status:
    1CA2 (Carbonic Anhydrase II)    — PASS (3His1Glu site)
    4TLN (Thermolysin)              — PASS (2His2Glu site)
    3CPA (Carboxypeptidase A)       — PASS (2His1Glu site)
    1CDO (Alcohol Dehydrogenase)    — PASS (multi-site)
    4MT2 (Metallothionein)          — PASS (multi-site Cys cluster)
    1ZNF (Zinc Finger, NMR)         — PASS (2His2Cys, Model 0 only)

Fix v1.1.0:
    - Restrict parsing to Model 0 only, preventing NMR ensemble duplication
      (e.g. 1ZNF has 37 models; v1.0 returned 148 residues instead of 4)

Known limitations:
    - Distance-based cutoff (5.0 Å default) may exclude long-bond
      coordination (>5.0 Å). Validated on primary coordination only.
    - Residue sequence numbers in PDB deposits may differ from those
      in original publications by ±1.
    - Bridging ligands in cluster-type sites (e.g. Metallothionein)
      may appear in multiple zinc_site_id records. This is expected.

Author: Fathima Farwin Mohamed Milhan
Project: pH-ZinCloud, MSc Big Data Analytics, RGU
"""

import os
import logging
from typing import Optional
from Bio.PDB import PDBParser
from Bio.PDB.NeighborSearch import NeighborSearch

# ── Logging setup ─────────────────────────────────────────────────────────────
os.makedirs("logs", exist_ok=True)
logging.basicConfig(
    filename="logs/parser.log",
    level=logging.INFO,
    format="%(asctime)s — %(levelname)s — %(message)s"
)

# ── Constants ──────────────────────────────────────────────────────────────────
ZINC_LIGAND_RESIDUES = {"CYS", "HIS", "ASP", "GLU"}
ZINC_RESNAME         = "ZN"
DEFAULT_CUTOFF       = 5.0


def get_all_atoms(model) -> list:
    """
    Extract all atoms from a single Biopython Model object.

    Args:
        model: A Biopython Model object (structure[0]).

    Returns:
        A flat list of all Atom objects in the model.
    """
    return list(model.get_atoms())


def find_zinc_atoms(model) -> list:
    """
    Find all zinc atoms in a single model.

    Args:
        model: A Biopython Model object (structure[0]).

    Returns:
        A list of Atom objects representing zinc ions.
    """
    zinc_atoms = []

    for chain in model:
        for residue in chain:
            if residue.resname == ZINC_RESNAME:
                atoms = list(residue.get_atoms())
                if atoms:
                    zinc_atoms.append(atoms[0])

    return zinc_atoms


def find_coordinating_residues(
    zinc_atom,
    all_atoms: list,
    cutoff: float = DEFAULT_CUTOFF
) -> list:
    """
    Finding all coordinating residues within cutoff Angstroms of a zinc atom.

    Using Biopython's NeighborSearch for efficient spatial lookup.
    Filters to only ZINC_LIGAND_RESIDUES (CYS, HIS, ASP, GLU) that are
    standard amino acids (not HETATM). Applies a secondary coordinating
    atom distance check to prevent distal atoms from triggering inclusion.

    Args:
        zinc_atom:  The zinc Atom object to search around.
        all_atoms:  All atoms in the model (used to build spatial index).
        cutoff:     Search radius in Angstroms.

    Returns:
        A list of Residue objects that coordinate the zinc ion.
    """
    COORD_ATOMS = {
        "HIS": ["NE2", "ND1"],
        "CYS": ["SG"],
        "ASP": ["OD1", "OD2"],
        "GLU": ["OE1", "OE2"],
    }

    ns = NeighborSearch(all_atoms)
    nearby_atoms = ns.search(zinc_atom.coord, cutoff, level="A")

    coordinating_residues = []
    seen_residue_ids = set()

    for atom in nearby_atoms:
        parent_residue = atom.get_parent()

        if parent_residue.resname == ZINC_RESNAME:
            continue
        if parent_residue.resname not in ZINC_LIGAND_RESIDUES:
            continue
        if parent_residue.id[0] != " ":
            continue

        residue_uid = (
            parent_residue.get_parent().id,
            parent_residue.id[1]
        )
        if residue_uid in seen_residue_ids:
            continue

        # Secondary coordinating atom distance check
        preferred = COORD_ATOMS.get(parent_residue.resname, [])
        coord_atom_distance = float("inf")

        for atom_name in preferred:
            if parent_residue.has_id(atom_name):
                dist = zinc_atom - parent_residue[atom_name]
                if dist < coord_atom_distance:
                    coord_atom_distance = dist

        if coord_atom_distance == float("inf"):
            if parent_residue.has_id("CA"):
                coord_atom_distance = zinc_atom - parent_residue["CA"]

        if coord_atom_distance > cutoff:
            continue

        seen_residue_ids.add(residue_uid)
        coordinating_residues.append(parent_residue)

    return coordinating_residues


def get_coordinating_atom(residue, zinc_atom) -> tuple:
    """
    Find the specific atom within a residue closest to zinc.

    Returns:
        A tuple of (atom_name, distance_in_angstroms).
    """
    preferred_atoms = {
        "HIS": ["NE2", "ND1"],
        "CYS": ["SG"],
        "ASP": ["OD1", "OD2"],
        "GLU": ["OE1", "OE2"],
    }

    candidates    = preferred_atoms.get(residue.resname, [])
    best_atom_name = "CA"
    best_distance  = float("inf")

    for atom_name in candidates:
        if residue.has_id(atom_name):
            distance = zinc_atom - residue[atom_name]
            if distance < best_distance:
                best_distance  = distance
                best_atom_name = atom_name

    if best_distance == float("inf"):
        if residue.has_id("CA"):
            best_distance  = zinc_atom - residue["CA"]
            best_atom_name = "CA"

    return best_atom_name, round(best_distance, 3)


def get_cb_coords(residue) -> Optional[tuple]:
    """
    Get CB atom coordinates for 3D visualisation.
    Falls back to CA for Glycine.
    """
    for atom_name in ["CB", "CA"]:
        if residue.has_id(atom_name):
            coord = residue[atom_name].coord
            return (round(float(coord[0]), 3),
                    round(float(coord[1]), 3),
                    round(float(coord[2]), 3))
    return None


def classify_site_type(residue_names: list) -> str:
    """
    Classify a zinc site by coordinating residue composition.
    e.g. ['HIS','HIS','HIS','GLU'] → '3His1Glu'
    """
    counts = {}
    short  = {"HIS": "His", "CYS": "Cys", "ASP": "Asp", "GLU": "Glu"}

    for name in residue_names:
        label = short.get(name, name)
        counts[label] = counts.get(label, 0) + 1

    order = ["His", "Cys", "Asp", "Glu"]
    parts = [f"{counts[k]}{k}" for k in order if k in counts]
    return "".join(parts) if parts else "Unknown"


def parse_zinc_sites(
    pdb_file: str,
    cutoff: float = DEFAULT_CUTOFF
) -> list:
    """
    Main function: parse a PDB file and return all zinc coordination sites.

    Always uses Model 0 only — handles both X-ray (single model) and
    NMR (multi-model ensemble) structures correctly.

    Args:
        pdb_file:  Path to a .pdb file.
        cutoff:    Distance cutoff in Angstroms. Default 5.0.

    Returns:
        List of dicts, one per coordinating residue. Empty list if no zinc.
    """
    pdb_id = os.path.splitext(os.path.basename(pdb_file))[0].upper()
    logging.info(f"Parsing {pdb_id} from {pdb_file}")

    try:
        parser    = PDBParser(QUIET=True)
        structure = parser.get_structure(pdb_id, pdb_file)
    except Exception as e:
        logging.error(f"Failed to parse {pdb_file}: {e}")
        return []

    # ── FIX: Always use Model 0 only ──────────────────────────────────────────
    # NMR structures contain multiple models (conformational ensemble).
    # Using all models would multiply every residue by the model count.
    # Model 0 is the representative/best model by PDB convention.
    try:
        model = structure[0]
    except KeyError:
        logging.error(f"{pdb_id}: No model 0 found in structure.")
        return []

    all_atoms  = get_all_atoms(model)

    if not all_atoms:
        logging.warning(f"{pdb_id}: No atoms found in model 0.")
        return []

    zinc_atoms = find_zinc_atoms(model)

    if not zinc_atoms:
        logging.info(f"{pdb_id}: No zinc atoms found in model 0.")
        return []

    logging.info(f"{pdb_id}: Found {len(zinc_atoms)} zinc atom(s) in model 0.")

    results = []

    for zn_atom in zinc_atoms:
        zn_chain     = zn_atom.get_parent().get_parent().id
        zn_seq       = zn_atom.get_parent().id[1]
        zinc_site_id = f"{pdb_id}_ZN_{zn_seq}_{zn_chain}"

        coord_residues = find_coordinating_residues(zn_atom, all_atoms, cutoff)

        if not coord_residues:
            logging.warning(
                f"{pdb_id}: Zinc at {zinc_site_id} has no coordinating "
                f"residues within {cutoff} Å. Skipping."
            )
            continue

        residue_names = [r.resname for r in coord_residues]
        site_type     = classify_site_type(residue_names)
        total_ligands = len(coord_residues)

        logging.info(f"{pdb_id}: Site {zinc_site_id} — {site_type} ({total_ligands} ligands)")

        for residue in coord_residues:
            coord_atom, distance = get_coordinating_atom(residue, zn_atom)
            cb_coords            = get_cb_coords(residue)

            results.append({
                "pdb_id":           pdb_id,
                "zinc_site_id":     zinc_site_id,
                "zinc_chain":       zn_chain,
                "zinc_seq_num":     zn_seq,
                "site_type":        site_type,
                "residue_name":     residue.resname,
                "residue_chain":    residue.get_parent().id,
                "residue_seq_num":  residue.id[1],
                "coord_atom":       coord_atom,
                "distance_to_zinc": distance,
                "cb_coords":        cb_coords,
                "total_ligands":    total_ligands,
            })

    return results


def summarise_site(site_data: list) -> dict:
    """
    Summarise residue-level records for one zinc site into a readable dict.
    """
    if not site_data:
        return {}

    residue_parts = [
        f"{r['residue_name']}{r['residue_seq_num']}"
        f"({r['coord_atom']}, {r['distance_to_zinc']}Å)"
        for r in site_data
    ]
    avg_dist = round(
        sum(r["distance_to_zinc"] for r in site_data) / len(site_data), 3
    )

    return {
        "zinc_site_id":    site_data[0]["zinc_site_id"],
        "site_type":       site_data[0]["site_type"],
        "total_ligands":   site_data[0]["total_ligands"],
        "residue_summary": ", ".join(residue_parts),
        "avg_distance":    avg_dist,
    }