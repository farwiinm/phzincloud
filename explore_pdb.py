import glob
import numpy as np
from Bio.PDB import PDBParser

def atom_distance(atom1, atom2):
    diff = atom1.coord - atom2.coord
    return np.sqrt(np.sum(diff ** 2))

parser = PDBParser(QUIET=True)

# 1. Locate all PDB files in the directory
pdb_files = glob.glob("data/raw/test_proteins/*.pdb")

# 2. Loop through the first 6 files found
for file_path in pdb_files[:6]:
    # Extract filename for the ID (e.g., '1CA2')
    struct_id = file_path.split("/")[-1].replace(".pdb", "")
    
    print(f"\n{'='*40}")
    print(f"PROCESSING PROTEIN: {struct_id}")
    print(f"{'='*40}")

    structure = parser.get_structure(struct_id, file_path)

    # --- Summary of Chains ---
    for model in structure:
        for chain in model:
            print(f"Chain: {chain.id}, Residues: {len(list(chain.get_residues()))}")

    # --- Find Non-standard Residues ---
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] != " ":  
                    print(f"Non-standard: {residue.resname}, Chain {chain.id}, Seq {residue.id[1]}")

    # --- Zinc Residue Identification ---
    zinc_residues = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.resname == "ZN":
                    zinc_residues.append(residue)
                    print(f"Found zinc in chain {chain.id} at position {residue.id[1]}")

    # --- Coordinate Printing ---
    for zn in zinc_residues:
        zn_atom = list(zn.get_atoms())[0]
        print(f"Zinc coordinates for {struct_id}: {zn_atom.coord}")

print(f"\nFinished processing {len(pdb_files[:6])} proteins.")
