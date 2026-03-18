# src/explore_propka_output.py
"""
Understand PROPKA's .pka output format before writing the parser for it.
"""

import subprocess
import os

PDB_FILE = "data/raw/test_proteins/1CA2.pdb"
PKA_FILE = PDB_FILE.replace(".pdb", ".pka")

# Run PROPKA
print(f"Running PROPKA on {PDB_FILE}...")
result = subprocess.run(["propka3", PDB_FILE], capture_output=True, text=True)

if result.returncode != 0:
    print("PROPKA failed:")
    print(result.stderr)
else:
    print("PROPKA succeeded.\n")

# Read and print the .pka file
if os.path.exists(PKA_FILE):
    with open(PKA_FILE, "r") as f:
        content = f.readlines()

    print("First 80 lines of PROPKA output:")
    print("-" * 60)
    for i, line in enumerate(content[:80]):
        print(f"{i+1:3d}: {line}", end="")

    # Find the SUMMARY section (most useful for us)
    print("\n\n" + "=" * 60)
    print("SUMMARY section:")
    in_summary = False
    for line in content:
        if "SUMMARY OF THIS PREDICTION" in line:
            in_summary = True
        if in_summary:
            print(line, end="")
else:
    print(f"Expected output file not found: {PKA_FILE}")