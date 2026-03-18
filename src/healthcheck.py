# src/healthcheck.py
"""
Container healthcheck for pH-ZinCloud.
Verifies all imports work and core maths are correct.
Exits 0 on success, 1 on any failure.
"""
import sys
import os

# Same path fix as pipeline.py
_SRC_DIR = os.path.dirname(os.path.abspath(__file__))
if _SRC_DIR not in sys.path:
    sys.path.insert(0, _SRC_DIR)

print("pH-ZinCloud container healthcheck starting...")

# ── Imports ───────────────────────────────────────────────────────────────────
checks = [
    ("Biopython",         "from Bio.PDB import PDBParser, NeighborSearch"),
    ("pandas",            "import pandas as pd"),
    ("numpy",             "import numpy as np"),
    ("parse_zinc_sites",  "from parse_zinc_sites import parse_zinc_sites, classify_site_type"),
    ("pka_lookup",        "from pka_lookup import get_pka"),
    ("scoring_engine",    "from scoring_engine import site_stability_score, interpret_score, henderson_hasselbalch"),
]

all_ok = True
for name, stmt in checks:
    try:
        exec(stmt)
        print(f"  [OK] {name}")
    except Exception as e:
        print(f"  [FAIL] {name}: {e}")
        all_ok = False

if not all_ok:
    print("\nHealthcheck failed at import stage.")
    sys.exit(1)

# ── Maths verification ────────────────────────────────────────────────────────
from scoring_engine import henderson_hasselbalch, site_stability_score

math_checks = [
    ("H-H: pH==pKa gives 0.5",      lambda: abs(henderson_hasselbalch(6.0, 6.0) - 0.5) < 1e-10),
    ("H-H: pH>>pKa gives ~1.0",     lambda: henderson_hasselbalch(9.0, 6.0) > 0.99),
    ("H-H: pH<<pKa gives ~0.0",     lambda: henderson_hasselbalch(3.0, 6.0) < 0.01),
    ("3-His site at pH 6.0 ~ 0.125",lambda: abs(
        site_stability_score([("HIS",6.0,3),("HIS",6.0,3),("HIS",6.0,3)], 6.0)["overall_score"] - 0.125
    ) < 0.01),
]

for name, check in math_checks:
    try:
        assert check(), f"assertion failed"
        print(f"  [OK] {name}")
    except Exception as e:
        print(f"  [FAIL] {name}: {e}")
        all_ok = False

if not all_ok:
    print("\nHealthcheck failed at maths verification stage.")
    sys.exit(1)

print("\nAll checks passed. Container is ready.")
sys.exit(0)