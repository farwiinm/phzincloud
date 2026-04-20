# pH-ZinCloud

> **A cloud-native pipeline for pH-dependent zinc-binding site stability analysis across the human proteome**

[![Python](https://img.shields.io/badge/Python-3.11-blue.svg)](https://python.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![GCP](https://img.shields.io/badge/Cloud-Google%20Cloud%20Platform-orange.svg)](https://cloud.google.com)
[![Streamlit](https://img.shields.io/badge/Dashboard-Streamlit-red.svg)](https://streamlit.io)

**MSc Big Data Analytics — Robert Gordon University**  
**Author:** Fathima Farwin Mohamed Milhan  
**Academic Year:** 2025–2026

---

## What is pH-ZinCloud?

Zinc stabilises roughly 10% of all human proteins. Whether those proteins stay intact depends not only on their structure but on the **pH of the environment** they are in — a fact that every existing zinc-site prediction tool ignores.

pH-ZinCloud fills this gap. It applies the **Henderson–Hasselbalch equation** as a probabilistic scoring framework to estimate, for every coordinating residue in a zinc-binding site, the probability it remains available for coordination at any given pH. The joint product of those probabilities is the **site stability score** — a number between 0 (completely disrupted) and 1 (fully stable).

Run on 989 human zinc-binding proteins from the RCSB PDB:

| Metric | Value |
|---|---|
| Proteins analysed | 989 |
| Unique zinc sites identified | 2,661 |
| Coordinating residue records | 9,998 |
| pH-switch candidates (33.2%) | 883 |
| Tier 2 pKa coverage (PROPKA) | 66.5% |
| Total GCP cost | $0.67 |

---

## Live Dashboard

**Try it now — no installation required:**

🔗 https://phzincloud.streamlit.app/

1. Enter any PDB ID (try `1CA2` or `3CPA`)
2. Click **Analyse**
3. Move the **pH slider** from 7.4 → 5.0 and watch the 3D viewer change from green to red

---

## Architecture

```
RCSB PDB (internet)
        ↓
  fetch_pdbs.py          ← downloads PDB files via REST API
        ↓
  gs://phzincloud-data/batch_proteins/   ← Cloud Storage
        ↓
  pipeline.py            ← Cloud Run Job (4 vCPU / 8 GiB / 3600s)
    ├── parse_zinc_sites.py    (Biopython NeighborSearch, 5Å cutoff)
    ├── pka_lookup.py          (PKAD-R → PROPKA → canonical fallback)
    └── scoring_engine.py      (Henderson–Hasselbalch joint probability)
        ↓
  gs://phzincloud-data/outputs/results_batch_cloud.csv
        ↓
  BigQuery: phzincloud.phzincloud_results.residue_scores_cloud
        ↓
  Streamlit Dashboard    ← live at public URL, no login required
```

**Four GCP services used:**
- **Artifact Registry** — stores the Docker container image
- **Cloud Storage** — holds 989 PDB files and output CSVs
- **Cloud Run Jobs** — executes the pipeline on demand, no persistent server
- **BigQuery** — stores and serves 9,998 result rows via SQL

---

## Repository Structure

```
phzincloud/
├── src/
│   ├── parse_zinc_sites.py    # PDB parser — finds zinc + coordinating residues
│   ├── pka_lookup.py          # Three-tier pKa assignment (PKAD-R/PROPKA/canonical)
│   ├── scoring_engine.py      # Henderson–Hasselbalch scoring engine
│   ├── pipeline.py            # Main orchestration — chains all modules
│   ├── gcs_utils.py           # Transparent local/GCS file routing
│   └── healthcheck.py         # Container startup verification (3 math assertions)
│
├── tests/
│   ├── test_parser.py         # 14 unit tests for the parser module
│   ├── test_edge_cases.py     # Robustness tests (missing files, empty PDB, etc.)
│   └── validate_parser.py     # Literature cross-validation (1CA2, 4TLN, 3CPA, 1CDO)
│
├── dashboard/
│   ├── app.py                 # Streamlit web dashboard
│   ├── bq_client.py           # BigQuery client + live pipeline fallback
│   └── requirements.txt       # Dashboard-specific dependencies
│
├── data/
│   └── reference/             # PKAD-R CSV (Tier 1 pKa source)
│
├── results/
│   └── results_cloudrun_pkad.csv   # Final batch results (9,998 rows, 22 columns)
│
├── Dockerfile                 # Two-stage build (builder + runtime)
├── .dockerignore
├── requirements.txt           # Pinned pipeline dependencies
├── .gitignore
├── PARSER_COMPLETE.md         # Parser validation findings
├── DEPLOYMENT_NOTES.md        # Cloud deployment decisions and workarounds
└── README.md                  # This file
```

---

## Quick Start — Use the Dashboard

The fastest way to explore results is the live Streamlit dashboard at the link above. 

**Suggested proteins to try:**

| PDB ID | Protein | Site Type | What to observe |
|---|---|---|---|
| `1CA2` | Carbonic anhydrase II | His₃ | High stability at pH 7.4 → disrupted at pH 5.0 |
| `3CPA` | Carboxypeptidase A | His₂Glu | Agree with ZincSight at pH 7.4 |
| `4TLN` | Thermolysin | His₂Glu | Agree with ZincSight at pH 7.4 |
| `1ZNF` | Zinc finger (classical) | Cys₂His₂ | Shows PROPKA Cys pKa limitation |

---

## Local Installation

### Prerequisites

- Python 3.11
- Git
- (Optional) Docker Desktop — for running the containerised pipeline

### Setup

```bash
# 1. Clone the repository
git clone https://github.com/farwiinm/phzincloud.git
cd phzincloud

# 2. Create and activate a virtual environment
python -m venv phzincloud
# Windows:
phzincloud\Scripts\activate
# macOS/Linux:
source phzincloud/bin/activate

# 3. Install dependencies
pip install -r requirements.txt

# 4. Run the test suite to verify the installation
python tests/test_parser.py
python tests/test_edge_cases.py
```

Expected output from tests:
```
T1 PASS: 1CA2 contains at least one zinc site
T2 PASS: 1CA2 has at least 3 coordinating His residues
...
T14 PASS: Malformed PDB file returns empty list without crash
All 14 tests passed.
```

### Run the pipeline locally on a single protein

```bash
# Download a test PDB file
python -c "
import requests
r = requests.get('https://files.rcsb.org/download/1CA2.pdb')
open('data/raw/1CA2.pdb', 'w').write(r.text)
"

# Run the pipeline
python src/pipeline.py
# Results saved to: results/results_batch.csv
```

### Run the dashboard locally

```bash
cd dashboard
pip install -r requirements.txt
streamlit run app.py
```
Then open `http://localhost:8501` in your browser.

---

## Running the Full Batch Pipeline

### Option A — Run locally (no GCP required)

```bash
# Download 989 proteins (takes ~15 minutes, ~3GB)
python src/fetch_pdbs.py

# Set environment variables
# Windows PowerShell:
$env:PDB_INPUT_DIR = "data/raw/batch_proteins/"
$env:OUTPUT_CSV = "results/results_batch.csv"

# Run the pipeline
python src/pipeline.py
```

### Option B — Run via Docker

```bash
# Build the container
docker build -t phzincloud:latest .

# Run with local paths
docker run \
  -e PDB_INPUT_DIR=/data/batch_proteins \
  -e OUTPUT_CSV=/results/output.csv \
  -v $(pwd)/data:/data \
  -v $(pwd)/results:/results \
  phzincloud:latest
```

### Option C — Run on Google Cloud (replicating the original analysis)

Requires GCP project with billing enabled, `gcloud` CLI installed, and Artifact Registry set up.

```bash
# Authenticate
gcloud auth login
gcloud config set project phzincloud

# Build and push container
gcloud auth print-access-token | docker login -u oauth2accesstoken \
  --password-stdin us-central1-docker.pkg.dev

docker build -t phzincloud:latest .
docker tag phzincloud:latest \
  us-central1-docker.pkg.dev/phzincloud/phzincloud-repo/pipeline:latest
docker push \
  us-central1-docker.pkg.dev/phzincloud/phzincloud-repo/pipeline:latest

# Execute the Cloud Run Job
gcloud run jobs execute phzincloud-scorer \
  --region=us-central1 \
  --wait
```

Expected output:
```
OK Running execution... 1 / 1 complete
Execution [phzincloud-scorer-XXXXX] has successfully completed.
```

---

## Results Schema

The output CSV / BigQuery table (`residue_scores_cloud`) has 22 columns:

| Column | Description |
|---|---|
| `pdb_id` | RCSB PDB identifier (e.g. `1CA2`) |
| `zinc_site_id` | Unique site identifier (e.g. `1CA2_ZN_263_A`) |
| `zinc_chain` | Chain containing the zinc atom |
| `zinc_seq_num` | Residue sequence number of zinc |
| `site_type` | Coordination motif (e.g. `3His`, `2His1Glu`) |
| `residue_name` | Residue type (`HIS`, `CYS`, `ASP`, `GLU`) |
| `residue_chain` | Chain of coordinating residue |
| `residue_seq` | Sequence number of coordinating residue |
| `coord_atom` | Specific coordinating atom (`NE2`, `SG`, `OE1`, etc.) |
| `distance_to_zinc` | Distance in Ångströms |
| `pka_value` | Assigned pKa value |
| `pka_tier` | Data source tier (1=experimental, 2=PROPKA, 3=canonical) |
| `pka_source` | Description of pKa source |
| `total_ligands` | Number of coordinating residues at this site |
| `pH_4_0_score` | Site stability score at pH 4.0 |
| `pH_5_0_score` | Site stability score at pH 5.0 |
| `pH_6_0_score` | Site stability score at pH 6.0 |
| `pH_7_0_score` | Site stability score at pH 7.0 |
| `pH_7_4_score` | Site stability score at pH 7.4 (physiological) |
| `pH_8_0_score` | Site stability score at pH 8.0 |
| `pH_9_0_score` | Site stability score at pH 9.0 |
| `is_ph_switch` | Boolean — pH-switch candidate flag |

---

## Scoring Method

**Henderson–Hasselbalch probability per residue:**
```
P(deprotonated) = 1 / (1 + 10^(pKa − pH))
```

**Site stability score (joint probability):**
```
P(site stable) = ∏ P(residue_i deprotonated)
```

**pH-switch detection criteria (all four must hold):**
- Score at pH 8.0 > 0.5 (stable at mildly alkaline conditions)
- Score at pH 6.0 < 0.5 (disrupted at endosomal/tumour microenvironment pH)
- Score drop (pH 8.0 → pH 6.0) > 0.4 (steep transition)
- Score at pH 7.0 > 0.3 (retains partial function near physiological pH)

**Three-tier pKa assignment:**

| Tier | Source | Coverage (this dataset) | Notes |
|---|---|---|---|
| 1 | PKAD-R experimental | 0% | Gold standard; no overlap with 989-protein dataset |
| 2 | PROPKA 3.5.0 | 66.5% | Protein-context-aware; known Cys overestimation |
| 3 | IUPAC canonical | 33.5% | His=6.0, Cys=8.3, Asp=3.9, Glu=4.1 |

---

## Known Limitations

**1. PROPKA Cys pKa overestimation** — PROPKA models the metal-free structure and cannot account for the pKa lowering of Cys upon zinc coordination (typical shift: 8.3 → 5.0–6.5). This causes Cys-rich sites (Cys₄, Cys₂His₂, Cys₃His) to appear falsely disrupted at physiological pH. Future enhancement: apply a metal-coordination correction factor (−2 to −3 pH units) for Cys residues identified as zinc ligands.

**2. Independence assumption** — The joint probability product treats each residue's protonation as independent. In reality, protonation states are electrostatically coupled. This is mitigated by using protein-context pKa values from PROPKA but not fully eliminated.

**3. Dataset scope** — Analysis restricted to human X-ray crystallographic structures (≤2.5 Å). AlphaFold2 structures excluded (require AlphaFill preprocessing). NMR and cryo-EM structures excluded.

**4. Threshold sensitivity** — The 33.2% pH-switch rate is model-dependent. The four threshold values (0.5, 0.5, 0.4, 0.3) are design choices, not empirically calibrated values.

---

## Validation Results

| Protein | Literature Site | pH-ZinCloud Site | Agreement |
|---|---|---|---|
| 1CA2 (carbonic anhydrase) | His94, His96, His119 | His94, His96, His119 | ✓ |
| 4TLN (thermolysin) | His142, His146, Glu166 | His142, His146, Glu166 | ✓ |
| 3CPA (carboxypeptidase A) | His69, Glu72, His196 | His69, Glu72, His196 | ✓ |
| 1CDO (alcohol dehydrogenase) | Cys46, His68, Cys175 + Cys98, Cys101, Cys104, Cys112 | All 7 residues found | ✓ |

**ZincSight benchmarking at pH 7.4:** Agreement for all His-dominated sites (1CA2, 3CPA, 4TLN, 2HF8). Divergence for Cys-rich sites (1ZNF, 1CDO, 3A43, 4MT2) — explained by PROPKA Cys pKa overestimation, not a coordination geometry disagreement.

---

## Reproducibility

```bash
# Verify the container produces identical output to local Python
python tests/compare_outputs.py \
  results/results_local.csv \
  results/results_cloudrun_pkad.csv
```
Expected: `PASS: Cloud Run output matches local output exactly. Rows: 66 | Columns: 22`

---

## Dependencies

**Core pipeline** (`requirements.txt`):

| Library | Version | Purpose |
|---|---|---|
| biopython | 1.83 | PDB parsing, NeighborSearch |
| propka | 3.5.0 | Tier 2 pKa estimation |
| pandas | 2.2.1 | Data manipulation, CSV I/O |
| numpy | 1.26.4 | Numerical computation |
| scipy | 1.13.0 | Statistical functions |
| google-cloud-storage | 2.16.0 | GCS bucket I/O |
| google-cloud-bigquery | 3.20.0 | BigQuery data loading |
| requests | 2.31.0 | PDB file HTTP downloads |

**Dashboard** (`dashboard/requirements.txt`):

| Library | Version | Purpose |
|---|---|---|
| streamlit | 1.32.0 | Web dashboard framework |
| plotly | 5.20.0 | Interactive titration curve charts |
| py3Dmol | 2.0.4 | 3D molecular viewer (WebGL) |
| google-cloud-bigquery | 3.20.0 | Live BigQuery queries |

---

## Cloud Infrastructure Reference

| Service | Resource | Configuration |
|---|---|---|
| Artifact Registry | `phzincloud-repo` | Region: us-central1 |
| Cloud Storage | `gs://phzincloud-data/` | Folders: `batch_proteins/`, `outputs/` |
| Cloud Run Job | `phzincloud-scorer` | 4 vCPU, 8 GiB RAM, 3600s timeout |
| BigQuery | `phzincloud.phzincloud_results` | Table: `residue_scores_cloud` (9,998 rows) |

**Successful execution:** `phzincloud-scorer-95xhp` — 2026-03-26 05:01:50 UTC — 1/1 complete

**Total GCP cost:** $0.67 (December 2025 – March 2026)

---

## Academic Context

This project was submitted in partial fulfilment of the requirements for the degree of MSc Big Data Analytics at Robert Gordon University (2025–2026).

**Research Questions:**
1. How effectively can a Henderson–Hasselbalch-based model identify pH-sensitive zinc sites compared to static geometry approaches?
2. What cloud-native architecture best supports high-throughput zinc-site analysis without HPC infrastructure?
3. Can large-scale analysis identify pH-switch proteins across the human zinc metalloproteome?

**Key finding:** 33.2% of analysed zinc sites (883/2,661) meet pH-switch criteria, suggesting pH-dependent zinc regulation is substantially more widespread than previously characterised by case studies.

---

## Citation

If you use pH-ZinCloud in your research, please cite:

```
Mohamed Milhan, F.F. (2026) pH-ZinCloud: A cloud-native pipeline for pH-dependent 
zinc-binding site stability analysis. MSc dissertation, Robert Gordon University. 
Available at: https://github.com/farwiinm/phzincloud
```

---

## Licence

MIT — see [LICENSE](LICENSE) for details.

Data sources: RCSB PDB (CC0), PKAD-R (CC-BY, Cai et al. 2025), PROPKA (open source).

---

## References

- Cai et al. (2025) PKAD-R. *J. Comput. Biophys. Chem.* doi:10.1142/S2737416525500164
- Hekkelman et al. (2025) ZincSight. *Protein Science* 34(11). doi:10.1002/pro.70350
- Jumper et al. (2021) AlphaFold. *Nature* 596, 583–589. doi:10.1038/s41586-021-03819-2
- Li et al. (2005) PROPKA. *Proteins* 61(4), 704–721. doi:10.1002/prot.20660
- Mechti et al. (2025) ZincSight. *Nucleic Acids Res.* 53(1). doi:10.1093/nar/gkae1099
