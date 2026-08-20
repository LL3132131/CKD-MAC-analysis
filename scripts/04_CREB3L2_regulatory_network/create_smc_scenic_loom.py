# =============================================================================
# Title: create_smc_scenic_loom.py
# Purpose: Create the SMC-enriched subset expression loom for pySCENIC
# Manuscript Figure/Table: Figure 4G and supplementary
# Release role: MAIN RELEASE WORKFLOW
# Required inputs: SMC-enriched subset `sce_exp.csv`
# Generated outputs: `SCENIC_SMC_LOOM`
# Upstream dependencies: `03-3_sub-SCENIC.R`
# Downstream consumers: `SCENIC.sh smc`
# Configuration keys: SCENIC_RESULTS_DIR, SCENIC_SMC_LOOM
# Expected environment: Python; package versions are listed in environment/python_packages.tsv
# Example run command: python scripts/04_CREB3L2_regulatory_network/create_smc_scenic_loom.py
# Reproducibility notes: This creates the input loom only and does not calculate
# regulons or AUCell values.
# =============================================================================
from pathlib import Path as _ReleasePath
import sys as _release_sys
_release_shared = _ReleasePath(__file__).resolve().parents[1] / "_shared"
_release_sys.path.insert(0, str(_release_shared))
from paths import configured_path, load_project_paths
PATHS = load_project_paths()
import os

import loompy as lp
import numpy as np
import scanpy as sc


work_dir = configured_path(PATHS, "SCENIC_RESULTS_DIR", "01.sub_allCell")
os.chdir(work_dir)
x = sc.read_csv(str(work_dir / "sce_exp.csv"))
row_attrs = {"Gene": np.array(x.var_names)}
col_attrs = {"CellID": np.array(x.obs_names)}
lp.create(str(PATHS["SCENIC_SMC_LOOM"]), x.X.transpose(), row_attrs, col_attrs)
