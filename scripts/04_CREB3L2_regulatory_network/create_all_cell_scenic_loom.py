# =============================================================================
# Title: create_all_cell_scenic_loom.py
# Purpose: `SCENIC.sh`
# Manuscript Figure/Table: Figure 4 and supplementary
# Release role: MAIN RELEASE WORKFLOW
# Required inputs: all-cell `sce_exp.csv`
# Generated outputs: `rSCENIC.loom`
# Upstream dependencies: all-cell `sce_exp.csv`
# Downstream consumers: `SCENIC.sh`
# Configuration keys: SCENIC_RESULTS_DIR, SCENIC_ALL_CELL_LOOM
# Expected environment: Python; package versions are listed in environment/python_packages.tsv
# Example run command: python scripts/04_CREB3L2_regulatory_network/create_all_cell_scenic_loom.py
# Reproducibility notes: Analysis parameters and biological logic are unchanged
# from clean_release_v2; only path/configuration wiring and comments were added.
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


work_dir = configured_path(PATHS, "SCENIC_RESULTS_DIR", "00.allCell", "7scenic")
os.chdir(work_dir)
x = sc.read_csv(str(work_dir / "sce_exp.csv"))
row_attrs = {"Gene": np.array(x.var_names)}
col_attrs = {"CellID": np.array(x.obs_names)}
lp.create(str(PATHS["SCENIC_ALL_CELL_LOOM"]), x.X.transpose(), row_attrs, col_attrs)
