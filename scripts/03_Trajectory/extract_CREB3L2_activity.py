# =============================================================================
# Title: extract_CREB3L2_activity.py
# Purpose: Extract CREB3L2(+) regulon AUC values from the SMC-enriched subset SCENIC AUCell loom
# Manuscript Figure/Table: Figure 4G
# Release role: MAIN RELEASE WORKFLOW
# Required inputs: `SCENIC_SMC_AUC_LOOM`
# Generated outputs: `CREB3L2_ACTIVITY_CSV`, indexed by CellID
# Upstream dependencies: `SCENIC.sh smc`
# Downstream consumers: `02_plot_monocle2_results.R`
# Configuration keys: SCENIC_SMC_AUC_LOOM, CREB3L2_ACTIVITY_CSV
# Expected environment: Python; package versions are listed in environment/python_packages.tsv
# Example run command: python scripts/03_Trajectory/extract_CREB3L2_activity.py
# Reproducibility notes: Reads existing CellID and CREB3L2(+) AUCell values;
# it does not run grn, ctx or aucell and does not alter activity values.
# =============================================================================
from pathlib import Path as _ReleasePath
import sys as _release_sys
_release_shared = _ReleasePath(__file__).resolve().parents[1] / "_shared"
_release_sys.path.insert(0, str(_release_shared))
from paths import load_project_paths
PATHS = load_project_paths()
import loompy
import pandas as pd


loom_path = PATHS["SCENIC_SMC_AUC_LOOM"]
output_path = PATHS["CREB3L2_ACTIVITY_CSV"]

with loompy.connect(loom_path, validate=False) as lf:
    cell_ids = lf.ca["CellID"]
    regulon_names = lf.ca["RegulonsAUC"].dtype.names
    print("Regulon names:", regulon_names)
    creb3l2_auc = lf.ca["RegulonsAUC"]["CREB3L2(+)"]

creb3l2_df = pd.DataFrame(
    {
        "CellID": cell_ids,
        "CREB3L2_activity": creb3l2_auc,
    }
)
creb3l2_df.set_index("CellID", inplace=True)
creb3l2_df.to_csv(output_path)
