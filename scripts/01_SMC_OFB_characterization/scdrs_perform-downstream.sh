#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# Title: scdrs_perform-downstream.sh
# Purpose: `plot.py`; supplementary tables
# Manuscript Figure/Table: Figure 2; supplementary where indicated
# Release role: MAIN RELEASE WORKFLOW
# Required inputs: AnnData and full scDRS scores
# Generated outputs: Cell-type group statistics and gene-level results
# Upstream dependencies: AnnData and full scDRS scores
# Downstream consumers: `plot.py`; supplementary tables
# Configuration keys: DATA_PROCESSED_DIR, DATA_EXTERNAL_DIR, SMC_RESULTS_DIR, GWAS_DIR, LD_REFERENCE_DIR
# Expected environment: POSIX-compatible shell with tools listed in environment/command_line_tools.tsv
# Example run command: bash scripts/01_SMC_OFB_characterization/scdrs_perform-downstream.sh
# Reproducibility notes: Analysis parameters and biological logic are unchanged
# from clean_release_v2; only path/configuration wiring and comments were added.
# =============================================================================
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
# shellcheck source=../_shared/paths.sh
source "${SCRIPT_DIR}/../_shared/paths.sh"
load_project_paths
echo "Start: `date`"
# Tool command is selected through config/paths.env.
h5ad_file="${SMC_RESULTS_DIR}/scdrs/7data.h5ad"
out_folder="${SMC_RESULTS_DIR}/scdrs"
score_file="${SMC_RESULTS_DIR}/scdrs/zscore.full_score.gz"
group_output="${out_folder}/zscore.scdrs_group.Celltype"
gene_output="${out_folder}/zscore.scdrs_gene"
require_input "${h5ad_file}" "atlas AnnData"
require_input "${score_file}" "scDRS full score"
mkdir -p "${out_folder}"
rm -f -- "${group_output}" "${gene_output}"

"${SCDRS_BIN}" perform-downstream \
    --h5ad-file "$h5ad_file"\
    --score-file "$score_file"\
    --out-folder "$out_folder"\
    --group-analysis Celltype\
    --gene-analysis\
    --flag-filter-data True\
    --flag-raw-count False

if [ ! -s "${group_output}" ]; then
    echo "scDRS Celltype group output was not generated or is empty: ${group_output}" >&2
    exit 1
fi
if [ ! -s "${gene_output}" ]; then
    echo "scDRS gene-level output was not generated or is empty: ${gene_output}" >&2
    exit 1
fi

echo "Done: `date`"
