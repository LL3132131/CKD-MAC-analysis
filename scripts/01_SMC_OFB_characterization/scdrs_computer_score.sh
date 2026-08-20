#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# Title: scdrs_computer_score.sh
# Purpose: Downstream group analysis
# Manuscript Figure/Table: Figure 2; supplementary where indicated
# Release role: MAIN RELEASE WORKFLOW
# Required inputs: `7data.h5ad`; disease gene set; de-identified covariates
# Generated outputs: Full scDRS score file
# Upstream dependencies: `7data.h5ad`; disease gene set; de-identified covariates
# Downstream consumers: Downstream group analysis
# Configuration keys: DATA_PROCESSED_DIR, DATA_EXTERNAL_DIR, SMC_RESULTS_DIR, GWAS_DIR, LD_REFERENCE_DIR
# Expected environment: POSIX-compatible shell with tools listed in environment/command_line_tools.tsv
# Example run command: bash scripts/01_SMC_OFB_characterization/scdrs_computer_score.sh
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
gs_file="${SMC_RESULTS_DIR}/scdrs/disease_geneset.gs"
out_folder="${SMC_RESULTS_DIR}/scdrs"
cov="${DATA_EXTERNAL_DIR}/scdrs/covariates.tsv"
score_output="${out_folder}/zscore.full_score.gz"
require_input "${h5ad_file}" "atlas AnnData"
require_input "${gs_file}" "scDRS gene set"
require_input "${cov}" "de-identified covariates"
mkdir -p "$out_folder"
rm -f -- "${score_output}"

"${SCDRS_BIN}" compute-score \
    --h5ad-file "$h5ad_file"\
    --h5ad-species human\
    --gs-file "$gs_file"\
    --gs-species human\
    --out-folder "$out_folder"\
    --cov-file "$cov"\
    --flag-filter-data True\
    --flag-raw-count False\
    --n-ctrl 1000\
    --flag-return-ctrl-raw-score False\
    --flag-return-ctrl-norm-score True

if [ ! -s "${score_output}" ]; then
    echo "scDRS full-score output was not generated or is empty: ${score_output}" >&2
    exit 1
fi

echo "Done: `date`"
