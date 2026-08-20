#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# Title: scdrs-munge-gs.sh
# Purpose: scDRS score
# Manuscript Figure/Table: Figure 2; supplementary where indicated
# Release role: MAIN RELEASE WORKFLOW
# Required inputs: MAGMA symbol/Z-score table
# Generated outputs: Top-1000 `disease_geneset.gs`
# Upstream dependencies: MAGMA symbol/Z-score table
# Downstream consumers: scDRS score
# Configuration keys: DATA_PROCESSED_DIR, DATA_EXTERNAL_DIR, SMC_RESULTS_DIR, GWAS_DIR, LD_REFERENCE_DIR
# Expected environment: POSIX-compatible shell with tools listed in environment/command_line_tools.tsv
# Example run command: bash scripts/01_SMC_OFB_characterization/scdrs-munge-gs.sh
# Reproducibility notes: Analysis parameters and biological logic are unchanged
# from clean_release_v2; only path/configuration wiring and comments were added.
# =============================================================================
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
# shellcheck source=../_shared/paths.sh
source "${SCRIPT_DIR}/../_shared/paths.sh"
load_project_paths
echo "Start: `date`"
# Tool command is selected through config/paths.env.
out_dir="${SMC_RESULTS_DIR}/scdrs"
zscore="${SMC_RESULTS_DIR}/magma/magma_zscore_symbol.tsv"
gs_output="${out_dir}/disease_geneset.gs"
require_input "${zscore}" "MAGMA symbol/Z-score table"
mkdir -p "${out_dir}"
rm -f -- "${gs_output}"

"${SCDRS_BIN}" munge-gs \
    --out-file "${gs_output}" \
    --zscore-file "$zscore" \
    --weight zscore \
    --n-max 1000

if [ ! -s "${gs_output}" ]; then
    echo "scDRS gene-set output was not generated or is empty: ${gs_output}" >&2
    exit 1
fi

echo "Done: `date`"
