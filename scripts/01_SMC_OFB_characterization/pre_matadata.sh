#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# Title: pre_matadata.sh
# Purpose: Prepare the original GCST90476130 summary statistics for MAGMA.
# Manuscript Figure/Table: Figure 2 and supplementary
# Release role: MAIN RELEASE WORKFLOW
# Required inputs: GCST90476130.h.tsv.gz
# Generated outputs: magma_input.txt
# Upstream dependencies: GCST90476130 GWAS summary statistics
# Downstream consumers: GCST90476130_CKD4.sh
# Configuration keys: GWAS_DIR, SMC_RESULTS_DIR, ZCAT_BIN, AWK_BIN
# Expected environment: Shell; versions are listed in environment/command_line_tools.tsv
# Example run command: bash scripts/01_SMC_OFB_characterization/pre_matadata.sh
# Reproducibility notes: Column selection is unchanged from clean_release_v2.
# Only configured input/output paths and corrupted comment line breaks were repaired.
# =============================================================================
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
# shellcheck source=../_shared/paths.sh
source "${SCRIPT_DIR}/../_shared/paths.sh"
load_project_paths

gwas_input="${GWAS_DIR}/GCST90476130.h.tsv.gz"
magma_output_dir="${SMC_RESULTS_DIR}/magma"
magma_input="${magma_output_dir}/magma_input.txt"
require_input "${gwas_input}" "GCST90476130 summary statistics"
mkdir -p "${magma_output_dir}"

magma_input_tmp="${magma_input}.tmp.$$"
cleanup_magma_input_tmp() {
    rm -f -- "${magma_input_tmp}"
}
trap cleanup_magma_input_tmp EXIT

"${ZCAT_BIN}" "${gwas_input}" | \
"${AWK_BIN}" -F'\t' '
BEGIN {
    OFS = "\t";
}
NR == 1 {
    for (i=1; i<=NF; i++) {
        if ($i == "rsid") rsid_col = i;
        if ($i == "chromosome") chr_col = i;
        if ($i == "base_pair_location") pos_col = i;
        if ($i == "p_value") pval_col = i;
        if ($i == "n") n_col = i;
    }
    missing = "";
    if (!rsid_col) {
        missing = missing (missing == "" ? "" : ", ") "rsid";
    }
    if (!chr_col) {
        missing = missing (missing == "" ? "" : ", ") "chromosome";
    }
    if (!pos_col) {
        missing = missing (missing == "" ? "" : ", ") "base_pair_location";
    }
    if (!pval_col) {
        missing = missing (missing == "" ? "" : ", ") "p_value";
    }
    if (!n_col) {
        missing = missing (missing == "" ? "" : ", ") "n";
    }
    if (missing != "") {
        print "Missing required GWAS columns: " missing > "/dev/stderr";
        exit 2;
    }
    print "SNP", "CHR", "POS", "P", "N";
    next;
}
NR > 1 {
    snp = $(rsid_col);
    chr = $(chr_col);
    pos = $(pos_col);
    p_val = $(pval_col);
    n = $(n_col);
    print snp, chr, pos, p_val, n
}
END {
    if (NR == 0) {
        print "Missing required GWAS header: input is empty" > "/dev/stderr";
        exit 2;
    }
}' > "${magma_input_tmp}"

mv -- "${magma_input_tmp}" "${magma_input}"
trap - EXIT
