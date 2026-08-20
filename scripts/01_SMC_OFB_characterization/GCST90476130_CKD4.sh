#!/usr/bin/env bash
set -euo pipefail
# =============================================================================
# Title: GCST90476130_CKD4.sh
# Purpose: Run the original MAGMA annotation and gene-level analysis.
# Manuscript Figure/Table: Figure 2 and supplementary
# Release role: MAIN RELEASE WORKFLOW
# Required inputs: magma_input.txt; NCBI38 gene locations; European LD reference
# Generated outputs: MAGMA annotation and gene-level result files
# Upstream dependencies: pre_matadata.sh
# Downstream consumers: mk_gs.py
# Configuration keys: SMC_RESULTS_DIR, GWAS_DIR, LD_REFERENCE_DIR, MAGMA_BIN
# Expected environment: Shell; versions are listed in environment/command_line_tools.tsv
# Example run command: bash scripts/01_SMC_OFB_characterization/GCST90476130_CKD4.sh
# Reproducibility notes: MAGMA arguments are unchanged from clean_release_v2.
# Only command/path selection and corrupted comment line breaks were repaired.
# =============================================================================
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
# shellcheck source=../_shared/paths.sh
source "${SCRIPT_DIR}/../_shared/paths.sh"
load_project_paths

out_dir="${SMC_RESULTS_DIR}/magma"
magma_bin="${MAGMA_BIN}"
snp_loc_dir="${SMC_RESULTS_DIR}/magma/magma_input.txt"
gene_loc_dir="${GWAS_DIR}/NCBI38/NCBI38.gene.loc"
gene_eas="${LD_REFERENCE_DIR}/g1000_eur"
annotation_output="${out_dir}/magma_annotation.genes.annot"
gene_output="${out_dir}/step2.genes.out"
require_input "${snp_loc_dir}" "MAGMA input"
require_input "${gene_loc_dir}" "NCBI38 gene locations"
require_input "${gene_eas}.bed" "European LD BED"
require_input "${gene_eas}.bim" "European LD BIM"
require_input "${gene_eas}.fam" "European LD FAM"
mkdir -p "${out_dir}"

# Step 1: gene annotation.
rm -f -- "${annotation_output}"
"${magma_bin}" \
    --annotate window=10,10 \
    --snp-loc "${snp_loc_dir}" \
    --gene-loc "${gene_loc_dir}" \
    --out "${out_dir}/magma_annotation"

# Step 2: verify the expected annotation output.
if [ ! -s "${annotation_output}" ]; then
    echo "MAGMA gene annotation output was not generated or is empty: ${annotation_output}" >&2
    exit 1
fi

# Step 3: gene-level analysis.
rm -f -- "${gene_output}"
"${magma_bin}" \
    --bfile "${gene_eas}" \
    --pval "${snp_loc_dir}" use='SNP,P' ncol='N' \
    --gene-annot "${annotation_output}" \
    --out "${out_dir}/step2"

if [ ! -s "${gene_output}" ]; then
    echo "MAGMA gene-level output was not generated or is empty: ${gene_output}" >&2
    exit 1
fi
