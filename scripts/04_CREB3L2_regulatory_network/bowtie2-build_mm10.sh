#!/usr/bin/env bash
# =============================================================================
# Title: bowtie2-build_mm10.sh
# Purpose: CUT&Tag alignment
# Manuscript Figure/Table: Figure 6 and supplementary
# Release role: MAIN RELEASE WORKFLOW
# Required inputs: mm10 FASTA
# Generated outputs: mm10 Bowtie2 index
# Upstream dependencies: mm10 FASTA
# Downstream consumers: CUT&Tag alignment
# Configuration keys: CUTTAG_FASTQ_DIR, CUTTAG_RESULTS_DIR, REFERENCE_DIR, DATA_PROCESSED_DIR
# Expected environment: POSIX-compatible shell with tools listed in environment/command_line_tools.tsv
# Example run command: bash scripts/04_CREB3L2_regulatory_network/bowtie2-build_mm10.sh
# Reproducibility notes: Analysis parameters and biological logic are unchanged
# from clean_release_v2; only path/configuration wiring and comments were added.
# =============================================================================
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
# shellcheck source=../_shared/paths.sh
source "${SCRIPT_DIR}/../_shared/paths.sh"
load_project_paths
require_input "${REFERENCE_DIR}/mm10/mm10.fa" "mm10 FASTA"
"${BOWTIE2_BUILD_BIN}" --threads 8 "${REFERENCE_DIR}/mm10/mm10.fa" "${REFERENCE_DIR}/mm10/mm10"
