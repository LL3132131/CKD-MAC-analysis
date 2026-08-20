#!/usr/bin/env bash
# =============================================================================
# Title: CREB3L2-P_trim_first.sh
# Purpose: alignment/peak calling
# Manuscript Figure/Table: Figure 6 and supplementary
# Release role: MAIN RELEASE WORKFLOW
# Required inputs: paired CREB3L2 CUT&Tag FASTQ
# Generated outputs: trimmed paired FASTQ
# Upstream dependencies: paired CREB3L2 CUT&Tag FASTQ
# Downstream consumers: alignment/peak calling
# Configuration keys: CUTTAG_FASTQ_DIR, CUTTAG_RESULTS_DIR, REFERENCE_DIR, DATA_PROCESSED_DIR
# Expected environment: POSIX-compatible shell with tools listed in environment/command_line_tools.tsv
# Example run command: bash scripts/04_CREB3L2_regulatory_network/CREB3L2-P_trim_first.sh
# Reproducibility notes: Analysis parameters and biological logic are unchanged
# from clean_release_v2; only path/configuration wiring and comments were added.
# =============================================================================
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
# shellcheck source=../_shared/paths.sh
source "${SCRIPT_DIR}/../_shared/paths.sh"
load_project_paths
require_input "${CUTTAG_FASTQ_DIR}/CREB3L2-P_1.fq.gz" "CREB3L2-P read 1"
require_input "${CUTTAG_FASTQ_DIR}/CREB3L2-P_2.fq.gz" "CREB3L2-P read 2"
"${TRIM_GALORE_BIN}" -q 25 --phred33 --length 36 --stringency 3 --paired "${CUTTAG_FASTQ_DIR}/CREB3L2-P_1.fq.gz" "${CUTTAG_FASTQ_DIR}/CREB3L2-P_2.fq.gz" --gzip -o "${CUTTAG_RESULTS_DIR}/trimmed"
