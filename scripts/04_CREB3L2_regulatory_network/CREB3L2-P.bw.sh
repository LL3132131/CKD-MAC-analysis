#!/usr/bin/env bash
# =============================================================================
# Title: CREB3L2-P.bw.sh
# Purpose: Figure 6 coverage panel
# Manuscript Figure/Table: Figure 6 and supplementary
# Release role: MAIN RELEASE WORKFLOW
# Required inputs: deduplicated BAM
# Generated outputs: normalized CREB3L2 bigWig
# Upstream dependencies: deduplicated BAM
# Downstream consumers: Figure 6 coverage panel
# Configuration keys: CUTTAG_FASTQ_DIR, CUTTAG_RESULTS_DIR, REFERENCE_DIR, DATA_PROCESSED_DIR
# Expected environment: POSIX-compatible shell with tools listed in environment/command_line_tools.tsv
# Example run command: bash scripts/04_CREB3L2_regulatory_network/CREB3L2-P.bw.sh
# Reproducibility notes: Analysis parameters and biological logic are unchanged
# from clean_release_v2; only path/configuration wiring and comments were added.
# =============================================================================
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
# shellcheck source=../_shared/paths.sh
source "${SCRIPT_DIR}/../_shared/paths.sh"
load_project_paths
require_input "${CUTTAG_RESULTS_DIR}/CREB3L2-P/CREB3L2-P_bowtie2_sort_picard.bam" "deduplicated CUT&Tag BAM"
"${BAMCOVERAGE_BIN}" --bam "${CUTTAG_RESULTS_DIR}/CREB3L2-P/CREB3L2-P_bowtie2_sort_picard.bam" -o "${CUTTAG_RESULTS_DIR}/CREB3L2-P/CREB3L2-P_bowtie2_sort_picard.bw" --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --extendReads
