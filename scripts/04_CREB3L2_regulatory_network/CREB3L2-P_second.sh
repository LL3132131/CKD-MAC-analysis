#!/usr/bin/env bash
# =============================================================================
# Title: CREB3L2-P_second.sh
# Purpose: coverage and annotation
# Manuscript Figure/Table: Figure 6 and supplementary
# Release role: MAIN RELEASE WORKFLOW
# Required inputs: trimmed FASTQ; mm10 index
# Generated outputs: aligned/deduplicated BAM, broadPeak and bedGraph
# Upstream dependencies: trimmed FASTQ; mm10 index
# Downstream consumers: coverage and annotation
# Configuration keys: CUTTAG_FASTQ_DIR, CUTTAG_RESULTS_DIR, REFERENCE_DIR, DATA_PROCESSED_DIR
# Expected environment: POSIX-compatible shell with tools listed in environment/command_line_tools.tsv
# Example run command: bash scripts/04_CREB3L2_regulatory_network/CREB3L2-P_second.sh
# Reproducibility notes: Analysis parameters and biological logic are unchanged
# from clean_release_v2; only path/configuration wiring and comments were added.
# =============================================================================
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
# shellcheck source=../_shared/paths.sh
source "${SCRIPT_DIR}/../_shared/paths.sh"
load_project_paths
require_input "${CUTTAG_RESULTS_DIR}/trimmed/CREB3L2-P_1_val_1.fq.gz" "trimmed read 1"
require_input "${CUTTAG_RESULTS_DIR}/trimmed/CREB3L2-P_2_val_2.fq.gz" "trimmed read 2"
mkdir -p "${CUTTAG_RESULTS_DIR}/CREB3L2-P"
#####要是clean的文件就不需要先trim_galore
"${BOWTIE2_BIN}" -p 100 -x "${REFERENCE_DIR}/mm10/mm10" -1 "${CUTTAG_RESULTS_DIR}/trimmed/CREB3L2-P_1_val_1.fq.gz" -2 "${CUTTAG_RESULTS_DIR}/trimmed/CREB3L2-P_2_val_2.fq.gz" -S "${CUTTAG_RESULTS_DIR}/CREB3L2-P/CREB3L2-P_bowtie2.sam"
"${SAMTOOLS_BIN}" view -@ 10 -b -S "${CUTTAG_RESULTS_DIR}/CREB3L2-P/CREB3L2-P_bowtie2.sam" > "${CUTTAG_RESULTS_DIR}/CREB3L2-P/CREB3L2-P_bowtie2.bam"
"${SAMTOOLS_BIN}" sort -@ 10 "${CUTTAG_RESULTS_DIR}/CREB3L2-P/CREB3L2-P_bowtie2.bam" -o "${CUTTAG_RESULTS_DIR}/CREB3L2-P/CREB3L2-P_bowtie2_sort.bam"
"${SAMTOOLS_BIN}" index "${CUTTAG_RESULTS_DIR}/CREB3L2-P/CREB3L2-P_bowtie2_sort.bam"
####picard2
"${PICARD_BIN}" MarkDuplicates INPUT="${CUTTAG_RESULTS_DIR}/CREB3L2-P/CREB3L2-P_bowtie2_sort.bam" OUTPUT="${CUTTAG_RESULTS_DIR}/CREB3L2-P/CREB3L2-P_bowtie2_sort_picard.bam" METRICS_FILE="${CUTTAG_RESULTS_DIR}/CREB3L2-P/CREB3L2-P_bowtie2_sort_picard.markDup1.metric" REMOVE_DUPLICATES=true ASSUME_SORTED=true CREATE_INDEX=true 
######不同物种参数不同
"${MACS2_BIN}" callpeak -t "${CUTTAG_RESULTS_DIR}/CREB3L2-P/CREB3L2-P_bowtie2_sort_picard.bam" -f BAMPE -g mm --slocal 4000 --keep-dup all --name "${CUTTAG_RESULTS_DIR}/CREB3L2-P/CREB3L2-P" --bdg --broad --nomodel --extsize 200 -p 0.000001
