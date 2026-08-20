#!/bin/bash
# =============================================================================
# Title: SCENIC.sh
# Purpose: Run either all-cell or SMC-enriched subset pySCENIC
# Manuscript Figure/Table: Figure 4D-E (all-cell), Figure 4G (SMC subset)
# Release role: MAIN RELEASE WORKFLOW
# Required inputs: scope-specific all-cell or SMC-enriched expression loom; human TF list; cisTarget databases
# Generated outputs: scope-specific adjacency CSV, regulon CSV and AUCell loom
# Upstream dependencies: Expression loom; human TF list; cisTarget databases
# Downstream consumers: all-cell SCENIC plots/targets or SMC CREB3L2 activity extraction
# Configuration keys: SCENIC_DB_DIR, SCENIC_TF_LIST, SCENIC_ALL_CELL_LOOM,
# SCENIC_ALL_CELL_PREFIX, SCENIC_SMC_LOOM, SCENIC_SMC_PREFIX
# Expected environment: POSIX-compatible shell with tools listed in environment/command_line_tools.tsv
# Example run command: bash scripts/04_CREB3L2_regulatory_network/SCENIC.sh all-cell
# Reproducibility notes: Both scopes use the retained grn, ctx and aucell
# parameters. Separate loom and prefix keys prevent output collisions.
# =============================================================================
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
# shellcheck source=../_shared/paths.sh
source "${SCRIPT_DIR}/../_shared/paths.sh"
load_project_paths
echo "Start: $(date)"
# Tool command is selected through config/paths.env.
tf=${SCENIC_TF_LIST}
scenic_scope="${1:-all-cell}"
case "${scenic_scope}" in
  all-cell)
    exprDat="${SCENIC_ALL_CELL_LOOM}"
    pfx="${SCENIC_ALL_CELL_PREFIX}"
    ;;
  smc)
    exprDat="${SCENIC_SMC_LOOM}"
    pfx="${SCENIC_SMC_PREFIX}"
    ;;
  *)
    echo "Usage: $0 [all-cell|smc]"
    exit 2
    ;;
esac

dbDir=${SCENIC_DB_DIR}
grnOut=${pfx}.adj.csv
ctxOut=${pfx}.regulons.csv
aucellOut=${pfx}.SCE-auc.loom
require_input "${tf}" "SCENIC TF list"
require_input "${exprDat}" "SCENIC expression loom"
require_input "${dbDir}/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl" "SCENIC motif annotations"
require_input "${dbDir}/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" "SCENIC ranking database"

# run grn
echo "Statt GRN"
"${PYSCENIC_BIN}" grn --num_workers 15 \
  --sparse \
  --method grnboost2 \
  --output $grnOut \
  $exprDat \
  $tf

# run ctx
echo "Run RcisTarget"
"${PYSCENIC_BIN}" ctx --num_workers 15 \
  --output $ctxOut \
  --mask_dropouts \
  --expression_mtx_fname $exprDat \
  --mode "dask_multiprocessing" \
  --min_genes 10 \
  --annotations_fname $dbDir/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
  $grnOut \
  $dbDir/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather

# run aucell
echo "Run AUCell"

"${PYSCENIC_BIN}" aucell --num_workers 15 \
  --output $aucellOut \
  $exprDat \
  $ctxOut

echo "Done: $(date)"
