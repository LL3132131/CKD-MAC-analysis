# =============================================================================
# Title: overlap.R
# Purpose: Record the historical exploratory overlap between CREB3L2 SCENIC
# targets and CUT&Tag-associated genes.
# Manuscript Figure/Table: Provenance only; not a default Figure 6 producer
# Release role: PROVENANCE / EXPLORATORY
# Required inputs: generated CREB3L2 target CSV and generated CUT&Tag peak
# annotation with the explicitly checked historical columns
# Generated outputs: in-memory/console exploratory overlap and mapping objects;
# no formal manuscript handoff file
# Upstream dependencies: `03-2_SCENIC_plot.R` and `anno-end.R`
# Downstream consumers: None in the default public Figure 6 workflow
# Configuration keys: CUTTAG_RESULTS_DIR, SCENIC_RESULTS_DIR
# Expected environment: R; package versions are listed in environment/r_packages.tsv
# Example run command: Optional provenance only; not a default release command
# Reproducibility notes: This script preserves historical exploratory
# cross-species mapping logic. It does not generate the curated external inputs
# used by `05_target.R` and must not be interpreted as a validated automated
# Figure 6 integration workflow.
# =============================================================================
source("scripts/_shared/paths.R")
PATHS <- load_project_paths()

library(biomaRt)
library(dplyr)
library(homologene)

output_dir <- configured_path(PATHS, "CUTTAG_RESULTS_DIR", "integration")
dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
if (!dir.exists(output_dir)) {
  stop(
    "Failed to create output directory: ",
    output_dir
  )
}
setwd(output_dir)
regulon <- read.csv(configured_path(
  PATHS,
  "SCENIC_RESULTS_DIR",
  "00.allCell",
  "7scenic",
  "CREB3L2_target_genes.csv"
))
required_regulon_columns <- c("TF", "gene")
missing_regulon_columns <- setdiff(
  required_regulon_columns,
  colnames(regulon)
)
if (length(missing_regulon_columns) > 0L) {
  stop(
    "CREB3L2 target table is missing required columns: ",
    paste(missing_regulon_columns, collapse = ", ")
  )
}
#########################
creb3l2_targets <- regulon %>%
  filter(TF == "CREB3L2") %>%
  pull(gene)  # 提取靶基因列

cut_tag_annotation_file <- configured_path(
  PATHS,
  "CUTTAG_RESULTS_DIR",
  "CREB3L2-P",
  "CREB3L2-P_peaks.broadPeak_annotation.xls"
)
cut_tag_anno <- read.delim(cut_tag_annotation_file, sep = "\t")
required_cut_tag_columns <- c(
  "Annotation",
  "Distance_to_TSS",
  "Nearest_Gene",
  "Gene_Name"
)
missing_cut_tag_columns <- setdiff(
  required_cut_tag_columns,
  colnames(cut_tag_anno)
)
if (length(missing_cut_tag_columns) > 0L) {
  stop(
    "CUT&Tag annotation table is missing required historical columns: ",
    paste(missing_cut_tag_columns, collapse = ", ")
  )
}
# 提取峰关联的靶基因（优先选择启动子区或近端调控区的基因）
# 筛选条件：Annotation包含"Promoter"或Distance_to_TSS < 2000（2kb内）
creb3l2_bound_genes <- cut_tag_anno %>%
  filter(
    grepl("Promoter", Annotation) | abs(Distance_to_TSS) <= 2000
  ) %>%
  pull(Nearest_Gene)  # 提取最近基因列（需确保基因名为Symbol）

mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")  # 人基因组
# 假设CUT&Tag注释的基因是Ensembl ID，需转为Symbol
creb3l2_bound_genes_symbol <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = creb3l2_bound_genes,
  mart = mart
) %>% pull(external_gene_name)

overlap_genes <- intersect(creb3l2_targets, creb3l2_bound_genes_symbol)
cat("共找到", length(overlap_genes), "个双证据支持的靶基因：\n")
print(overlap_genes)
p_value <- phyper(
  q = length(overlap_genes) - 1,
  m = length(creb3l2_targets),
  n = 20000 - length(creb3l2_targets),
  k = length(creb3l2_bound_genes_symbol),
  lower.tail = FALSE
)
cat("超几何检验p值：", p_value)  # p < 0.05表示重叠显著
##############################################################

mouse_genes <- unique(cut_tag_anno$Gene_Name) 
# 小鼠（NCBI taxonomy ID=10090）→人（ID=9606）同源基因映射
homologs <- homologene(mouse_genes, 10090, 9606)
# 筛选一对一同源基因（避免多对多或无同源的情况）
human_genes <- homologs %>% 
  filter(one2one == TRUE) %>%  # 仅保留严格一对一同源
  pull(human_gene)  # 获取人同源基因列表

