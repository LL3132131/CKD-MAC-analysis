# =============================================================================
# Title: 03-2_SCENIC_plot.R
# Purpose: target integration
# Manuscript Figure/Table: Figure 4 and supplementary
# Release role: MAIN RELEASE WORKFLOW
# Required inputs: all-cell AUCell loom and `<SCENIC_ALL_CELL_PREFIX>.regulons.csv`
# Generated outputs: CREB3L2 regulon target data
# Upstream dependencies: the single manuscript all-cell SCENIC run
# Downstream consumers: target integration
# Configuration keys: DATA_PROCESSED_DIR, SCENIC_DB_DIR, SCENIC_TF_LIST, SCENIC_RESULTS_DIR
# Expected environment: R; package versions are listed in environment/r_packages.tsv
# Example run command: Rscript scripts/04_CREB3L2_regulatory_network/03-2_SCENIC_plot.R
# Reproducibility notes: CREB3L2 targets are extracted from the all-cell
# regulon output; the independent SMC-enriched run supplies Figure 4G activity.
# =============================================================================
source("scripts/_shared/paths.R")
PATHS <- load_project_paths()

# Historical environment note: activate a compatible R environment before use.
##可视化
library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)#
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(pheatmap)
library(ggplot2)
library(stringr)
library(pheatmap)

setwd(configured_path(PATHS, "SCENIC_RESULTS_DIR", "00.allCell", "7scenic"))

sce_SCENIC <- open_loom(paste0(PATHS$SCENIC_ALL_CELL_PREFIX, ".SCE-auc.loom"))
regulons_incidMat <- get_regulons(sce_SCENIC, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
class(regulons)

regulonAUC <- get_regulons_AUC(sce_SCENIC, column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(sce_SCENIC)

################RSS分析，查看细胞类型特异性转录因子，需要先加载seurat对象，提取metadata信息，并进行分析！默认是点图！
human_data <- readRDS(configured_path(PATHS, "DATA_PROCESSED_DIR", "7data_umap.rds"))
required_metadata <- c("sample_name", "calcification_group", "Celltype")
missing_metadata <- setdiff(required_metadata, colnames(human_data@meta.data))
if (length(missing_metadata) > 0) {
  stop(
    "7data_umap.rds is missing required metadata: ",
    paste(missing_metadata, collapse = ", ")
  )
}

expected_samples <- c(
  "Mild_1", "Mild_2", "Mild_3", "Mild_4",
  "Severe_1", "Severe_2", "Severe_3"
)
observed_samples <- unique(as.character(human_data$sample_name))
if (!setequal(observed_samples, expected_samples)) {
  stop(
    "7data_umap.rds must contain exactly the seven documented MAC samples. ",
    "Observed: ", paste(sort(observed_samples), collapse = ", ")
  )
}

expected_groups <- c("Mild", "Severe")
observed_groups <- unique(as.character(human_data$calcification_group))
if (!setequal(observed_groups, expected_groups)) {
  stop(
    "calcification_group must contain exactly Mild and Severe. ",
    "Observed: ", paste(sort(observed_groups), collapse = ", ")
  )
}

atlas_celltype_levels <- c(
  "Neutrophil",
  "Fibroblast",
  "non-OFB-SMC",
  "Endothelial cell",
  "CD8 T cell",
  "OFB-SMC",
  "CD4 T cell",
  "NK Cell",
  "Pericyte",
  "Mixed myeloid cell",
  "Macrophage",
  "ProMacs",
  "B cell",
  "Eosinophil",
  "Mast cell"
)
observed_celltypes <- unique(as.character(human_data$Celltype))
prohibited_7data_celltypes <- intersect(
  observed_celltypes,
  c("Dendritic cell", "DC", "unKnown")
)
if (length(prohibited_7data_celltypes) > 0) {
  stop(
    "7data_umap.rds contains prohibited seven-sample Celltype labels: ",
    paste(sort(prohibited_7data_celltypes), collapse = ", ")
  )
}
if (!"Pericyte" %in% observed_celltypes) {
  stop("7data_umap.rds is missing the required Pericyte Celltype.")
}
if (!setequal(observed_celltypes, atlas_celltype_levels)) {
  stop(
    "Celltype must contain exactly the documented 15 seven-sample cell types. ",
    "Missing: ",
    paste(setdiff(atlas_celltype_levels, observed_celltypes), collapse = ", "),
    "; unexpected: ",
    paste(setdiff(observed_celltypes, atlas_celltype_levels), collapse = ", ")
  )
}

cellinfo <- human_data@meta.data[, c(
  "Celltype",
  "sample_name",
  "calcification_group",
  "nFeature_RNA",
  "nCount_RNA"
), drop = FALSE]
cellinfo$celltype <- cellinfo$Celltype

cellTypes <-  as.data.frame(subset(cellinfo,select = 'celltype'))
selectedResolution <- "celltype"
sub_regulonAUC <- regulonAUC

rss <- calcRSS(AUC=getAUC(sub_regulonAUC),
               cellAnnotation=cellTypes[colnames(sub_regulonAUC),
                                        selectedResolution])

rss=na.omit(rss)
rssPlot <-
  plotRSS(rss,cluster_columns = FALSE,order_rows = TRUE,thr=0.1,varName = "cellType",col.low = '#330066',col.mid = '#66CC66',col.high = '#FFCC33')
rssPlot
##################################################
library(readr)
library(dplyr)
library(tidyr) 
input_file <- paste0(PATHS$SCENIC_ALL_CELL_PREFIX, ".regulons.csv")
output_file <- configured_path(PATHS, "SCENIC_RESULTS_DIR", "00.allCell", "7scenic", "CREB3L2_target_genes.csv")

target_pair_pattern <- "\\(['\"]?([A-Za-z0-9.-]+)['\"]?,\\s*([+-]?(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:[eE][+-]?[0-9]+)?)\\)"
weight_threshold_pattern <- "weight>([0-9]+(?:\\.[0-9]+)?)%"

# Non-scientific parser self-checks. These fixed strings do not read or analyze
# project data; they only protect gene symbols containing digits, hyphens or dots.
parser_self_test <- c(
  "('MYO1D', 0.12345)",
  "('PLXDC2', 0.456)",
  "('HLA-DRA', 0.77)"
)
parser_self_test_match <- str_match(parser_self_test, target_pair_pattern)
if (
  !identical(parser_self_test_match[, 2], c("MYO1D", "PLXDC2", "HLA-DRA")) ||
  !isTRUE(all.equal(
    as.numeric(parser_self_test_match[, 3]),
    c(0.12345, 0.456, 0.77),
    tolerance = 0
  ))
) {
  stop("Internal SCENIC target-pair parser self-check failed.")
}

regulons <- read_csv(
  input_file,
  skip = 1,                  # 跳过第1行无效表头
  show_col_types = FALSE,
  name_repair = "unique"     # 修复重复列名
)
 
creb3l2_data <- regulons %>% 
  filter(`...1` == "CREB3L2") %>%  # 筛选转录因子为CREB3L2的行
  select(
    TF = `...1`, 
    Module = `...2`, 
    Context,                      # 新增：用于提取权重阈值的列
    TargetGenes
  ) %>% 
  # 从Context列提取纯numeric权重阈值（如"weight>90.0%" → 90.0）
  mutate(
    WeightThreshold = as.numeric(
      str_match(Context, weight_threshold_pattern)[, 2]
    )
  )

if (anyNA(creb3l2_data$WeightThreshold)) {
  stop("CREB3L2 regulon rows contain missing or invalid WeightThreshold values.")
}

# 2. 提取基因名和权重，并按阈值分组
creb3l2_final <- creb3l2_data %>% 
  # 提取所有 ('基因名',权重) 对
  mutate(gene_weight_pairs = str_extract_all(TargetGenes, target_pair_pattern)) %>% 
  unnest(cols = gene_weight_pairs) %>% 
  # 使用同一个带捕获组的表达式分别提取完整gene symbol和逗号后的weight
  mutate(
    gene = str_match(gene_weight_pairs, target_pair_pattern)[, 2],
    weight = as.numeric(
      str_match(gene_weight_pairs, target_pair_pattern)[, 3]
    )
  ) %>% 
  # 保留关键列，按模块和阈值分组排序
  select(TF, Module, WeightThreshold, gene, weight) %>% 
  arrange(Module, desc(WeightThreshold), desc(weight))  # 先按模块，再按阈值（高→低），再按权重

if (nrow(creb3l2_final) == 0) {
  stop("No CREB3L2 target gene-weight pairs were parsed.")
}
if (anyNA(creb3l2_final$gene) || any(trimws(creb3l2_final$gene) == "")) {
  stop("Parsed CREB3L2 target genes contain NA or empty symbols.")
}
if (anyNA(creb3l2_final$weight)) {
  stop("Parsed CREB3L2 target weights contain NA values.")
}
if (!is.numeric(creb3l2_final$WeightThreshold) ||
    anyNA(creb3l2_final$WeightThreshold)) {
  stop("CREB3L2 WeightThreshold must be complete numeric values.")
}

print(creb3l2_final)
write_csv(creb3l2_final, output_file)
