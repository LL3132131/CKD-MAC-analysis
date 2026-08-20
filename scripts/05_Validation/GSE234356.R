# =============================================================================
# Title: GSE234356.R
# Purpose: Figure 6 target intersection
# Manuscript Figure/Table: Figure 6 and supplementary; Figure 5 is non-code wet-lab validation
# Release role: MAIN RELEASE WORKFLOW
# Required inputs: EXTERNAL_INPUT `GSE234356_DIR/GSE234356_gene_count.txt`
# Generated outputs: External-dataset DEG tables
# Upstream dependencies: GSE234356 count tables and six-sample grouping
# Downstream consumers: Figure 6 target intersection
# Configuration keys: DATA_PROCESSED_DIR, GSE234356_DIR, VALIDATION_RESULTS_DIR
# Expected environment: R; package versions are listed in environment/r_packages.tsv
# Example run command: Rscript scripts/05_Validation/GSE234356.R
# Reproducibility notes: Analysis parameters and biological logic are unchanged
# from clean_release_v2; only path/configuration wiring and comments were added.
# =============================================================================
source("scripts/_shared/paths.R")
PATHS <- load_project_paths()

library(GEOquery)
library(stringr)
library(dplyr)
library(tidyr)
library(dplyr)
library(tidyr)
library(limma)  # 用于归一化
library(readxl)             # 加载包
library(data.table) #多核读取文件
library(DESeq2)

output_dir <- PATHS$VALIDATION_RESULTS_DIR
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
options(timeout = 600)
Sys.setenv("VROOM_CONNECTION_SIZE" = 262144)  # 262144字节 = 256KB

a1 <- fread(configured_path(PATHS, "GSE234356_DIR", "GSE234356_gene_count.txt"),
              header = T,data.table = F)#载入counts，第一列设置为列名
colnames(a1)
################
dup_flag <- duplicated(a1$gene_name)  # 返回逻辑向量：TRUE=重复项，FALSE=首次出现项
a1_unique <- a1[!dup_flag, ]  # 保留gene_name首次出现的行（无重复）
counts <- a1_unique[,2:7] #截取样本基因表达量的counts部分作为counts 
rownames(counts) <- a1_unique$gene_name

sample_info <- data.frame(
  sample = colnames(counts),
  condition = c(rep("Calcification", 3), rep("Dendrophenol", 3)),  # 样本分组
  stringsAsFactors = FALSE
)
rownames(sample_info) <- sample_info$sample  # 行名与count_data列名一致
any(duplicated(a1_unique$gene_name))  # 应输出FALSE（无重复）
########
sample_info$condition <- factor(
  sample_info$condition,
  levels = c("Dendrophenol", "Calcification")  # 调整水平顺序
)
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = sample_info,
  design = ~ condition  # 实验设计：比较condition（Calcification vs Dendrophenol）
)
dds <- DESeq(dds)  # 标准化+差异分析
res <- results(dds)  # 返回包含log2FoldChange、pvalue等的结果数据框

res[which(rownames(res) == 'CREB3L2'),]
############################################
library(tibble)
library(dplyr)
# ----------------------
# 表格1：log2FC < -0.586 且 pvalue < 0.05
# ----------------------
table1 <- res %>%
  as.data.frame() %>%  # 将DESeq2结果转换为数据框
  rownames_to_column("gene_name") %>%  # 将行名（基因名）转为列
  filter(
    log2FoldChange > 0.586,  # log2FC < -0.586（对应fold change < 0.667）
    pvalue < 0.05,            # 原始p值 < 0.05
    !is.na(pvalue)            # 排除p值为NA的基因
  ) %>%
  arrange(log2FoldChange) %>%  # 按log2FC升序排列（下调程度最大的在前）
  select(gene_name, log2FoldChange, pvalue, padj)  # 保留关键列
# ----------------------
# 表格2：log2FC < -1 且 pvalue < 0.05
# ----------------------
table2 <- res %>%
  as.data.frame() %>%
  rownames_to_column("gene_name") %>%
  filter(
    log2FoldChange > 1,      # log2FC < -1（对应fold change < 0.5）
    pvalue < 0.05,
    !is.na(pvalue)
  ) %>%
  arrange(log2FoldChange) %>%
  select(gene_name, log2FoldChange, pvalue, padj)
# ----------------------
# 保存表格为CSV文件（保存在当前工作目录）
# ----------------------
write.csv(table1, "log2FC_more_than_0.586_pval_0.05.csv", row.names = FALSE)
write.csv(table2, "log2FC_more_than_1_pval_0.05.csv", row.names = FALSE)


