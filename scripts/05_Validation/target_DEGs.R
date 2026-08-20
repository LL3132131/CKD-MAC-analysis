# =============================================================================
# Title: target_DEGs.R
# Purpose: Figure 6 / supplementary
# Manuscript Figure/Table: Figure 6 and supplementary; Figure 5 is non-code wet-lab validation
# Release role: MAIN RELEASE WORKFLOW
# Required inputs: GSE234356 expression; final SMC object; atlas object
# Generated outputs: Target-expression plots/tables
# Upstream dependencies: GSE234356 expression; final SMC object; atlas object
# Downstream consumers: Figure 6 / supplementary
# Configuration keys: DATA_PROCESSED_DIR, GSE234356_DIR, VALIDATION_RESULTS_DIR
# Expected environment: R; package versions are listed in environment/r_packages.tsv
# Example run command: Rscript scripts/05_Validation/target_DEGs.R
# Reproducibility notes: Analysis parameters and biological logic are unchanged
# from clean_release_v2; only path/configuration wiring and comments were added.
# =============================================================================
source("scripts/_shared/paths.R")
PATHS <- load_project_paths()

library(Seurat)
library(tidyverse)
library(data.table)
library(broom)
library(ggrepel)
library(dplyr)
library(scop)
library(ggsci)
library(ggpubr)  # 用于添加显著性检验
library(SCP)
output_dir <- configured_path(PATHS, "VALIDATION_RESULTS_DIR", "boxplot")
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
#########Bulk seq

fpkm_matrix = fread(configured_path(PATHS, "GSE234356_DIR", "GSE234356_gene_fpkm.csv"))
setnames(fpkm_matrix, "gene_id", "Name")  # 若原代码用"Name"作为基因列名
#################
sample_metadata <- data.frame(
  SUBJID = c("Calcification1", "Calcification2", "Calcification3",
             "Dendrophenol1", "Dendrophenol2", "Dendrophenol3"),
  Diseased = c(1, 1, 1, 0, 0, 0)  # 1=疾病组，0=健康组
)
fpkm_long = fpkm_matrix %>%
  # 步骤1：筛选目标基因
  filter(gene_name %in% target_genes) %>%
  # 步骤2：仅保留基因名和样本列（排除所有注释列，避免类型冲突）
  select(gene_name, Calcification1:Dendrophenol3) %>%
  # 步骤3：转换为长格式（样本名列统一为"Sample"，表达量为"FPKM"）
  pivot_longer(
    cols = -gene_name,  # 排除gene_name列，仅转换样本列
    names_to = "Sample",
    values_to = "FPKM"
  ) %>%
  # 步骤4：合并样本分组信息（显式指定列名映射：Sample -> SUBJID）
  left_join(sample_metadata, by = c("Sample" = "SUBJID")) %>%
  # 步骤5：创建分组标签（健康组在前，疾病组在后）
  mutate(Disease = factor(ifelse(Diseased == 1, "Diseased", "Treatment"), 
                          levels = c("Treatment", "Diseased")))
p = fpkm_long %>%
  ggplot(aes(x = Disease, y = log2(FPKM + 1), color = Disease)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(size = 2, position = position_jitter(width = 0.2)) +
  # 新增：组间差异显著性检验（Wilcoxon秩和检验，适用于小样本）
  stat_compare_means(
    comparisons = list(c("Treatment", "Diseased")),  # 明确比较健康组vs疾病组
    method = "t.test",  # 非参数检验（样本量n=3，推荐用Wilcoxon）
    label = "p.signif",      # 显示星号（ns:不显著, *:p<0.05, **:p<0.01, ***:p<0.001）
    label.x = 1,           # 星号水平位置（1.5对应x轴中间）
    size = 3                 # 星号大小
  ) +
  facet_wrap(~gene_name, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c("Treatment" = "#00468BFF", "Diseased" = "#ED0000FF")) +
  labs(x = "Group", y = "log2(FPKM + 1)", title = "Differential Expression with Significance") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    strip.text = element_text(size = 12, face = "italic"),
    legend.position = "none"
  )
ggsave(p, file= 'GSE234356_target.png')

#############scRNA seq
obj <- readRDS(configured_path(PATHS, "DATA_PROCESSED_DIR", "smc_subclusters_final_seurat.rds"))
required_smc_metadata <- c("Celltype", "sample_name", "calcification_group")
missing_smc_metadata <- setdiff(required_smc_metadata, colnames(obj@meta.data))
if (length(missing_smc_metadata) > 0) {
  stop(
    "smc_subclusters_final_seurat.rds is missing required metadata: ",
    paste(missing_smc_metadata, collapse = ", ")
  )
}
DefaultAssay(obj) <- 'RNA'
Idents(obj) <- "Celltype"
sub <- subset(obj, idents = "OFB-SMC")
Idents(sub) <- "calcification_group"
genes <- c("MYO1D", "PKD1", "SMAD7", "PLXDC2", "ANXA4", "AIG1", "PGRMC2", "SLC17A5", "TMEM50B")

p <- FeatureStatPlot(sub, stat.by = genes,group.by = c("calcification_group"), plot_type = "box",comparisons = list(c("Mild","Severe")),
    pairwise_method = "t.test",sig_label = "p.signif", sig_labelsize = 3
)
ggsave(p,file = 'targer_boxplot.pdf',dpi = 600)

for (i in genes){
p <- FeatureDimPlot(obj,features  = i,reduction = "UMAP",theme_use = "theme_blank")
ggsave(p,file = paste0('sub_targer-',i,'_umap.pdf'),dpi = 600)
}

obj <- readRDS(configured_path(PATHS, "DATA_PROCESSED_DIR", "7data_umap.rds"))
for (i in genes){
 p <- FeatureDimPlot(obj,features  = i,reduction = "UMAP",theme_use = "theme_blank")
 ggsave(p,file = paste0('all_targer-',i,'_umap.pdf'),dpi = 600)
}
