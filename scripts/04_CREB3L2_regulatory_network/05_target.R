# =============================================================================
# Title: 05_target.R
# Purpose: Figure 6
# Manuscript Figure/Table: Figure 6 and supplementary
# Release role: MAIN RELEASE WORKFLOW
# Required inputs: Atlas object; EXTERNAL_INPUT curated target/intersection tables
# Generated outputs: OFB-SMC DEG files, pathway table and four-set Venn plots
# Upstream dependencies: Atlas object; SCENIC/CUT&Tag/external target tables
# Downstream consumers: Figure 6
# Configuration keys: CUTTAG_FASTQ_DIR, CUTTAG_RESULTS_DIR, REFERENCE_DIR, DATA_PROCESSED_DIR
# Expected environment: R; package versions are listed in environment/r_packages.tsv
# Example run command: Rscript scripts/04_CREB3L2_regulatory_network/05_target.R
# Reproducibility notes: Analysis parameters and biological logic are unchanged
# from clean_release_v2; only path/configuration wiring and comments were added.
# =============================================================================
source("scripts/_shared/paths.R")
PATHS <- load_project_paths()

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggrepel)
library(cowplot)
library(SCP)
library(scop)

output_dir <- configured_path(PATHS, "CUTTAG_RESULTS_DIR", "target_integration")
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
figure6_output_dir <- configured_path(
  PATHS,
  "RESULTS_DIR",
  "figures",
  "Figure6"
)
dir.create(
  figure6_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
if (!dir.exists(figure6_output_dir)) {
  stop(
    "Failed to create output directory: ",
    figure6_output_dir
  )
}
setwd(output_dir)
obj = readRDS(configured_path(PATHS, "DATA_PROCESSED_DIR", "7data_umap.rds"))
required_metadata <- c("sample_name", "calcification_group", "Celltype")
missing_metadata <- setdiff(required_metadata, colnames(obj@meta.data))
if (length(missing_metadata) > 0) {
  stop(
    "7data_umap.rds is missing required metadata: ",
    paste(missing_metadata, collapse = ", ")
  )
}
if (!"RNA" %in% names(obj@assays)) {
  stop("7data_umap.rds is missing the required RNA assay.")
}
DefaultAssay(obj) <- "RNA"

observed_groups <- unique(as.character(obj$calcification_group))
expected_groups <- c("Mild", "Severe")
if (!setequal(observed_groups, expected_groups)) {
  stop(
    "calcification_group must contain exactly Mild and Severe. ",
    "Observed: ", paste(sort(observed_groups), collapse = ", ")
  )
}

Idents(obj) <- "Celltype"
ofb <- subset(obj, idents = "OFB-SMC")
DefaultAssay(ofb) <- "RNA"
Idents(ofb) <- "calcification_group"

Idents(obj) <- "Celltype"
ofb_smc_markers <- FindMarkers(obj, ident.1 = "OFB-SMC",  logfc.threshold = 0.5, min.pct = 0.1,min.diff.pct = 0.1,  test.use = "wilcox")

ofb_smc_markers2 <- FindMarkers(ofb, ident.1 = "Severe",test.use = "wilcox",only.pos = T,logfc.threshold = 0.58)

write.csv(ofb_smc_markers,file = 'OFB-SMC_degs.csv')
write.csv(ofb_smc_markers2,file = 'OFB-SMC-Severe-VS-Mild_degs.csv')

#########################GO
library("tidyverse")
library ("org.Hs.eg.db")
library("enrichplot")
library("GOplot")
library("plyr")
library("clusterProfiler")
library("DOSE")
diff1 <- read.csv(configured_path(PATHS, "DATA_EXTERNAL_DIR", "manuscript_inputs", "figue6", "CREB3L2_target.csv"),header = T)
gene.df1 <- bitr(diff1$SYMBOL,
                fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)
head(gene.df1)
Enrich <- gene.df1 %>%
  inner_join(diff1, by = "SYMBOL")
head(Enrich)

gene <- Enrich$ENTREZID
ego <- enrichGO(gene = gene,
                OrgDb = org.Hs.eg.db,
                pvalueCutoff =0.05,
                qvalueCutoff = 0.05,
                ont="all",
                readable =T)
head(ego)
write.csv(ego,"Cuttag+scenic+degs_pathway.csv",row.names=FALSE)
barplot(ego, showCategory = 20,color = "pvalue")
##############################韦恩图
library(ggVennDiagram)
library(ggplot2)
library(dplyr)
gene <- read.csv(configured_path(PATHS, "DATA_EXTERNAL_DIR", "manuscript_inputs", "figue6", "4data.csv"),header = T)
gene_list <- list(
  Cuttag = gene$Cuttag,
  pySCENIC = gene$pySCENIC,
  OFB_SMC.DEGs = gene$OFB_SMC.DEGs,
  GSE234356 = gene$GSE234356
)

# 绘制韦恩图
p <- ggVennDiagram(gene_list,label_alpha=0,label = "count", set_color = c('#322b5d','#a31975','#e8685b','#f8ca72')) + scale_fill_gradient(low="white",high = "white")
  #scale_fill_manual(c("#B2DF8A","#1F78B4","#FDBF6F","#FF7F00")) + 
print(p)
ggsave(configured_path(PATHS, "RESULTS_DIR", "figures", "Figure6", "venn_diagram_4columns.pdf"), p, width = 8, height = 7, dpi = 600)

library(ggvenn)
p_ggvenn <- ggvenn(
  data = gene_list,         # 数据列表
  columns = NULL,           # 对选中的列名绘图，最多选择4个，NULL为默认全选
  show_elements = F,        # 当为TRUE时，显示具体的交集情况，而不是交集个数
  label_sep = "\n",         # 当show_elements = T时生效，分隔符 \n 表示的是回车的意思
  show_percentage = T,      # 显示每一组的百分比
  digits = 1,               # 百分比的小数点位数
  fill_color = c("#E41A1C", "#1E90FF", "#FF8C00", "#80FF00"), # 填充颜色
  fill_alpha = 0.5,         # 填充透明度
  stroke_color = "white",   # 边缘颜色
  stroke_alpha = 0.5,       # 边缘透明度
  stroke_size = 0.5,        # 边缘粗细
  stroke_linetype = "solid", # 边缘线条 # 实线：solid  虚线：twodash longdash 点：dotdash dotted dashed  无：blank
  set_name_color = "black", # 组名颜色
  set_name_size = 6,        # 组名大小
  text_color = "black",     # 交集个数颜色
  text_size = 4             # 交集个数文字大小
)
# 保存图形（可选）
ggsave(
  filename = "venn_diagram_4columns.png",
  plot = p_ggvenn,
  width = 8,
  height = 7,
  dpi = 300
)
