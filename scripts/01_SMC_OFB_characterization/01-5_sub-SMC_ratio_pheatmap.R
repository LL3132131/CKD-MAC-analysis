# =============================================================================
# Title: 01-5_sub-SMC_ratio_pheatmap.R
# Purpose: Figure 2 / supplementary
# Manuscript Figure/Table: Figure 2; supplementary where indicated
# Release role: SUPPLEMENTARY
# Required inputs: Final SMC object; marker/order tables
# Generated outputs: SMC heatmap, composition plots and source table
# Upstream dependencies: Final SMC object; marker/order tables
# Downstream consumers: Figure 2 / supplementary
# Configuration keys: DATA_PROCESSED_DIR, DATA_EXTERNAL_DIR, SMC_RESULTS_DIR, GWAS_DIR, LD_REFERENCE_DIR
# Expected environment: R; package versions are listed in environment/r_packages.tsv
# Example run command: Rscript scripts/01_SMC_OFB_characterization/01-5_sub-SMC_ratio_pheatmap.R
# Reproducibility notes: Analysis parameters and biological logic are unchanged
# from clean_release_v2; only path/configuration wiring and comments were added.
# =============================================================================
source("scripts/_shared/paths.R")
PATHS <- load_project_paths()

library(data.table)
library(dplyr)
library(Seurat)
library(ggplot2)
library(viridis)
library(scop)
library(pheatmap)

output_dir <- configured_path(PATHS, "SMC_RESULTS_DIR", "05-sub_SMC-ratio")
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
#obj <- readRDS(configured_path(PATHS, "DATA_PROCESSED_DIR", "smc_subclusters_initial_seurat.rds"))
obj <- readRDS(configured_path(PATHS, "DATA_PROCESSED_DIR", "smc_subclusters_final_seurat.rds"))
DefaultAssay(obj) <- 'RNA'
Idents(obj) <- obj$Celltype
difgene <-  FindAllMarkers(obj, test.use = "wilcox",only.pos=T)
head(difgene)
write.csv(difgene,paste0('all-gene-2','.csv'),quote = FALSE)
#difgene <- read.csv(all-gene-2)
co <- c('#3361A5','#3164AB','#3067B1','#2F6AB7','#2E6EBE','#2C71C4','#2B74CA','#2A78D1','#297BD7','#287EDD','#2682E4','#2585EA','#2488F0','#238CF3','#218FF4','#2092F5','#1F96F6','#1E99F7','#1C9CF8','#1B9FF9','#1AA3FA','#18A6FB','#17A9FC','#16ADFD','#14B0FE','#16B3FE','#1FB5FD','#29B7FC','#32BAFA','#3BBCF9','#45BEF8','#4EC0F6','#57C2F5','#61C5F4','#6AC7F3','#74C9F1','#7DCBF0','#86CDEF','#8CCEED','#90CFEC','#95CFEA','#99D0E9','#9ED0E7','#A3D1E5','#A7D1E4','#ACD2E2','#B0D3E1','#B5D3DF','#BAD4DE','#BED4DC','#C2D4D9','#C5D4D3','#C9D4CE','#CCD4C8','#CFD4C2','#D3D4BD','#D6D3B7','#D9D3B2','#DDD3AC','#E0D3A7','#E3D3A1','#E7D39B','#EAD295','#EBD08B','#EDCD81','#EECA77','#F0C86D','#F1C563','#F3C359','#F4C04F','#F6BD44','#F8BB3A','#F9B830','#FBB626','#FCB31C','#FBAA1A','#F99F1C','#F7941D','#F5891E','#F37E20','#F17321','#EF6822','#ED5D24','#EB5225','#E94726','#E73B27','#E53029','#E22929','#DC2828','#D72727','#D22626','#CD2525','#C72424','#C22323','#BD2222','#B82121','#B22020','#AD1F1F','#A81E1E','#A31D1D',"#800026","#800026","#800026","#800026","#800026","#800026","#800026")
cluster.averages <- AverageExpression(obj)
data = data.frame(cluster.averages$RNA)
marker = read.csv(configured_path(PATHS, "DATA_EXTERNAL_DIR", "manuscript_inputs", "figue2", "pheatmap", "pheatmap-marker.csv"),header = F)
data$V1 = rownames(data)
data = merge(marker,data,by = 'V1',sort =F)
data = as.data.frame(data)
rownames(data) = data$V1
data = data[,9:length(data)]###删除V1-V8

p = pheatmap(data,cluster_cols = T,cluster_rows = T,scale = "row",color =co)
p
mat_cluster <- data[p$tree_row$order, p$tree_col$order]
#write.csv(mat_cluster,'order.csv',quote = F)
#####排序聚类
or = read.csv(configured_path(PATHS, "DATA_EXTERNAL_DIR", "manuscript_inputs", "figue2", "pheatmap", "pheatmap-ordergene.csv"))
l <- colnames(or[,2:length(colnames(or))])
dat = data[,l]#l <- gsub("_", ".", l),seuratV5
p = pheatmap(dat,cluster_cols = F,cluster_rows = T,scale = "row",color =co)
p
mat_cluster <- dat[p$tree_row$order,]
#write.csv(mat_cluster,'order1.csv',quote = F)
or2 = read.csv(configured_path(PATHS, "DATA_EXTERNAL_DIR", "manuscript_inputs", "figue2", "pheatmap", "order2.csv"))
da = dat[or2$X,]
p = pheatmap(da,
    cluster_cols = F,
    cluster_rows = F,
    scale = "row",
    color = colorRampPalette(c('MediumTurquoise ','white','OrangeRed'))(100),
    angle_col = 45,display_numbers = F,fontsize_col = 10,gaps_col = NULL)
p
ggsave(plot = p,filename = 'Celltype_pheatmap.pdf',width = 5.8,height = 20.8)
ggsave(plot = p,filename = 'Celltype_pheatmap.png',width = 5.8,height = 20.8,dpi = 300)
########################################
library(grid)  # 加载grid包，提供gpar函数
m <- marker$V1
m<- m[-1]
m <- unique(m)
cols = c('#B2DF8A','#1F78B4','#FDBF6F','#FF7F00','#A6CEE3','#33A02C')

ht8 <- GroupHeatmap(obj,
  features = m, group.by = "Celltype",group_palcolor = cols,
  add_reticle = TRUE, #heatmap_palette = "viridis",
  nlabel = 0, show_row_names = TRUE,
  ht_params = list(row_gap = unit(0, "mm"), row_names_gp = gpar(fontsize = 10,fontfamily = "Arial",fontface = "bold"))
)
p <- ht8$plot
ggsave(p, file ='scp_Celltype_pheatmap.png',width = 5,height = 9,dpi = 900)
