# =============================================================================
# Title: 04_CellChat-All.R
# Purpose: Figure 3 and supplementary
# Manuscript Figure/Table: Figure 3 and supplementary
# Release role: MAIN RELEASE WORKFLOW
# Required inputs: `7data_umap.rds`
# Generated outputs: CellChat RDS; interaction/pathway/bubble plots; pathway table
# Upstream dependencies: `7data_umap.rds`
# Downstream consumers: Figure 3 and supplementary
# Configuration keys: DATA_PROCESSED_DIR, CELLCHAT_RESULTS_DIR
# Expected environment: R; package versions are listed in environment/r_packages.tsv
# Example run command: Rscript scripts/02_Cell_communication/04_CellChat-All.R
# Reproducibility notes: Analysis parameters and biological logic are unchanged
# from clean_release_v2; only path/configuration wiring and comments were added.
# =============================================================================
source("scripts/_shared/paths.R")
PATHS <- load_project_paths()

library(ggplot2)
library(Seurat)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

# Here we load a scRNA-seq data matrix and its associated cell meta data
# This is a combined data from two biological conditions: normal and diseases
obj <- readRDS(configured_path(PATHS, "DATA_PROCESSED_DIR", "7data_umap.rds"))
required_metadata <- c("sample_name", "calcification_group", "Celltype")
missing_metadata <- setdiff(required_metadata, colnames(obj@meta.data))
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
observed_samples <- unique(as.character(obj$sample_name))
if (!setequal(observed_samples, expected_samples)) {
  stop(
    "7data_umap.rds must contain exactly the seven documented MAC samples. ",
    "Observed: ", paste(sort(observed_samples), collapse = ", ")
  )
}

expected_groups <- c("Mild", "Severe")
observed_groups <- unique(as.character(obj$calcification_group))
if (!setequal(observed_groups, expected_groups)) {
  stop(
    "calcification_group must contain exactly Mild and Severe. ",
    "Observed: ", paste(sort(observed_groups), collapse = ", ")
  )
}

major_celltype_levels <- c(
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
observed_celltypes <- unique(as.character(obj$Celltype))
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
if (!setequal(observed_celltypes, major_celltype_levels)) {
  stop(
    "Celltype must contain exactly the documented 15 seven-sample cell types. ",
    "Missing: ",
    paste(setdiff(major_celltype_levels, observed_celltypes), collapse = ", "),
    "; unexpected: ",
    paste(setdiff(observed_celltypes, major_celltype_levels), collapse = ", ")
  )
}

Idents(obj) <- "Celltype"
output_dir <- PATHS$CELLCHAT_RESULTS_DIR
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
figure3_output_dir <- configured_path(
  PATHS,
  "RESULTS_DIR",
  "figures",
  "Figure3"
)
dir.create(
  figure3_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
if (!dir.exists(figure3_output_dir)) {
  stop(
    "Failed to create output directory: ",
    figure3_output_dir
  )
}
setwd(output_dir)
obj1 = obj

data.input <- obj1[["RNA"]]@data
Idents(obj1) <- "Celltype"
labels <- Idents(obj1)
meta <- data.frame(labels = labels, row.names = names(labels))

cellchat <- createCellChat(object = data.input, group.by = "Celltype", meta = obj1@meta.data)
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)

dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact"), key = "annotation")
# set the used database in the object
cellchat@DB <- CellChatDB.use
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
ptm = Sys.time()
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))
# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.human)

ptm = Sys.time()
cellchat <- computeCommunProb(cellchat, type = "triMean")

#Users can filter out the cell-cell communication if there are only few cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

#Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

#Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))
ptm = Sys.time()
groupSize <- as.numeric(table(cellchat@idents))

pdf('Number_of_interactions.pdf',width =6,height = 6)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

pdf('Interaction_weights.pdf',width =6,height = 6)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

pdf('Celltype_weights.pdf',width =12,height = 12)
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()
# 逐个展示细胞亚群的互作信号强度-状图
pdf(paste0('every-subcelltype_barplot.pdf'))
lapply(1:nrow(mat),function(i){
    df = data.frame(cell_type = colnames(mat),weight = as.numeric(mat[i,]))
    p = ggplot(df, aes(x = cell_type, y = weight, fill =cell_type)) +
    geom_bar(stat = "identity",fill = "steelblue",width=0.3)+
    labs(title = rownames(mat)[i], x = 'subcellType', y = "Siganl Strength")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5),axis.text = element_text(size = 4)) + theme_bw()
    return(p)
})
dev.off()

path = unique(cellchat@netP$pathways)
lapply(path,function(pathways.show){
    try({
        par(mfrow=c(1,1))
        pathway_pdf <- paste0(
          "pathways-",
          pathways.show,
          ".pdf"
        )
        pdf(
          pathway_pdf,
          width = 10,
          height = 10
        )
        netVisual_aggregate(cellchat,signaling = pathways.show,layout = 'circle')
        p1 = netVisual_heatmap(cellchat,signaling = pathways.show,color.heatmap = "Reds")
        print(p1)
        p2 = netAnalysis_contribution(cellchat,signaling = pathways.show)
        print(p2)
        extractEnrichedLR(cellchat,signaling = pathways.show,geneLR.return = FALSE)
        cellchat = netAnalysis_computeCentrality(cellchat,slot.name = "netP")
        netAnalysis_signalingRole_network(cellchat,signaling = pathways.show,width=15,height=2.5,font.size=10)
        dev.off()
    },silent =F)
})

p <- netVisual_bubble(cellchat, sources.use = "OFB-SMC", remove.isolate = FALSE)
ggsave(p,file = 'netVisual_bubble.png',width = 6,height = 28)

saveRDS(cellchat,paste0('no_sample_all','_cellchat.rds'))
##################################################
df.net <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.net, "net_pathway.csv")

#需要指定受体细胞和配体细胞
#levels(cellchat@idents)
pathways <- cellchat@netP$pathways
print(pathways)  # 输出通路名称（如"TGFb"、"VEGF"、"NOTCH"等）
path <- c('COLLAGEN','FN1','THBS','VISFATIN','SPP1','JAM','CD99','CypA','PARs','NOTCH','MK','ANGPT','GRN','GALECTIN')
p = netVisual_bubble(cellchat, signaling = path, remove.isolate = T
,targets.use = "OFB-SMC"
)
range(p$data$prob,na.rm = T)
summary(p$data$prob,na.rm = T)

p1 <- p + scale_size_continuous(range = c(4, 8), guide = "none") +   # 调整气泡大小范围
  scale_color_gradientn(
    colours = c("#2760a9", "white", "#e50f20"),  # 定义颜色向量
    values = scales::rescale(c(0.2, 0.35, 0.5)),  # 定义颜色映射的范围
    name = "Commun. Prob." ) +  # 图例标题 
  xlab(label = NULL) + 
  ylab(label = NULL) + 
#  geom_vline(xintercept = seq(1.5, length(unique(df.net$source)) - 0.5, 1)[1:6],lwd = 0.5) + ## 根据 netVisual_bubble 函数的源码，修改格子线的粗细
 # geom_hline(yintercept = seq(1.5, length(unique(df.net$interaction_name_2)) - 0.5, 1)[1:3], lwd = 0.5) + 
  theme(axis.title.x = element_text(size = 16),  # 设置 x 轴标题字体大小
    axis.title.y = element_text(size = 16),  # 设置 y 轴标题字体大小
    axis.text.x = element_text(size = 14),  # 设置 x 轴刻度标签字体大小
    axis.text.y = element_text(size = 14, face = "italic"),   # 设置 y 轴刻度标签字体大小
    panel.border = element_rect(color = "black", fill=NA, size=1),  # 设置四周边框的颜色和粗细
    legend.key.size = unit(0.6, 'cm'),  # 设置图例键的大小
    legend.text = element_text(size = 12),  # 设置图例文本的大小
    legend.title = element_text(size = 13)
  ) 
p1
ggsave(p1,file = configured_path(PATHS, "RESULTS_DIR", "figures", "Figure3", "netVisual_bubble.png"),width = 8,height = 28)
