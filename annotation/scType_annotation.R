
###########################
#### scType annotation ####
###########################

# Packages ####
library(dplyr)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(HGNChelper)
library(readxl)

#####


# scType annotation with MCA #####


# LOAD AND CLUSTER THE DATA
setwd("C:/Users/Usuario/Desktop/master/tfm/HCL_data")

# 2) Then, we load the h5Seurat file into a Seurat object
mca <- LoadH5Seurat("MCA1.1_adata_Peripheral.h5seurat")
mca


# normalize data
mca[["percent.mt"]] <- PercentageFeatureSet(mca, pattern = "^MT-")
mca <- NormalizeData(mca, normalization.method = "LogNormalize", scale.factor = 10000)
mca <- FindVariableFeatures(mca, selection.method = "vst", nfeatures = 2000)

# scale and run PCA
mca <- ScaleData(mca, features = rownames(mca))
mca <- RunPCA(mca, features = VariableFeatures(object = mca))

# Check number of PC components (we selected 10 PCs for downstream analysis, based on Elbow plot)
ElbowPlot(mca)

# cluster and visualize
mca <- FindNeighbors(mca, dims = 1:10)
mca <- FindClusters(mca, resolution = 0.8)
mca <- RunUMAP(mca, dims = 1:10)
DimPlot(mca, reduction = "umap")


# CELL TYPE ASSIGNMENT
# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")


# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Immune system" # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)


# get cell-type by cell matrix

start_time <- Sys.time()
es.max = sctype_score(scRNAseqData = mca[["RNA"]]@scale.data, scaled = TRUE,
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)


# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix.
# In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

# merge by cluster
start_time <- Sys.time()
cL_resutls = do.call("rbind", lapply(unique(mca@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(mca@meta.data[mca@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(mca@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
end_time <- Sys.time()
time_mca2 <- end_time - start_time

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])


mca@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,];
  mca@meta.data$customclassif[mca@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}
end_time <- Sys.time()
time_mca <- end_time - start_time

DimPlot(mca, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')


# load libraries
lapply(c("ggraph","igraph","tidyverse", "data.tree"), library, character.only = T)

# prepare edges
cL_resutls=cL_resutls[order(cL_resutls$cluster),]; edges = cL_resutls; edges$type = paste0(edges$type,"_",edges$cluster); edges$cluster = paste0("cluster ", edges$cluster); edges = edges[,c("cluster", "type")]; colnames(edges) = c("from", "to"); rownames(edges) <- NULL

# prepare nodes
nodes_lvl1 = sctype_scores[,c("cluster", "ncells")]; nodes_lvl1$cluster = paste0("cluster ", nodes_lvl1$cluster); nodes_lvl1$Colour = "#f1f1ef"; nodes_lvl1$ord = 1; nodes_lvl1$realname = nodes_lvl1$cluster; nodes_lvl1 = as.data.frame(nodes_lvl1); nodes_lvl2 = c();
ccolss= c("#5f75ae","#92bbb8","#64a841","#e5486e","#de8e06","#eccf5a","#b5aa0f","#e4b680","#7ba39d","#b15928","#ffff99", "#6a3d9a","#cab2d6","#ff7f00","#fdbf6f","#e31a1c","#fb9a99","#33a02c","#b2df8a","#1f78b4","#a6cee3")
for (i in 1:length(unique(cL_resutls$cluster))){
  dt_tmp = cL_resutls[cL_resutls$cluster == unique(cL_resutls$cluster)[i], ]; nodes_lvl2 = rbind(nodes_lvl2, data.frame(cluster = paste0(dt_tmp$type,"_",dt_tmp$cluster), ncells = dt_tmp$scores, Colour = ccolss[i], ord = 2, realname = dt_tmp$type))
}
nodes = rbind(nodes_lvl1, nodes_lvl2); nodes$ncells[nodes$ncells<1] = 1;
files_db = openxlsx::read.xlsx(db_)[,c("cellName","shortName")]; files_db = unique(files_db); nodes = merge(nodes, files_db, all.x = T, all.y = F, by.x = "realname", by.y = "cellName", sort = F)
nodes$shortName[is.na(nodes$shortName)] = nodes$realname[is.na(nodes$shortName)]; nodes = nodes[,c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]

mygraph <- graph_from_data_frame(edges, vertices=nodes)

# Make the graph
gggr<- ggraph(mygraph, layout = 'circlepack', weight=I(ncells)) +
  geom_node_circle(aes(filter=ord==1,fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) + geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +
  theme_void() + geom_node_text(aes(filter=ord==2, label=shortName, colour=I("#ffffff"), fill="white", repel = !1, parse = T, size = I(log(ncells,25)*1.5)))+ geom_node_label(aes(filter=ord==1,  label=shortName, colour=I("#000000"), size = I(3), fill="white", parse = T), repel = !0, segment.linetype="dotted")

scater::multiplot(DimPlot(mca, reduction = "umap", label = TRUE, repel = TRUE, cols = ccolss), gggr, cols = 2)



# load auto-detection function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/auto_detect_tissue_type.R")

# guess a tissue type
tissue_guess = auto_detect_tissue_type(path_to_db_file = db_, scaled = TRUE, seuratObject = mca, assay = "RNA")



# COMPARE ANNOTATION

pred <- mca$seurat_clusters
pred <- as.numeric(pred)

pred[pred == 12] <- "Naive CD8+ T cells"
pred[pred == 15] <- "Pre-B cells"
pred[pred == 1] <- "Naive B cells"
pred[pred == 9] <- "Natural killer  cells"
pred[pred == 5] <- "Naive CD8+ T cells"
pred[pred == 3] <- "Naive CD8+ T cells"
pred[pred == 8] <- "Intermediate monocytes"
pred[pred == 4] <- "Neutrophils"
pred[pred == 16] <- "Macrophages"
pred[pred == 2] <- "Macrophages"
pred[pred == 13] <- "Erythroid-like and erythroid precursor cells"
pred[pred == 6] <- "Neutrophils"
pred[pred == 11] <- "Pro-B cells"
pred[pred == 10] <- "Neutrophils"
pred[pred == 14] <- "Basophils"
pred[pred == 7] <- "Neutrophils"


MCA_cell_info <- read_excel("C:/Users/Usuario/Desktop/master/tfm/HCL_data/MCA1.1_cell_info.xlsx")
MCA_cell_info_Peripherial <- MCA_cell_info[grep("PeripheralBlood",MCA_cell_info$cellnames),]


true <- MCA_cell_info_Peripherial$celltype


table(true, pred)




# scType annotation with PBMCs ####

setwd("C:/Users/Usuario/Desktop/master/tfm/PBMC_data")
load("pbmc_processed_all.rda")
pbmc <- pbmcAll
rm(pbmcAll)

# CELL TYPE ASSIGNMENT
# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")



# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Immune system" # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)


# get cell-type by cell matrix
start_time <- Sys.time()
es.max = sctype_score(scRNAseqData = pbmc[["RNA"]]@scale.data, scaled = TRUE,
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix.
# In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(pbmc@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(pbmc@meta.data[pbmc@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(pbmc@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])


pbmc@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,];
  pbmc@meta.data$customclassif[pbmc@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}
end_time <- Sys.time()
time_pbmc <- end_time - start_time


DimPlot(pbmc, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')


# load libraries
lapply(c("ggraph","igraph","tidyverse", "data.tree"), library, character.only = T)

# prepare edges
cL_resutls=cL_resutls[order(cL_resutls$cluster),]; edges = cL_resutls; edges$type = paste0(edges$type,"_",edges$cluster); edges$cluster = paste0("cluster ", edges$cluster); edges = edges[,c("cluster", "type")]; colnames(edges) = c("from", "to"); rownames(edges) <- NULL

# prepare nodes
nodes_lvl1 = sctype_scores[,c("cluster", "ncells")]; nodes_lvl1$cluster = paste0("cluster ", nodes_lvl1$cluster); nodes_lvl1$Colour = "#f1f1ef"; nodes_lvl1$ord = 1; nodes_lvl1$realname = nodes_lvl1$cluster; nodes_lvl1 = as.data.frame(nodes_lvl1); nodes_lvl2 = c();
ccolss= c("#5f75ae","#92bbb8","#64a841","#e5486e","#de8e06","#eccf5a","#b5aa0f","#e4b680","#7ba39d","#b15928","#ffff99", "#6a3d9a","#cab2d6","#ff7f00","#fdbf6f","#e31a1c","#fb9a99","#33a02c","#b2df8a","#1f78b4","#a6cee3")
for (i in 1:length(unique(cL_resutls$cluster))){
  dt_tmp = cL_resutls[cL_resutls$cluster == unique(cL_resutls$cluster)[i], ]; nodes_lvl2 = rbind(nodes_lvl2, data.frame(cluster = paste0(dt_tmp$type,"_",dt_tmp$cluster), ncells = dt_tmp$scores, Colour = ccolss[i], ord = 2, realname = dt_tmp$type))
}
nodes = rbind(nodes_lvl1, nodes_lvl2); nodes$ncells[nodes$ncells<1] = 1;
files_db = openxlsx::read.xlsx(db_)[,c("cellName","shortName")]; files_db = unique(files_db); nodes = merge(nodes, files_db, all.x = T, all.y = F, by.x = "realname", by.y = "cellName", sort = F)
nodes$shortName[is.na(nodes$shortName)] = nodes$realname[is.na(nodes$shortName)]; nodes = nodes[,c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]

mygraph <- graph_from_data_frame(edges, vertices=nodes)

# Make the graph
gggr<- ggraph(mygraph, layout = 'circlepack', weight=I(ncells)) +
  geom_node_circle(aes(filter=ord==1,fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) + geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +
  theme_void() + geom_node_text(aes(filter=ord==2, label=shortName, colour=I("#ffffff"), fill="white", repel = !1, parse = T, size = I(log(ncells,25)*1.5)))+ geom_node_label(aes(filter=ord==1,  label=shortName, colour=I("#000000"), size = I(3), fill="white", parse = T), repel = !0, segment.linetype="dotted")

scater::multiplot(DimPlot(pbmc, reduction = "umap", label = TRUE, repel = TRUE, cols = ccolss), gggr, cols = 2)


# load auto-detection function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/auto_detect_tissue_type.R")

# guess a tissue type
tissue_guess = auto_detect_tissue_type(path_to_db_file = db_, scaled = T, seuratObject = pbmc, assay = "RNA")

meta <- read_table2("meta.txt")
meta <- meta[-1,]

meta <- meta[grep("^pbmc1_10x_v2",meta$NAME),]
meta <- meta[meta$NAME %in% colnames(GetAssayData(pbmc)),]

pbmc <- pbmc[,meta$NAME]

true <- meta$CellType

table(true)


pred <- pbmc$seurat_clusters
pred <- as.numeric(pred)


pred[pred == 2] <- "Classical Monocytes"
pred[pred == 12] <- "Myeloid Dendritic cells"
pred[pred == 13] <- "Classical Monocytes"
pred[pred == 11] <- "Non-classical monocytes"
pred[pred == 3] <- "Natural killer  cells"
pred[pred == 15] <- "Plasmacytoid Dendritic cells"
pred[pred == 6] <- "Naive B cells"
pred[pred == 9] <- "Naive B cells"
pred[pred == 14] <- "Platelets"
pred[pred == 4] <- "Memory CD4+ T cells"
pred[pred == 5] <- "Effector CD8+ T cells"
pred[pred == 1] <- "Effector CD8+ T cells"
pred[pred == 10] <- "T cells"
pred[pred == 7] <- "Naive CD8+ T cells"
pred[pred == 16] <- "Progenitor cells"
pred[pred == 8] <- "Naive CD4+ T cells"
pred[pred == 17] <- "Basophils"


table(true, pred)

