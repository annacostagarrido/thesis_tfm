
############################
#### SingleR annotation ####
############################


# Packages ####
library(readr)
library(readxl)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(SingleR)
######


# SingleR PBMCs - Monaco with all cell types ####
setwd("C:/Users/Usuario/Desktop/master/tfm/PBMC_data")

load("pbmc_processed_all.rda")
ref <- celldex::MonacoImmuneData()

# Per each cell
start_time <- Sys.time()
pred_clust <- SingleR(method = "single",
                      test=as.SingleCellExperiment(pbmcAll),
                      ref=ref,
                      labels=ref$label.main)
end_time <- Sys.time()
time_cells <- end_time - start_time

write.csv(pred_clust$pruned.labels,
          "./results/SingleR_eachCell_Monaco_All.csv",
          row.names = F)
write.csv(as.numeric(time_cells),
          "./results/SingleR_time_eachCell_Monaco_All.csv",
          row.names = F)


# Per clusters
start_time <- Sys.time()
pred_clust <- SingleR(test=as.SingleCellExperiment(pbmcAll), ref=ref,
                      clusters = pbmcAll$seurat_clusters,
                      labels=ref$label.main)
end_time <- Sys.time()
time_clusters <- end_time - start_time


pred <- pbmcAll$seurat_clusters
pred <- as.numeric(pred)

pred[pred == 1] <- "T cells"
pred[pred == 2] <- "Monocytes"
pred[pred == 3] <- "NK cells"
pred[pred == 4] <- "CD4+ T cells"
pred[pred == 5] <- "T cells"
pred[pred == 6] <- "B cells"
pred[pred == 7] <- "T cells"
pred[pred == 8] <- "CD8+ T cells"
pred[pred == 9] <- "B cells"
pred[pred == 10] <- "T cells"
pred[pred == 11] <- "Monocytes"
pred[pred == 12] <- "Dendritic cells"
pred[pred == 13] <- NA
pred[pred == 14] <- "NK cells"
pred[pred == 15] <- "Dendritic cells"
pred[pred == 16] <- NA
pred[pred == 17] <- "Progenitors"

write.csv(pred,
          file = "./results/SingleR_clusters_Monaco_scaled_All.csv",
          row.names = F)
write.csv(as.numeric(time_clusters),
          file = "./results/SingleR_time_clusters_Monaco_scaled_All.csv",
          row.names = F)





# SingleR PBMCs - Monaco with equal cell types ####
load("pbmc_processed_select.rda")
load("refMonacoSelect.rda")

# Per each cell

start_time <- Sys.time()
pred_clust <- SingleR(method = "single",
                      test=as.SingleCellExperiment(pbmcSelect),
                      ref=refMonaco,
                      labels=refMonaco$label.main)
end_time <- Sys.time()
time_cells <- end_time - start_time

write.csv(pred_clust$pruned.labels, "./results/SingleR_eachCell_Monaco_Select.csv", row.names = F)
write.csv(as.numeric(time_cells), "./results/SingleR_time_eachCell_Monaco_Select.csv", row.names = F)


# Per clusters
start <- Sys.time()
pred_clust <- SingleR(test=as.SingleCellExperiment(pbmcSelect),
                      ref=refMonaco,
                      clusters = pbmcSelect$seurat_clusters,
                      labels=refMonaco$label.main)
end_time <- Sys.time()
time_clusters <- end_time - start_time


pred <- pbmcSelect$seurat_clusters
pred <- as.numeric(pred)

pred[pred == 1] <- "CD4+ T cells"
pred[pred == 2] <- "Monocytes"
pred[pred == 3] <- "NK cells"
pred[pred == 4] <- "B cells"
pred[pred == 5] <- "CD8+ T cells"
pred[pred == 6] <- "Monocytes"
pred[pred == 7] <- "B cells"
pred[pred == 8] <- "CD8+ T cells"
pred[pred == 9] <- "NK cells"
pred[pred == 10] <- "Monocytes"
pred[pred == 11] <- "Dendritic cells"
pred[pred == 12] <- "Monocytes"
pred[pred == 13] <- "NK cells"

write.csv(pred, file = "./results/SingleR_clusters_Monaco_scaled_Select.csv", row.names = F)
write.csv(as.numeric(time_clusters), "./results/SingleR_time_clusters_Monaco_Select.csv", row.names = F)





# SingleR MCA - ImmGen with all cell types ####

setwd("C:/Users/Usuario/Desktop/master/tfm/HCL_data")
load("mca_processed_all.rda")
ref <- celldex::ImmGenData()

# Per each cell
start_time <- Sys.time()
pred_clust <- SingleR(method = "single",
                      test=as.SingleCellExperiment(mcaAll),
                      ref=ref,
                      labels=ref$label.main)
end_time <- Sys.time()
time_cells <- end_time - start_time

write.csv(pred_clust$pruned.labels,
          "./results/SingleR_eachCell_ImmGen_scaled_All.csv",
          row.names = F)
write.csv(as.numeric(time_cells*60),
          "./results/SingleR_time_eachCell_ImmGen_scaled_All.csv",
          row.names = F)


# Per clusters

start_time <- Sys.time()
pred_clust <- SingleR(test=as.SingleCellExperiment(mcaAll), ref=ref,
                      clusters = mcaAll$seurat_clusters,
                      labels=ref$label.main)
end_time <- Sys.time()
time_clusters <- end_time - start_time

pred <- mcaAll$seurat_clusters
pred <- as.numeric(pred)

pred[pred == 1] <- "T cells"
pred[pred == 2] <- "B cells"
pred[pred == 3] <- "Monocytes"
pred[pred == 4] <- "Neutrophils"
pred[pred == 5] <- "Neutrophils"
pred[pred == 6] <- "Stem cell"
pred[pred == 7] <- "Neutrophils"
pred[pred == 8] <- "T cells"
pred[pred == 9] <- "NK cells"
pred[pred == 10] <- "Neutrophils"
pred[pred == 11] <- "Monocytes"
pred[pred == 12] <- "B cells, pro"
pred[pred == 13] <- "T cells"
pred[pred == 14] <- "B cells"
pred[pred == 15] <- "Stem cells"
pred[pred == 16] <- "B cells"
pred[pred == 17] <- "Monocytes"
pred[pred == 18] <- "Monocytes"
pred[pred == 19] <- "Monocytes"
pred[pred == 20] <- "Basophils"
pred[pred == 21] <- "B cells"

write.csv(pred,
          "./results/SingleR_clusters_ImmGen_scaled_All.csv",
          row.names = F)
write.csv(as.numeric(time_clusters),
          "./results/SingleR_time_clusters_ImmGen_scaled_All.csv",
          row.names = F)


# SingleR MCA - ImmGen with with equal cell types ####


load("mca_processed_select.rda")
load("refImmScaledSelect.rda")


# Per each cell

start_time <- Sys.time()
pred_clust <- SingleR(method = "single",
                      test=as.SingleCellExperiment(mcaSelect),
                      ref=refSelect,
                      labels=refSelect$label.main)
end_time <- Sys.time()
time_cells <- end_time - start_time

write.csv(pred_clust$pruned.labels,
          "./results/SingleR_eachCell_ImmGen_scaled_Select.csv",
          row.names = F)
write.csv(as.numeric(time_cells),
          "./results/SingleR_time_eachCell_ImmGen_scaled_Select.csv",
          row.names = F)


# Per clusters

start_time <- Sys.time()
pred_clust <- SingleR(test=as.SingleCellExperiment(mcaSelect), ref=refSelect,
                      clusters = mcaSelect$seurat_clusters,
                      labels=refSelect$label.main)
end_time <- Sys.time()
time_clusters <- end_time - start_time

pred <- mcaSelect$seurat_clusters
pred <- as.numeric(pred)

pred[pred == 1] <- "B cells"
pred[pred == 2] <- "T cells"
pred[pred == 3] <- "Neutrophils"
pred[pred == 4] <- "T cells"
pred[pred == 5] <- "Neutrophils"
pred[pred == 6] <- "T cells"
pred[pred == 7] <- "Neutrophils"
pred[pred == 8] <- "Neutrophils"
pred[pred == 9] <- "Neutrophils"
pred[pred == 10] <- "Neutrophils"
pred[pred == 11] <- "B cells"
pred[pred == 12] <- "T cells"
pred[pred == 13] <- "B cells"
pred[pred == 14] <- "Neutrophils"
pred[pred == 15] <- "Neutrophils"

write.csv(pred_clust$pruned.labels,
          "./results/SingleR_clusters_ImmGen_scaled_Select.csv",
          row.names = F)
write.csv(as.numeric(time_clusters),
          "./results/SingleR_time_clusters_ImmGen_scaled_Select.csv",
          row.names = F)
