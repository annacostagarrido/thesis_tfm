
###################
# Processing data #
###################


#### Packages ####

library(Matrix)
library(readr)
library(readxl)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(SingleR)
library(tidyverse)
library(batchelor)
library(scran)
library(SingleCellExperiment)

##################


#### DATA PREPARATION - PBMC ####

setwd("C:/Users/Usuario/Desktop/master/tfm/PBMC_data")

# Import dataset in MM format
pbmc.data <- readMM("counts.umi.txt")
dim(pbmc.data) # 33694 x 44433

# Import the names of cells and gens
cells <- read_csv("cells.umi.new.txt", col_names = FALSE) # 44432
genes <- read_delim("genes.umi.txt", delim = "_",
                    escape_double = FALSE, col_names = FALSE,
                    trim_ws = TRUE) # 33694

colnames(pbmc.data) <- cells$X1
rownames(pbmc.data) <- genes$X2


# Subset 10x Chromium method
pbmc.data <- pbmc.data[,grep("^pbmc1_10x_v2",colnames(pbmc.data))]


# Create Seurat object
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc",
                           min.cells = 3, min.features = 200)



# Obtain the true cell types

meta <- read_table2("meta.txt")
meta <- meta[-1,]

meta <- meta[grep("^pbmc1_10x_v2",meta$NAME),]
meta <- meta[meta$NAME %in% colnames(GetAssayData(pbmc)),]


metaSelect <- filter(meta,
                     CellType == "B" |
                       CellType == "Dendritic" |
                       CellType == "CD14+" |
                       CellType == "CD16+" |
                       CellType == "Natural" |
                       CellType == "CD4+" )




pbmcAll <- pbmc[,meta$NAME]
pbmcSelect <- pbmc[,metaSelect$NAME]



# Processing all cells
pbmcAll <- NormalizeData(pbmcAll)
pbmcAll <- FindVariableFeatures(pbmcAll, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmcAll)
pbmcAll <- ScaleData(pbmcAll, features = all.genes)
pbmcAll <- RunPCA(pbmcAll, features = VariableFeatures(object = pbmcAll))
pbmcAll <- RunUMAP(pbmcAll, reduction="pca", dims=1:20)
pbmcAll <- FindNeighbors(pbmcAll, reduction="pca", dims=1:20)
pbmcAll <- FindClusters(pbmcAll, resolution=0.8)

save(pbmcAll, file = "pbmc_processed_all.rda")
write.csv(GetAssayData(pbmcAll, slot = "scale.data"), file = "pbmc_processed_all.csv", row.names = T)
write.csv(meta$CellType, file = "labels_pbmc_all.csv", row.names = F)


# Processing Select cells
pbmcSelect <- NormalizeData(pbmcSelect)
pbmcSelect <- FindVariableFeatures(pbmcSelect, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmcSelect)
pbmcSelect <- ScaleData(pbmcSelect, features = all.genes)
pbmcSelect <- RunPCA(pbmcSelect, features = VariableFeatures(object = pbmcSelect))
pbmcSelect <- RunUMAP(pbmcSelect, reduction="pca", dims=1:20)
pbmcSelect <- FindNeighbors(pbmcSelect, reduction="pca", dims=1:20)
pbmcSelect <- FindClusters(pbmcSelect, resolution=0.8)


save(pbmcSelect, file = "pbmc_processed_select.rda")
write.csv(GetAssayData(pbmcSelect, slot = "scale.data"), file = "pbmc_processed_select.csv", row.names = T)
write.csv(metaSelect$CellType, file = "labels_pbmc_select.csv", row.names = F)


# ImmGen reference

ref <- celldex::ImmGenData()
table(ref$label.main)


cell_names <- data.frame(ref$label.main, colnames(assay(ref)))

reference_cell_Imm <- cell_names %>% filter(ref.label.main ==  "B cells" |
                                              ref.label.main ==  "T cells" |
                                              ref.label.main ==  "DC" |
                                              ref.label.main ==  "Monocytes" |
                                              ref.label.main ==  "NK cells"

)

refImm <- ref[,reference_cell_Imm$colnames.assay.ref..]

save(refImm, file = "refImmSelect.rda")
write.csv(ScaleData(assay(refImm)), file = "refImmScaledSelect.csv", row.names = T)
write.csv(reference_cell_Imm$ref.label.main, file = "labels_refImm_select.csv", row.names = F)


# Monaco reference

ref <- celldex::MonacoImmuneData()
table(ref$label.main)


cell_names <- data.frame(ref$label.main, colnames(assay(ref)))

reference_cell_Monaco <- cell_names %>% filter(ref.label.main ==  "B cells" |
                                                 ref.label.main ==  "CD4+ T cells" |
                                                 ref.label.main ==  "CD8+ T cells" |
                                                 ref.label.main ==  "NK cells" |
                                                 ref.label.main ==  "Dendritic cells" |
                                                 ref.label.main ==  "Monocytes"

)

refMonaco <- ref[,reference_cell_Monaco$colnames.assay.ref..]

save(refMonaco, file = "refMonacoSelect.rda")
write.csv(ScaleData(assay(refMonaco)), file = "refMonacoScaledSelect.csv", row.names = T)
write.csv(reference_cell_Monaco$ref.label.main, file = "labels_refMonaco_select.csv", row.names = F)


ref <- celldex::ImmGenData()
write.csv(ScaleData(assay(ref)), file = "refImmScaled.csv", row.names = T)
write.csv(ref$label.main, file = "labels_refImm.csv", row.names = F)


ref <- celldex::MonacoImmuneData()
write.csv(ScaleData(assay(ref)), file = "refMonacoScaled.csv", row.names = T)
write.csv(ref$label.main, file = "labels_refMonaco.csv", row.names = F)

#####


#### DATA PROCESSING - MCA ####

setwd("C:/Users/Usuario/Desktop/master/tfm/HCL_data/")

# 1) Convert the AnnData file to an h5Seurat file using the Convert function;
Convert("MCA1.1_adata_Peripheral.h5ad", dest = "h5seurat", overwrite = T)


# 2) Then, we load the h5Seurat file into a Seurat object
mca <- LoadH5Seurat("MCA1.1_adata_Peripheral.h5seurat")


# Processing all cells

mcaAll <- NormalizeData(mca)
mcaAll <- FindVariableFeatures(mcaAll, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(mcaAll)
mcaAll <- ScaleData(mcaAll, features = all.genes)
mcaAll <- RunPCA(mcaAll, features = VariableFeatures(object = mcaAll))
mcaAll <- RunUMAP(mcaAll, reduction="pca", dims=1:20)
mcaAll <- FindNeighbors(mcaAll, reduction="pca", dims=1:20)
mcaAll <- FindClusters(mcaAll, resolution=0.8)

MCA_cell_info <- read_excel("C:/Users/Usuario/Desktop/master/tfm/HCL_data/MCA1.1_cell_info.xlsx")
MCA_cell_info_Peripherial <- MCA_cell_info[grep("PeripheralBlood",MCA_cell_info$cellnames),]

save(mcaAll, file = "mca_processed_all.rda")
write.csv(GetAssayData(mcaAll, slot = "scale.data"), file = "mca_processed_all.csv", row.names = T)
write.csv(MCA_cell_info_Peripherial$celltype, file = "labels_mca_all.csv", row.names = F)



# Processing select cells

test_cells <- filter(MCA_cell_info_Peripherial,
                     celltype == "B cell" |
                       celltype == "Cd4+ T cell" |
                       celltype == "Cd8+ T cell" |
                       celltype == "Dendritic cell" |
                       celltype == "Macrophage" |
                       celltype == "Mast cell" |
                       celltype == "Neutrophil" |
                       celltype == "T cell")


mcaSelect <- mca[,colnames(mca@assays$RNA@data) %in% test_cells$cellnames]

mcaSelect <- NormalizeData(mcaSelect)
mcaSelect <- FindVariableFeatures(mcaSelect, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(mcaSelect)
mcaSelect <- ScaleData(mcaSelect, features = all.genes)
mcaSelect <- RunPCA(mcaSelect, features = VariableFeatures(object = mcaSelect))
mcaSelect <- RunUMAP(mcaSelect, reduction="pca", dims=1:20)
mcaSelect <- FindNeighbors(mcaSelect, reduction="pca", dims=1:20)
mcaSelect <- FindClusters(mcaSelect, resolution=0.8)

save(mcaSelect, file = "mca_processed_select.rda")
write.csv(GetAssayData(mcaSelect, slot = "scale.data"), file = "mca_processed_select.csv", row.names = T)
write.csv(test_cells$celltype, file = "labels_mca_select.csv", row.names = F)


# ImmGen reference

ref <- celldex::ImmGenData()
refSelect <- ref[,ref$label.main == "B cells" |
                   ref$label.main == "DC" |
                   ref$label.main == "Macrophages" |
                   ref$label.main == "Mast cells" |
                   ref$label.main == "Neutrophils" |
                   ref$label.main == "T cells"]

save(refSelect, file = "refImmScaledSelect.rda")
write.csv(ScaleData(assay(refSelect)), file = "refImmScaledSelect.csv", row.names = T)
write.csv(refSelect$label.main, file = "labels_refImm_select.csv", row.names = F)

#####


##### DATA PROCESSING with Mutual Nearest Neighbour (MNN) algorithm ####

# Get assay data
pbmc.data <- GetAssayData(pbmcAll, slot = "counts")
mca.data <- GetAssayData(mcaAll, slot = "counts")
row.names(mca.data) <- toupper(row.names(mca.data))

# Obtain the same features (genes) between sets
mca.data <- mca.data[(intersect(row.names(mca.data), row.names(pbmc.data))),]
pbmc.data <- pbmc.data[(intersect(row.names(mca.data), row.names(pbmc.data))),]

# Transpose data (cells x features)
mca.data <- t(mca.data)
pbmc.data <- t(pbmc.data)

# Sort by features
mca.data <- mca.data[,sort(colnames(mca.data))]
pbmc.data <- pbmc.data[,sort(colnames(pbmc.data))]

# Obtain a SingleCellExperiment object
sce1 <- as.data.frame(mca.data)
sce1 <- SingleCellExperiment(assays =  list(counts = as.matrix(sce1)))

sce2 <- as.data.frame(pbmc.data)
sce2 <- SingleCellExperiment(assays =  list(counts = as.matrix(sce2)))

# Compute log-normalized expression values using library size-derived size factors
# (with batchelor)
out <- multiBatchNorm(sce1, sce2)
sce1 <- out[[1]]
sce2 <- out[[2]]

# obtain top 5000 genes with the largest biological components of their variance
# (with scran)
dec1 <- modelGeneVar(sce1)
dec2 <- modelGeneVar(sce2)
combined.dec <- combineVar(dec1, dec2)
chosen.hvgs <- getTopHVGs(combined.dec, n=5000)


# MNN identification and correction
combined <- correctExperiments(A=sce1, B=sce2, PARAM=FastMnnParam())

# Visualitzation
combined <- runPCA(combined, subset_row=chosen.hvgs)
combined <- runTSNE(combined, dimred="corrected")
plotTSNE(combined, colour_by="batch")

