
# Processing of JArribas dataset
#################################

# Packages ####
library(Seurat)
library(SingleR)
library(dplyr)
#####

# Obtaining directories ####
baseDir <- file.path("C:/Users/Usuario/Desktop/master/tfm/scripts_SingleR_merce")
dataDir <- file.path(baseDir, "rawData")
analysisDir <- file.path(baseDir, "analysis")
resultsDir <- file.path(baseDir, "results")
######

## LOAD DATA FORM H5 FILES ####
dw4Data <- Read10X_h5(file=file.path(dataDir, "A0007_4w/filtered_feature_bc_matrix.h5"))
dw9Data <- Read10X_h5(file=file.path(dataDir, "A0008_9w/filtered_feature_bc_matrix.h5"))
di9Data <- Read10X_h5(file=file.path(dataDir, "A0009_9w_ind/filtered_feature_bc_matrix.h5"))
#####


# Processing ####
# WEEK 4
# creates Seurat Object
dw4 <- CreateSeuratObject(counts=dw4Data, project = 'dw4') #, min.cells=3, min.features=200)
dw4[["percent.mt"]] <- PercentageFeatureSet(dw4, pattern = "^mt-")


jpeg(file=file.path(resultsDir, paste("QC", "dw4", "VlnPlot", "jpeg", sep=".")), width=20, height=15, units="in", res=300)
VlnPlot(dw4, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))
dev.off()

jpeg(file=file.path(resultsDir, paste("QC", "dw4", "FeatureScatter", "jpeg", sep=".")), width=20, height=15, units="in", res=300)
plot1 <- FeatureScatter(dw4, feature1="nCount_RNA", feature2="percent.mt")
plot2 <- FeatureScatter(dw4, feature1="nCount_RNA", feature2="nFeature_RNA")
CombinePlots(plots=list(plot1,plot2))
dev.off()

# subseting the counts of Feature/RNA
dw4 <- subset(dw4, subset = nFeature_RNA > 200 & nFeature_RNA < 7500)

# Normalize data
dw4 <- NormalizeData(dw4)

# Find the most variable features (method = vst??)
dw4 <- FindVariableFeatures(dw4, selection.method="vst", nfeatures=2000)
top20.dw4 <- head(VariableFeatures(dw4), 20)


jpeg(file=file.path(resultsDir, paste("QC", "dw4","VariableFeatures", "Top10", "jpeg", sep=".")), width=20, height=15, units="in", res=300)
plot1 <- VariableFeaturePlot(dw4)
plot2 <- LabelPoints(plot=plot1, points=top20.dw4, repel=TRUE)
CombinePlots(plots=list(plot1,plot2))
dev.off()

saveRDS(dw4, file=file.path(resultsDir, paste("dw4", "Log", "rds", sep=".")))

# WEEK 9
dw9 <- CreateSeuratObject(counts=dw9Data, project = 'dw9') #, min.cells=3, min.features=200)
dw9[["percent.mt"]] <- PercentageFeatureSet(dw9, pattern = "^mt-")
jpeg(file=file.path(resultsDir, paste("QC", "dw9", "VlnPlot", "jpeg", sep=".")), width=20, height=15, units="in", res=300)
VlnPlot(dw9, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))
dev.off()

jpeg(file=file.path(resultsDir, paste("QC", "dw9", "FeatureScatter", "jpeg", sep=".")), width=20, height=15, units="in", res=300)
plot1 <- FeatureScatter(dw9, feature1="nCount_RNA", feature2="percent.mt")
plot2 <- FeatureScatter(dw9, feature1="nCount_RNA", feature2="nFeature_RNA")
CombinePlots(plots=list(plot1,plot2))
dev.off()

dw9 <- subset(dw9, subset = nFeature_RNA > 200 & nFeature_RNA < 7500)

dw9 <- NormalizeData(dw9)
dw9 <- FindVariableFeatures(dw9, selection.method="vst", nfeatures=2000)
top20.dw9 <- head(VariableFeatures(dw9), 20)

jpeg(file=file.path(resultsDir, paste("QC", "dw9","VariableFeatures", "Top10", "jpeg", sep=".")), width=20, height=15, units="in", res=300)
plot1 <- VariableFeaturePlot(dw9)
plot2 <- LabelPoints(plot=plot1, points=top20.dw9, repel=TRUE)
CombinePlots(plots=list(plot1,plot2))
dev.off()

saveRDS(dw9, file=file.path(resultsDir, paste("dw9", "Log", "rds", sep=".")))

# WEEK 9 + TREATMENT
di9 <- CreateSeuratObject(counts=di9Data, project = 'di9') #, min.cells=3, min.features=200)
di9[["percent.mt"]] <- PercentageFeatureSet(di9, pattern = "^mt-")
jpeg(file=file.path(resultsDir, paste("QC", "di9", "VlnPlot", "jpeg", sep=".")), width=20, height=15, units="in", res=300)
VlnPlot(di9, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))
dev.off()

jpeg(file=file.path(resultsDir, paste("QC", "di9", "FeatureScatter", "jpeg", sep=".")), width=20, height=15, units="in", res=300)
plot1 <- FeatureScatter(di9, feature1="nCount_RNA", feature2="percent.mt")
plot2 <- FeatureScatter(di9, feature1="nCount_RNA", feature2="nFeature_RNA")
CombinePlots(plots=list(plot1,plot2))
dev.off()

di9 <- subset(di9, subset = nFeature_RNA > 200 & nFeature_RNA < 7500)

di9 <- NormalizeData(di9)
di9 <- FindVariableFeatures(di9, selection.method="vst", nfeatures=2000)
top20.di9 <- head(VariableFeatures(di9), 20)

jpeg(file=file.path(resultsDir, paste("QC", "di9","VariableFeatures", "Top10", "jpeg", sep=".")), width=20, height=15, units="in", res=300)
plot1 <- VariableFeaturePlot(di9)
plot2 <- LabelPoints(plot=plot1, points=top20.di9, repel=TRUE)
CombinePlots(plots=list(plot1,plot2))
dev.off()


saveRDS(di9, file=file.path(resultsDir, paste("di9", "Log", "rds", sep=".")))


## Dimentionality reduction of individual samples

# genes
allGenes.dw4 <- rownames(dw4)
allGenes.dw9 <- rownames(dw9)
allGenes.di9 <- rownames(di9)

# scale
dw4 <- ScaleData(dw4, features=allGenes.dw4)
dw9 <- ScaleData(dw9, features=allGenes.dw9)
di9 <- ScaleData(di9, features=allGenes.di9)

# run PCA
dw4 <- RunPCA(dw4, features=VariableFeatures(object=dw4))
dw9 <- RunPCA(dw9, features=VariableFeatures(object=dw9))
di9 <- RunPCA(di9, features=VariableFeatures(object=di9))


# Integration anchors
dwAnchors <- FindIntegrationAnchors(object.list=c(dw4, dw9, di9),
                                    dims=1:20)
dwIntegrated <- IntegrateData(anchorset=dwAnchors, dims=1:20)
dwIntegrated <- ScaleData(dwIntegrated)
dwIntegrated <- RunPCA(dwIntegrated, npcs=20)
dwIntegrated <- RunUMAP(dwIntegrated, reduction="pca", dims=1:20)
dwIntegrated <- FindNeighbors(dwIntegrated, reduction="pca", dims=1:20)
dwIntegrated <- FindClusters(dwIntegrated, resolution=0.8)
dwIntegrated <- RunUMAP(dwIntegrated, dims=1:20)

jpeg(file=file.path(resultsDir, paste("QC", "Integrated", "UMAP", "20D", "samples", "jpeg", sep=".")), width=8, height=6, units="in", res=300)
DimPlot(dwIntegrated, reduction="umap", group.by="orig.ident")
dev.off()

jpeg(file=file.path(resultsDir, paste("QC", "Integrated", "UMAP", "20D", "samples", "split", "jpeg", sep=".")), width=18, height=6, units="in", res=300)
DimPlot(dwIntegrated, reduction="umap", split.by="orig.ident")
dev.off()

jpeg(file=file.path(resultsDir, paste("QC", "Integrated", "UMAP", "20D", "clusters", "jpeg", sep=".")), width=8, height=6, units="in", res=300)
DimPlot(dwIntegrated, reduction="umap", label=TRUE)
dev.off()

jpeg(file=file.path(resultsDir, paste("QC", "Integrated", "20D", "Elbow", "jpeg", sep=".")), width=8, height=6, units="in", res=300)
ElbowPlot(dwIntegrated, ndims=20)
dev.off()

jpeg(file=file.path(resultsDir, paste("QC", "Integrated", "20D", "DimHeatmap", "jpeg", sep=".")), width=8, height=6, units="in", res=300)
DimHeatmap(dwIntegrated, dims=1:20, cells=500, balanced=TRUE)
dev.off()
#####


# Export rds object (for SingleR annoations) ####
saveRDS(dwIntegrated, file=file.path(resultsDir, "dwIntegrated.20D.rds"))



# Export JArribas and reference data as csv (for SVM annotations) ####
write.csv(GetAssayData(dwIntegrated, slot = "scale.data"), file = "JArribas_processed_select.csv", row.names = T)

DefaultAssay(dwIntegrated) <- "RNA"
write.csv(GetAssayData(dwIntegrated, slot = "scale.data"), file = "JArribas_processed_all.csv", row.names = T)


ref <- celldex::ImmGenData()
write.csv(ScaleData(assay(ref)), file = "refImmScaled.csv", row.names = T)
write.csv(ref$label.main, file = "labels_main_refImm.csv", row.names = F)
write.csv(ref$label.fine, file = "labels_specific_refImm.csv", row.names = F)
