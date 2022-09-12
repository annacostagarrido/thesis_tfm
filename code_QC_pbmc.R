
##############################
#### Quality control PBMC ####
##############################


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



# Get Single cell violin plot for nFeature_RNA and nCount_RNA
pbmc <- PercentageFeatureSet(pbmc, "^MT-", col.name = "percent_mito")
pbmc <- PercentageFeatureSet(pbmc, "^RP[SL]", col.name = "percent_ribo")
# Percentage hemoglobin genes - includes all genes starting with HB except HBP.
pbmc <- PercentageFeatureSet(pbmc, "^HB[^(P)]", col.name = "percent_hb")

pbmc <- PercentageFeatureSet(pbmc, "PECAM1|PF4", col.name = "percent_plat")


feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")
VlnPlot(pbmc, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) +
  NoLegend()



# Get the scatter plot for these variables
FeatureScatter(pbmc, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)

