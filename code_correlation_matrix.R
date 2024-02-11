
##########################################
#### Get cell type correlation matrix ####
##########################################


# Packages ####
library(corrplot)
library(ggplot2)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
####################

# change this line

# Get correlation matrix for PBMC #####

setwd("C:/Users/Usuario/Desktop/master/tfm/PBMC_data")
load("pbmc_processed_all.rda")

labels_pbmc_all <- read.csv("pbmc_labels.csv")

pbmc.data <- as.data.frame(pbmcAll@assays$RNA@counts)
cell_populations <- names(table(labels_pbmc_all))
mean_cell = NULL

for(i in 1:length(cell_populations)){
  mean_cell[[i]] <- apply(pbmc.data[,labels_pbmc_all == cell_populations[i]], 1, mean)
}

mean_cell_df <- matrix(as.vector(unlist(mean_cell)), nrow = 19128, ncol = 9)

corr_matrix <- cor(mean_cell_df)


colnames(corr_matrix) <- cell_populations
rownames(corr_matrix) <- cell_populations

mean(apply(corr_matrix, 1, function(x) sort(x, TRUE)[2])) # complexity

plot <- corrplot(corr_matrix,  method = 'circle', type = 'lower', insig='blank',
                 addCoef.col ='black', number.cex = 0.5, order = 'hclust', diag=TRUE,
                 tl.col = "black")


ggsave("correlation_plot.pdf",plot, width = 10,height=11, dpi=300, units = "in")

#####



# Get correlation matrix for MCA #####

setwd("C:/Users/Usuario/Desktop/master/tfm/HCL_data/")

load("mca_processed_all.rda")
labels_mca_all <- read.csv("labels_mca_all.csv")

mca.data <- as.data.frame(mcaAll@assays$RNA@counts)
cell_populations <- names(table(labels_mca_all))
cell_populations <- cell_populations[cell_populations != "Cd4+ T cell" &
                                       cell_populations != "Macrophage"]
mean_cell = NULL

for(i in 1:length(cell_populations)){
  mean_cell[[i]] <- apply(mca.data[,labels_mca_all == cell_populations[i]], 1, mean)
}

mean_cell_df <- matrix(as.vector(unlist(mean_cell)), nrow = 34947, ncol = 10)
corr_matrix <- cor(mean_cell_df)

mean(apply(corr_matrix, 1, function(x) sort(x, TRUE)[2]))


colnames(corr_matrix) <- cell_populations
rownames(corr_matrix) <- cell_populations

plot <- corrplot(corr_matrix,  method = 'circle', type = 'lower', insig='blank',
                 addCoef.col ='black', number.cex = 0.5, order = 'hclust', diag=TRUE,
                 tl.col = "black")


ggsave("correlation_plot.pdf",plot, width = 10,height=11, dpi=300, units = "in")

###############
