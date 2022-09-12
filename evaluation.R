
#######################################
#### EVALUATION of all experiments ####
#######################################

# Packages ####
library(readr)
library(readxl)
library(tidyverse)
library(crfsuite)
library(udpipe)
library(ggplot2)
library(scales)
library(ggpubr)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(pals)
library(grid)
library(gridExtra)
######


# Function for the evaluation

evaluation <- function(true, pred, cells, n_correct, n_cells) {
  # Get evaluation in terms of f1 score and global % correctly classified
  evaluation <- crf_evaluation(pred = pred, obs = true)

  # Remove those labels not seen in the test data
  evaluation <- evaluation$bylabel[evaluation$bylabel$support != 0,]

  # Change na f1 values to 0
  evaluation$f1[is.na(evaluation$f1)] <- 0

  # Compute the number of cells correctly classified and global percentage
  n <- sum(pred == true, na.rm = T)
  perc <- sum(pred == true, na.rm = T) / length(pred)

  # Get the global F1 score
  globalF1 <- sum(evaluation$f1 * evaluation$support, na.rm = T) / sum(evaluation$support)

  # Get the specific F1 score for each cell seen in the test data
  specificF1 <- sprintf("%.2f",evaluation$f1)
  names(specificF1) <- evaluation$label
  # Check the order of the cells
  specificF1 <- specificF1[match(cells,names(specificF1))]


  # get the Final evaluation
  f1_evaluation <- c("Num. cells correctly classified (%)" = paste0(n, " (",round(perc * 100, 2), "%)"),
                     "F1 - score weighted" = round(globalF1, 2),
                     "Cell types in test data and reference",
                     specificF1[1:n_correct],
                     "Cell types in test data but not in the reference",
                     specificF1[(n_correct +1) :n_cells])



  # Computation of the number of cells and % correctly classified for each cell

  # Global confusion matrix
  taula <- table(true, pred)

  # Confusion matrix with only those cells seen in both test and reference data
  taula_select <- taula[rownames(taula) %in% cells[1:n_correct],]
  taula_select <- taula_select[match(cells[1:n_correct], rownames(taula_select)),]

  # final n and % for each of the cells
  not_found <- rownames(taula_select)[!rownames(taula_select) %in% colnames(taula_select)]
  taula_select2 <- data.frame(as.data.frame.matrix(taula_select))
  colnames(taula_select2) <- colnames(taula_select)
  empty <- data.frame(matrix(NA, nrow = nrow(taula_select),
                             ncol = length(not_found)))

  taula_select2 <- data.frame(taula_select2, empty)
  colnames(taula_select2) <- c(colnames(taula_select), not_found)

  taula_select <- taula_select2

  # Confusion matrix with the percentages
  prop_taula <- (taula_select / table(true)[match(cells[1:n_correct],names(table(true)))])*100

  n <- unlist(sapply(rownames(taula_select), function(x) taula_select[x,x]))
  perc <- unlist(sapply(rownames(prop_taula), function(x) prop_taula[x,x]))

  n_perc <- paste0(n, " (", sprintf("%.2f", perc), "%)" )
  names(n_perc) <- rownames(taula_select)

  n_perc_evaluation <- n_perc[match(cells,names(n_perc))]
  names(n_perc_evaluation) <- cells

  return(list(F1 = f1_evaluation, n_perc  = n_perc_evaluation))
}

#####

# MCA - ImmGen classification ####
##################################

setwd("C:/Users/Usuario/Desktop/master/tfm/HCL_data")

true_all <- read.csv("labels_mca_all.csv")
true_select <- read.csv("labels_mca_select.csv")

true_all <- txt_recode(true_all$x, from = names(table(true_all)),
                       to = c("B cells", "B cell(Plasmocyte)", "T cells", "T cells",
                              "DC", "Erythroid cell", "Lymphocyte", "Macrophages",
                              "Mast cells", "Myoloid cell", "Neutrophils", "T cells"))

true_select <- txt_recode(true_select$x, from = names(table(true_select)),
                          to = c("B cells", "T cells", "T cells",
                                 "DC", "Macrophages",
                                 "Mast cells", "Neutrophils", "T cells"))

# 1) SVM with all cell types of test data and reference data
pred_all <- read.csv("./results/SVM_linear_ImmGen_scaled_All.csv")

SVM_all <- evaluation(true_all, pred_all$X0,
                      cells = c("T cells", "B cells", "DC",
                                "Neutrophils", "Macrophages",
                                "Mast cells", "B cell(Plasmocyte)",
                                "Erythroid cell", "Lymphocyte",
                                "Myoloid cell"),
                      n_correct = 6, n_cells = 10)



# 2) SVM with only the cell type that are in both datasets (test and reference)
pred_select <- read.csv("./results/SVM_linear_ImmGen_scaled_Select.csv")


SVM_select <- evaluation(true_select, pred_select$X0,
                         cells = c("T cells", "B cells", "DC",
                                   "Neutrophils", "Macrophages",
                                   "Mast cells"),
                         n_correct = 6, n_cells = 10)


# 3) SingleR (clusters) with all cell types of test data and reference data
pred_clusters_all <- read.csv("./results/SingleR_clusters_ImmGen_scaled_All.csv")

SingleR_clusters_all <- evaluation(true_all, pred_clusters_all$x,
                                   cells = c("T cells", "B cells", "DC",
                                             "Neutrophils", "Macrophages",
                                             "Mast cells", "B cell(Plasmocyte)",
                                             "Erythroid cell", "Lymphocyte",
                                             "Myoloid cell"),
                                   n_correct = 6, n_cells = 10)


# 4) SingleR (per each cell) with all cell types of test data and reference data
pred_cells_all <- read.csv("./results/SingleR_eachCell_ImmGen_scaled_All.csv")

SingleR_cells_all <- evaluation(true_all, pred_cells_all$x,
                                cells = c("T cells", "B cells", "DC",
                                          "Neutrophils", "Macrophages",
                                          "Mast cells", "B cell(Plasmocyte)",
                                          "Erythroid cell", "Lymphocyte",
                                          "Myoloid cell"),
                                n_correct = 6, n_cells = 10)


# 5) SingleR (clusters) with only the cell type that are in both datasets (test and reference)
pred_clusters_select <- read.csv("./results/SingleR_clusters_ImmGen_scaled_Select.csv")

SingleR_clusters_select <- evaluation(true_select, pred_clusters_select$x,
                                      cells = c("T cells", "B cells", "DC",
                                                "Neutrophils", "Macrophages",
                                                "Mast cells"),
                                      n_correct = 6, n_cells = 10)


# 6) SingleR (per each cell) with only the cell type that are in both datasets (test and reference)
pred_cells_select <- read.csv("./results/SingleR_eachCell_ImmGen_scaled_Select.csv")

SingleR_cells_select <- evaluation(true_select, pred_cells_select$x,
                                   cells = c("T cells", "B cells", "DC",
                                             "Neutrophils", "Macrophages",
                                             "Mast cells"),
                                   n_correct = 6, n_cells = 10)


# 7) MCA prediction using scPred

setwd("C:/Users/Usuario/Desktop/master/tfm/HCL_data")

true_all <- read.csv("labels_mca_all.csv")

setwd("C:/Users/Usuario/Desktop/master/tfm/sctype/")


pred_sctype <- read.csv("mca/pred_true.csv")

pred_sctype <- txt_recode(pred_sctype$pred, from = names(table(pred_sctype$pred)),
                          to = c("Basophils", "Erythroid cell", "Intermediate monocytes",
                                 "Macrophage", "B cell", "Cd8+ T cell", "Natural killer  cells",
                                 "Neutrophil", "B cell", "B cell"))


sctype_mca <- evaluation(true_all$x, pred_sctype, names(table(true_all$x)), n_cells = length(names(table(true_all$x))), length(names(table(true_all$x))))

openxlsx::write.xlsx(data.frame(cbind(sctype_mca$F1[1: 15], c(NA,NA, NA, sctype_mca$n_perc))),
                     file = "prediction_sctype_mca.xlsx", rowNames = TRUE)


F1_perc_global1 <- cbind(SVM_all$F1[1:2],
                         SingleR_clusters_all$F1[1:2],
                         SingleR_cells_all$F1[1:2])

F1_perc_global2 <- cbind(SVM_select$F1[1:2],
                         SingleR_clusters_select$F1[1:2],
                         SingleR_cells_select$F1[1:2])


F1_perc_global <- cbind(F1_perc_global1, F1_perc_global2)
colnames(F1_perc_global) <- c("SVM with linear kernel",		"SingleR - clusters",		"SingleR - each cell",
                              "SVM with linear kernel",		"SingleR - clusters",		"SingleR - each cell")


F1_perc_specific1 <- cbind(SVM_all$n_perc[1:6],
                           SVM_all$F1[4:9],
                           SingleR_clusters_all$n_perc[1:6],
                           SingleR_clusters_all$F1[4:9],
                           SingleR_cells_all$n_perc[1:6],
                           SingleR_cells_all$F1[4:9])

F1_perc_specific2 <- cbind(
  SVM_select$n_perc[1:6],
  SVM_select$F1[4:9],
  SingleR_clusters_select$n_perc[1:6],
  SingleR_clusters_select$F1[4:9],
  SingleR_cells_select$n_perc[1:6],
  SingleR_cells_select$F1[4:9])

F1_perc_specific <- cbind(F1_perc_specific1, F1_perc_specific2)
colnames(F1_perc_specific) <- c("SVM with linear kernel", NA,		"SingleR - clusters",	NA,	"SingleR - each cell", NA,
                                "SVM with linear kernel",	NA,	"SingleR - clusters",	NA,	"SingleR - each cell", NA)

openxlsx::write.xlsx(as.data.frame(F1_perc_global), file = "./results/table_F1_perc_global.xlsx", rowNames = T)
openxlsx::write.xlsx(data.frame(F1_perc_specific), file = "./results/table_F1_perc_specific.xlsx", rowNames = T)

#####


# PBMC - ImmGen / Monaco classification ####
############################################

setwd("C:/Users/Usuario/Desktop/master/tfm/PBMC_data/")

true_all <- read.csv("labels_pbmc_all.csv") # n = 6003
true_select <- read.csv("labels_pbmc_select.csv") # n = 3790

true_all_Monaco <- txt_recode(true_all$x, from = names(table(true_all)),
                              to = c("B cells", "Monocytes", "Monocytes",
                                     "CD4+ T cells", "Cytotoxic", "Dendritic cells",
                                     "Megakaryocyte", "NK cells", "Plasmacytoid"))

true_all_ImmGen <- txt_recode(true_all$x, from = names(table(true_all)),
                              to = c("B cells", "Monocytes", "Monocytes",
                                     "T cells", "Cytotoxic", "DC",
                                     "Megakaryocyte", "NK cells", "Plasmacytoid"))

true_select_Monaco <- txt_recode(true_select$x, from = names(table(true_select)),
                                 to = c("B cells", "Monocytes", "Monocytes",
                                        "CD4+ T cells", "Dendritic cells",
                                        "NK cells"))

true_select_ImmGen <- txt_recode(true_select$x, from = names(table(true_select)),
                                 to = c("B cells", "Monocytes", "Monocytes",
                                        "T cells", "DC", "NK cells"))





# 1) SVM with all cell types of test data and reference data - ImmGen

pred_all_ImmGen <- read.csv("./results/SVM_linear_Imm_scaled_All.csv")


SVM_all_ImmGen <- evaluation(true_all_ImmGen, pred_all_ImmGen$X0, cells = c("B cells", "DC", "Monocytes", "NK cells","T cells", "Cytotoxic",
                                                                            "Megakaryocyte", "Plasmacytoid"),
                             n_correct = 5, n_cells = 8)




# 2) SVM with all cell types of test data and reference data - Monaco
pred_all_Monaco <- read.csv("./results/SVM_linear_Monaco_scaled_All.csv")

pred_all_Monaco$X0[pred_all_Monaco$X0 == "T cells"] <- "CD4+ T cells"
SVM_all_Monaco <- evaluation(true_all_Monaco, pred_all_Monaco$X0, cells = c("B cells", "Dendritic cells",  "Monocytes",
                                                                            "NK cells", "CD4+ T cells",
                                                                            "Cytotoxic", "Megakaryocyte", "Plasmacytoid"),
                             n_correct = 5, n_cells = 8)


# 3) SVM with specific cell types found in both test data and reference data ImmGen
pred_select_ImmGen <- read.csv("./results/SVM_linear_Imm_scaled_Select.csv")


SVM_select_ImmGen <- evaluation(true_select_ImmGen, pred_select_ImmGen$X0, cells = c("B cells", "Monocytes",
                                                                                     "T cells", "DC", "NK cells"),
                                n_correct = 5, n_cells = 8)

# 4) SVM specific cell types found in both test data and reference data Monaco

pred_select_Monaco <- read.csv("./results/SVM_linear_Monaco_scaled_Select.csv")


SVM_select_Monaco <- evaluation(true_select_Monaco, pred_select_Monaco$X0, cells =  c("B cells", "Dendritic cells",  "Monocytes",
                                                                                      "NK cells", "CD4+ T cells"),
                                n_correct = 5, n_cells = 8)

# 5) SingleR (clusters) PBMCs - Monaco with all cell types
pred_clusters_all <- read.csv("./results/SingleR_clusters_Monaco_scaled_All.csv")
pred_clusters_all$x[pred_clusters_all$x == "T cells"] <- "CD4+ T cells"

SingleR_cluster_all <- evaluation(true_all_Monaco, pred_clusters_all$x, cells =  c("B cells", "Dendritic cells",  "Monocytes",
                                                                                   "NK cells", "CD4+ T cells",
                                                                                   "Cytotoxic", "Megakaryocyte", "Plasmacytoid"),
                                  n_correct = 5, n_cells = 8)


# 6) SingleR (per each cell) PBMCs - Monaco with all cell types
pred_cells_all <- read.csv("./results/SingleR_eachCell_Monaco_All.csv")
pred_cells_all$x[pred_cells_all$x == "T cells"] <- "CD4+ T cells"

SingleR_cell_all <- evaluation(true_all_Monaco, pred_cells_all$x, cells =  c("B cells", "Dendritic cells",  "Monocytes",
                                                                             "NK cells", "CD4+ T cells",
                                                                             "Cytotoxic", "Megakaryocyte", "Plasmacytoid"),
                               n_correct = 5, n_cells = 8)


# 7) SingleR (clusters) PBMCs - Monaco with equal cell types
pred_clusters_select <- read.csv("./results/SingleR_clusters_Monaco_scaled_Select.csv")
pred_clusters_select$x[pred_clusters_select$x == "T cells"] <- "CD4+ T cells"


SingleR_cluster_select <- evaluation(true_select_Monaco, pred_clusters_select$x, cells =  c("B cells", "Dendritic cells",  "Monocytes",
                                                                                            "NK cells", "CD4+ T cells"),
                                     n_correct = 5, n_cells = 8)


# 8) SingleR (per each cell) PBMCs - Monaco with equal cell types
pred_cells_select <- read.csv("./results/SingleR_eachCell_Monaco_Select.csv")
pred_cells_select$x[pred_cells_select$x == "T cells"] <- "CD4+ T cells"

SingleR_cell_select <- evaluation(true_select_Monaco, pred_cells_select$x, cells =  c("B cells", "Dendritic cells",  "Monocytes",
                                                                                      "NK cells", "CD4+ T cells"),
                                  n_correct = 5, n_cells = 8)

# 7) PBMC prediction using scPred

setwd("C:/Users/Usuario/Desktop/master/tfm/PBMC_data/")

true_all <- read.csv("labels_pbmc_all.csv")
true_all$x[true_all$x == "CD14+"] <- "Monocytes"
true_all$x[true_all$x == "CD16+"] <- "Monocytes"

setwd("C:/Users/Usuario/Desktop/master/tfm/sctype/")

pred_sctype <- read.csv("pbmc/pred_true.csv")

pred_sctype <- txt_recode(pred_sctype$pred, from = names(table(pred_sctype$pred)),
                          to = c("Basophils", "Monocytes", "Effector CD8+ T cells",
                                 "CD4+", "Dendritic", "B", "CD4+", "Naive CD8+ T cells",
                                 "Natural", "Monocytes", "Dendritic", "Platelets",
                                 "Progenitor cells", "CD4+"))


sctype_pbmc <- evaluation(true_all$x, pred_sctype, names(table(true_all$x)), n_cells = length(names(table(true_all$x))), length(names(table(true_all$x))))

openxlsx::write.xlsx(data.frame(cbind(sctype_pbmc$F1[1: 11], c(NA,NA, NA, sctype_pbmc$n_perc))),
                     file = "prediction_sctype_pbmc.xlsx", rowNames = TRUE)

#####


# Table PBMC - Monaco ####
##########################

F1_perc_global1 <- cbind(SVM_all_Monaco$F1[1:2],
                         SingleR_cluster_all$F1[1:2],
                         SingleR_cell_all$F1[1:2])

F1_perc_global2 <- cbind(SVM_select_Monaco$F1[1:2],
                         SingleR_cluster_select$F1[1:2],
                         SingleR_cell_select$F1[1:2])


F1_perc_global <- cbind(F1_perc_global1, F1_perc_global2)
colnames(F1_perc_global) <- c("SVM with linear kernel",		"SingleR - clusters",		"SingleR - each cell",
                              "SVM with linear kernel",		"SingleR - clusters",		"SingleR - each cell")


F1_perc_specific1 <- cbind(SVM_all_Monaco$n_perc[1:5],
                           SVM_all_Monaco$F1[4:8],
                           SingleR_cluster_all$n_perc[1:5],
                           SingleR_cluster_all$F1[4:8],
                           SingleR_cell_all$n_perc[1:5],
                           SingleR_cell_all$F1[4:8])

F1_perc_specific2 <- cbind(
  SVM_select_Monaco$n_perc[1:5],
  SVM_select_Monaco$F1[4:8],
  SingleR_cluster_select$n_perc[1:5],
  SingleR_cluster_select$F1[4:8],
  SingleR_cell_select$n_perc[1:5],
  SingleR_cell_select$F1[4:8])

F1_perc_specific <- cbind(F1_perc_specific1, F1_perc_specific2)
colnames(F1_perc_specific) <- c("SVM with linear kernel", NA,		"SingleR - clusters",	NA,	"SingleR - each cell", NA,
                                "SVM with linear kernel",	NA,	"SingleR - clusters",	NA,	"SingleR - each cell", NA)

openxlsx::write.xlsx(as.data.frame(F1_perc_global), file = "./results/table_F1_perc_global_Monaco.xlsx", rowNames = T)
openxlsx::write.xlsx(data.frame(F1_perc_specific), file = "./results/table_F1_perc_specific_Monaco.xlsx", rowNames = T)

#####


# Table PBMC - ImmGen ####
##########################

F1_perc_global_ImmGen <- cbind(SVM_all_ImmGen$F1[1:2],
                               SVM_select_ImmGen$F1[1:2])

colnames(F1_perc_global_ImmGen) <- c("SVM with linear kernel",
                                     "SVM with linear kernel")

F1_perc_specific_ImmGen <- cbind(SVM_all_ImmGen$n_perc[1:5],
                                 SVM_all_ImmGen$F1[4:8],
                                 SVM_select_ImmGen$n_perc[1:5],
                                 SVM_select_ImmGen$F1[4:8])
colnames(F1_perc_specific_ImmGen) <- c("SVM with linear kernel", NA,
                                       "SVM with linear kernel", NA)


openxlsx::write.xlsx(data.frame(F1_perc_global_ImmGen), file = "./results/table_F1_perc_global_ImmGen.xlsx", rowNames = T)
openxlsx::write.xlsx(data.frame(F1_perc_specific_ImmGen), file = "./results/table_F1_perc_select_ImmGen.xlsx", rowNames = T)

#####



# EVALUATION SVM with rejection option #####
############################################


setwd("C:/Users/Usuario/Desktop/master/tfm/analysis_2sc")

true_mca <- read.csv("labels_mca.csv")

true_mca <- txt_recode(true_mca$x, from = names(table(true_mca)),
                       to = c("B", "Unknown", "CD4+",
                              "CD4+", "Dendritic", "Unknown",
                              "Unknown", "Unknown", "Unknown",
                              "Unknown", "Unknown", "CD4+"))

true_pbmc <- read.csv("labels_pbmc.csv")

true_pbmc <- txt_recode(true_pbmc$x, from = names(table(true_pbmc)),
                        to = c("B cell", "Unknown", "Unknown", "Cd8+ T cell",
                               "Unknown", "Dendritic cell", "Unknown",
                               "Unknown", "Unknown"))


# Raw counts - NOT MNN (not aligned) ####

# SVM without rej pbmc as reference, MCA as test

pred_SVM_PBMC_train <- read.csv("./results/SVM_PBMC_train_MCA_test_notMNN.csv")

SVM_PBMC_train <- evaluation(true_mca, pred_SVM_PBMC_train$X0,
                             cells = c("B","CD4+", "Dendritic"),
                             n_correct = 3, n_cells = 3)

# SVM with rej pbmc as reference, MCA as test
pred_SVM_rej_PBMC_train <- read.csv("./results/SVM_rej70_PBMC_train_MCA_test_notMNN.csv")

SVM_rej_PBMC_train <- evaluation(true_mca, pred_SVM_rej_PBMC_train$X0,
                                 cells = c("B","CD4+", "Dendritic",
                                           "Unknown"),
                                 n_correct = 4, n_cells = 4)

# SVM with rej MCA as reference, pbmc as test
pred_SVM_MCA_train <- read.csv("./results/SVM_MCA_train_PBMC_test_notMNN.csv")

SVM_MCA_train <- evaluation(true_pbmc, pred_SVM_MCA_train$X0,
                            cells = c("B cell","Cd8+ T cell", "Dendritic cell"),
                            n_correct = 3, n_cells = 3)

# SVM with rej MCA as reference, pbmc as test
pred_SVM_rej_MCA_train <- read.csv("./results/SVM_rej70_MCA_train_PBMC_test_notMNN.csv")

SVM_rej_MCA_train <- evaluation(true_pbmc, pred_SVM_rej_MCA_train$X0,
                                cells = c("B cell","Cd8+ T cell", "Dendritic cell",
                                          "Unknown"),
                                n_correct = 4, n_cells = 4)
#####

# Seurat processing ####

# SVM without rej pbmc as reference, MCA as test

pred_SVM_PBMC_train_seurat <- read.csv("./results/SVM_PBMC_train_MCA_test_seurat.csv")

SVM_PBMC_train_seurat <- evaluation(true_mca, pred_SVM_PBMC_train_seurat$X0,
                                    cells = c("B","CD4+", "Dendritic"),
                                    n_correct = 3, n_cells = 3)

# SVM with rej pbmc as reference, MCA as test
pred_SVM_rej_PBMC_train_seurat <- read.csv("./results/SVM_rej70_PBMC_train_MCA_test_seurat.csv")

SVM_rej_PBMC_train_seurat <- evaluation(true_mca, pred_SVM_rej_PBMC_train_seurat$X0,
                                        cells = c("B","CD4+", "Dendritic",
                                                  "Unknown"),
                                        n_correct = 4, n_cells = 4)

# SVM with rej MCA as reference, pbmc as test
pred_SVM_MCA_train_seurat <- read.csv("./results/SVM_MCA_train_PBMC_test_seurat.csv")

SVM_MCA_train_seurat <- evaluation(true_pbmc, pred_SVM_MCA_train_seurat$X0,
                                   cells = c("B cell","Cd8+ T cell", "Dendritic cell"),
                                   n_correct = 3, n_cells = 3)

# SVM with rej MCA as reference, pbmc as test
pred_SVM_rej_MCA_train_seurat <- read.csv("./results/SVM_rej70_MCA_train_PBMC_test_seurat.csv")

SVM_rej_MCA_train_seurat <- evaluation(true_pbmc, pred_SVM_rej_MCA_train_seurat$X0,
                                       cells = c("B cell","Cd8+ T cell", "Dendritic cell",
                                                 "Unknown"),
                                       n_correct = 4, n_cells = 4)

#####

# MNN (aligned) ####

# SVM without rej pbmc as reference, MCA as test
pred_SVM_PBMC_train_MNN <- read.csv("./results/SVM_PBMC_train_MCA_test_MNN.csv")

SVM_PBMC_train_MNN <- evaluation(true_mca, pred_SVM_PBMC_train_MNN$X0,
                                 cells = c("B","CD4+", "Dendritic"),
                                 n_correct = 3, n_cells = 3)


# SVM with rej pbmc as reference, MCA as test
pred_SVM_rej_PBMC_train_MNN <- read.csv("./results/SVM_rej70_PBMC_train_MCA_test_MNN.csv")

SVM_rej_PBMC_train_MNN <- evaluation(true_mca, pred_SVM_rej_PBMC_train_MNN$X0,
                                     cells = c("B","CD4+", "Dendritic",
                                               "Unknown"),
                                     n_correct = 4, n_cells = 4)


# SVM without rej MCA as reference, pbmc as test

pred_SVM_MCA_train_MNN <- read.csv("./results/SVM_MCA_train_PBMC_test_MNN.csv")

SVM_MCA_train_MNN <- evaluation(true_pbmc, pred_SVM_MCA_train_MNN$X0,
                                cells = c("B cell","Cd8+ T cell", "Dendritic cell"),
                                n_correct = 3, n_cells = 3)


# SVM with rej MCA as reference, pbmc as test

pred_SVM_rej_MCA_train_MNN <- read.csv("./results/SVM_rej70_MCA_train_PBMC_test_MNN.csv")

SVM_rej_MCA_train_MNN <- evaluation(true_pbmc, pred_SVM_rej_MCA_train_MNN$X0,
                                    cells = c("B cell","Cd8+ T cell", "Dendritic cell",
                                              "Unknown"),
                                    n_correct = 4, n_cells = 4)

#####


# Table PBMC - train, MCA -test ####



F1_perc_global1 <- cbind(c(SVM_PBMC_train$F1[1:2],NA),
                         c(SVM_rej_PBMC_train$F1[1:2],
                           paste0(sum(pred_SVM_rej_PBMC_train == "Unknown"),
                                  " (",
                                  round(sum(pred_SVM_rej_PBMC_train$X0 == "Unknown") /
                                          length(pred_SVM_rej_PBMC_train$X0) * 100,2), "%)")))

F1_perc_global2 <- cbind(c(SVM_PBMC_train_seurat$F1[1:2],NA),
                         c(SVM_rej_PBMC_train_seurat$F1[1:2],
                           paste0(sum(pred_SVM_rej_PBMC_train_seurat == "Unknown"),
                                  " (",
                                  round(sum(pred_SVM_rej_PBMC_train_seurat$X0 == "Unknown") /
                                          length(pred_SVM_rej_PBMC_train_seurat$X0) * 100,2), "%)")))



F1_perc_global3 <- cbind(c(SVM_PBMC_train_MNN$F1[1:2],NA),
                         c(SVM_rej_PBMC_train_MNN$F1[1:2],
                           paste0(sum(pred_SVM_rej_PBMC_train_MNN == "Unknown"),
                                  " (",
                                  round(sum(pred_SVM_rej_PBMC_train_MNN$X0 == "Unknown") /
                                          length(pred_SVM_rej_PBMC_train_MNN$X0) * 100,2), "%)")))


F1_perc_global <- cbind(F1_perc_global1, F1_perc_global2, F1_perc_global3)
colnames(F1_perc_global) <- c("SVM with linear kernel",
                              "SVM with linear kernel and rejection option",
                              "SVM with linear kernel",
                              "SVM with linear kernel and rejection option",
                              "SVM with linear kernel",
                              "SVM with linear kernel and rejection option")


F1_perc_specific1 <- cbind(SVM_PBMC_train$n_perc[1:4],
                           c(SVM_PBMC_train$F1[4:6], NA),
                           SVM_rej_PBMC_train$n_perc[1:4],
                           SVM_rej_PBMC_train$F1[4:7])

F1_perc_specific2 <- cbind(SVM_PBMC_train_seurat$n_perc[1:4],
                           c(SVM_PBMC_train_seurat$F1[4:6], NA),
                           SVM_rej_PBMC_train_seurat$n_perc[1:4],
                           SVM_rej_PBMC_train_seurat$F1[4:7])


F1_perc_specific3 <- cbind(SVM_PBMC_train_MNN$n_perc[1:4],
                           c(SVM_PBMC_train_MNN$F1[4:6], NA),
                           SVM_rej_PBMC_train_MNN$n_perc[1:4],
                           SVM_rej_PBMC_train_MNN$F1[4:7])


F1_perc_specific <- cbind(F1_perc_specific1, F1_perc_specific2, F1_perc_specific3)
colnames(F1_perc_specific) <- c("SVM with linear kernel", NA,		"SVM with linear kernel and rejection option",	NA,
                                "SVM with linear kernel", NA,		"SVM with linear kernel and rejection option",	NA,
                                "SVM with linear kernel", NA,		"SVM with linear kernel and rejection option",	NA)

openxlsx::write.xlsx(as.data.frame(F1_perc_global), file = "./results/table_F1_perc_PBMC_train.xlsx", rowNames = T)
openxlsx::write.xlsx(data.frame(F1_perc_specific), file = "./results/table_F1_perc_specific_PBMC_train.xlsx", rowNames = T)

#####


# Table MCA - train, PBMC -test ####


F1_perc_global1 <- cbind(c(SVM_MCA_train$F1[1:2],NA),
                         c(SVM_rej_MCA_train$F1[1:2],
                           paste0(sum(pred_SVM_rej_MCA_train == "Unknown"),
                                  " (",
                                  round(sum(pred_SVM_rej_MCA_train$X0 == "Unknown") /
                                          length(pred_SVM_rej_MCA_train$X0) * 100,2), "%)")))

F1_perc_global2 <- cbind(c(SVM_MCA_train_seurat$F1[1:2],NA),
                         c(SVM_rej_MCA_train_seurat$F1[1:2],
                           paste0(sum(pred_SVM_rej_MCA_train_seurat == "Unknown"),
                                  " (",
                                  round(sum(pred_SVM_rej_MCA_train_seurat$X0 == "Unknown") /
                                          length(pred_SVM_rej_MCA_train_seurat$X0) * 100,2), "%)")))


F1_perc_global3 <- cbind(c(SVM_MCA_train_MNN$F1[1:2],NA),
                         c(SVM_rej_MCA_train_MNN$F1[1:2],
                           paste0(sum(pred_SVM_rej_MCA_train_MNN == "Unknown"),
                                  " (",
                                  round(sum(pred_SVM_rej_MCA_train_MNN$X0 == "Unknown") /
                                          length(pred_SVM_rej_MCA_train_MNN$X0) * 100,2), "%)")))


F1_perc_global <- cbind(F1_perc_global1, F1_perc_global2, F1_perc_global3)
colnames(F1_perc_global) <- c("SVM with linear kernel",
                              "SVM with linear kernel and rejection option",
                              "SVM with linear kernel",
                              "SVM with linear kernel and rejection option",
                              "SVM with linear kernel",
                              "SVM with linear kernel and rejection option")


F1_perc_specific1 <- cbind(SVM_MCA_train$n_perc[1:4],
                           c(SVM_MCA_train$F1[4:6], NA),
                           SVM_rej_MCA_train$n_perc[1:4],
                           SVM_rej_MCA_train$F1[4:7])

F1_perc_specific2 <- cbind(SVM_MCA_train_seurat$n_perc[1:4],
                           c(SVM_MCA_train_seurat$F1[4:6], NA),
                           SVM_rej_MCA_train_seurat$n_perc[1:4],
                           SVM_rej_MCA_train_seurat$F1[4:7])


F1_perc_specific3 <- cbind(SVM_MCA_train_MNN$n_perc[1:4],
                           c(SVM_MCA_train_MNN$F1[4:6], NA),
                           SVM_rej_MCA_train_MNN$n_perc[1:4],
                           SVM_rej_MCA_train_MNN$F1[4:7])


F1_perc_specific <- cbind(F1_perc_specific1, F1_perc_specific2, F1_perc_specific3)
colnames(F1_perc_specific) <- c("SVM with linear kernel", NA,		"SVM with linear kernel and rejection option",	NA,
                                "SVM with linear kernel", NA,		"SVM with linear kernel and rejection option",	NA,
                                "SVM with linear kernel", NA,		"SVM with linear kernel and rejection option",	NA)

openxlsx::write.xlsx(as.data.frame(F1_perc_global), file = "./results/table_F1_perc_MCA_train.xlsx", rowNames = T)
openxlsx::write.xlsx(data.frame(F1_perc_specific), file = "./results/table_F1_perc_specific_MCA_train.xlsx", rowNames = T)

#####


# General figure: UMAP plot and misclassification ####
######################################################

# Function to get misclassifications
data_misclassification <- function(true, pred, cells, all_cells){
  taula <- table(true, pred)

  taula_select <- taula[rownames(taula) %in% cells, ]
  n <- apply(taula_select, 1, sum)
  taula_select/n * 100

  perc_cell_type_not_equal <- list()
  for(i in 1:nrow(taula_select)){
    cell_type <- taula_select[i,]
    perc_cell_type_not_equal[[i]] <- paste0(round(cell_type[cell_type != 0]/n[i] *100, 2))
    names(perc_cell_type_not_equal[[i]]) <- names(cell_type[cell_type != 0]) }

  names(perc_cell_type_not_equal) <- rownames(taula_select)

  all <- NULL

  for(i in 1:length(cells)){
    all[[i]] <- sapply(perc_cell_type_not_equal[[i]], '[', seq(max(sapply(perc_cell_type_not_equal[[i]], length)))) %>% t() %>% as.data.frame()
    all[[i]] <- as.numeric(all[[i]])
    names(all[[i]]) <- names(perc_cell_type_not_equal[[i]])
  }

  for(i in 1:length(cells)){
    all_cells_not <- all_cells[!all_cells %in% names(all[[i]])]
    all[[i]] <- c(all[[i]], rep(0, length(all_cells_not)))
    names(all[[i]]) <- c(names(all[[i]])[1:(length(all[[i]]) - length(all_cells_not))], all_cells_not)
    all[[i]] <- all[[i]][sort(names(all[[i]]))]}

  misclassification <- matrix(unlist(all), nrow = length(cells), byrow = TRUE)
  rownames(misclassification) <- cells
  colnames(misclassification) <- names(all[[1]])

  return(misclassification)
}

#####

# Function to build a ggplot ####

get_ggplot<- function(true, pred, cells, all_cells, labels_x_title,
                      color = NULL){

  all <- NULL

  for (i in 1:length(cells)){
    all[[i]] <- data_misclassification(true, pred, cells, all_cells)[i,]
    all[[i]] <- data.frame(names(all[[i]]), all[[i]])
    colnames(all[[i]]) <- c(cells[i], "Percentage")
  }

  all <- lapply(all, function(dat) {
    dat$type <- colnames(dat)[1]
    colnames(dat)[1] <- "variable"
    dat
  })

  all <- lapply(all, function(dat){
    dat$variable <- factor(dat$variable, levels = all_cells)
    return(dat)
  })
  all <- do.call(rbind, all)


  if(labels_x_title == "only title"){

    return(ggplot(all, aes(variable, Percentage/100, fill = variable)) +
             geom_col() +
             theme_bw() +
             theme(axis.text.x=element_blank(),
                   axis.ticks.x=element_blank(),
                   axis.title.x = element_blank(),
                   legend.position="none",
                   plot.margin=unit(c(0.6,0.5,0,1), "cm")) +
             scale_y_continuous(labels=scales::percent, limits = c(0,1)) +
             facet_wrap(~ type, scales = "free_x", ncol = length(cells)) +
             scale_fill_manual(name = "Cell type",
                               values = color) +
             labs(y = "Percentage") )
  }

  if(labels_x_title == "no title and labels"){

    return(ggplot(all, aes(variable, Percentage/100, fill = variable)) +
             geom_col() +
             theme_bw() +
             theme(axis.text.x=element_blank(),
                   axis.ticks.x=element_blank(),
                   axis.title.x = element_blank(),
                   legend.position="none") +
             scale_y_continuous(labels=scales::percent, limits = c(0,1)) +
             facet_wrap(~ type, scales = "free_x", ncol = length(cells)) +
             theme(strip.text.x = element_blank(),
                   plot.margin=unit(c(0.6,0.5,0,1), "cm")) +
             scale_fill_manual(name = "Cell type",
                                 values = color) +
             xlab("Cell type") + ylab("Percentage") )
  }

  if(labels_x_title == "only labels"){

    return(ggplot(all, aes(variable, Percentage/100, fill = variable)) +
             geom_col() +
             theme_bw() +
             theme(axis.text.x=element_blank(),
                   axis.ticks.x=element_blank(),
                   axis.title.x = element_blank(),
                   legend.position="none") +
             scale_y_continuous(labels=scales::percent, limits = c(0,1)) +
             facet_wrap(~ type, scales = "free_x", ncol = length(cells)) +
             theme(strip.text.x = element_blank(),
                   plot.margin=unit(c(0.6,0.5,0,1), "cm")) +
             scale_fill_manual(name = "Cell type",
                               values = color) +
             xlab("Cell type") + ylab("Percentage") )


  }
}


#####

# General Figure PBMCs - Monaco #####

setwd("C:/Users/Usuario/Desktop/master/tfm/PBMC_data/")


true_all <- read.csv("labels_pbmc_all.csv")
true_all_Monaco <- txt_recode(true_all$x, from = names(table(true_all)),
                              to = c("B cells", "Monocytes", "Monocytes",
                                     "CD4+ T cells", "Cytotoxic", "Dendritic cells",
                                     "Megakaryocyte", "NK cells", "Plasmacytoid"))

pred_all_Monaco <- read.csv("./results/SVM_linear_Monaco_scaled_All.csv")
pred_clusters_all <- read.csv("./results/SingleR_clusters_Monaco_scaled_All.csv")
pred_cells_all <- read.csv("./results/SingleR_eachCell_Monaco_All.csv")
pred_sctype <- read.csv("C:/Users/Usuario/Desktop/master/tfm/sctype/pbmc/pred_true.csv")
pred_sctype <- txt_recode(pred_sctype$pred, from = names(table(pred_sctype$pred)),
                          to = c("Basophils", "Monocytes", "CD8+ T cells",
                                 "CD4+ T cells", "Dendritic cells", "B cells",
                                 "CD4+ T cells", "CD8+ T cells",
                                 "NK cells", "Monocytes", "Dendritic cells", "Platelets",
                                 "Progenitors", "CD4+ T cells"))



# Set order of cell types and color

only_reference <- sort(unique(
  c(pred_all_Monaco$X0[!pred_all_Monaco$X0 %in% names(table(true_all_Monaco))],
    pred_clusters_all$x[!pred_clusters_all$x %in% names(table(true_all_Monaco))],
    pred_cells_all$x[!pred_cells_all$x %in% names(table(true_all_Monaco))])
))

only_sctype <- sort(unique(
  pred_sctype[!pred_sctype %in% names(table(true_all_Monaco))]
))

only_sctype <- only_sctype[3]


order_color <- data.frame(
  cell_type = c("B cells",
                "CD4+ T cells",
                "Dendritic cells",
                "Monocytes",
                "NK cells",
                "Cytotoxic",
                "Megakaryocyte",
                "Plasmacytoid",
                only_reference,
                only_sctype),
  color = cols25(13))


# UMAP plot - All cell types
load("pbmc_processed_all.rda")
true_all <- factor(true_all_Monaco, levels = order_color$cell_type)
Idents(pbmcAll) <- true_all
true <- UMAPPlot(pbmcAll, reduction="umap", cols = order_color$color[order_color$cell_type %in% levels(true_all)]
) +
  labs(title = "True annotation") + xlab("UMAP dimension 1") +
  ylab("UMAP dimension 2") + theme(legend.position="none")


pred_all <- factor(pred_all_Monaco$X0, levels = order_color$cell_type[order_color$cell_type %in% unique(pred_all_Monaco$X0)])
Idents(pbmcAll) <- pred_all
svm <- DimPlot(pbmcAll, reduction="umap", cols = order_color$color[order_color$cell_type %in% unique(pred_all)]) +
  labs(title = "SVM") + xlab("UMAP dimension 1") +
  ylab("UMAP dimension 2") + theme(legend.position="none")


pred_clusters_all <- factor(pred_clusters_all$x,
                            levels = order_color$cell_type[order_color$cell_type %in% unique(pred_clusters_all$x)])
Idents(pbmcAll) <- pred_clusters_all
cluster <- DimPlot(pbmcAll, reduction="umap", cols = order_color$color[order_color$cell_type %in% unique(pred_clusters_all)]) +
  labs(title = "SingleR - clusters") + xlab("UMAP dimension 1") +
  ylab("UMAP dimension 2") + theme(legend.position="none")


pred_cells_all <- factor(pred_cells_all$x,
                         levels = order_color$cell_type[order_color$cell_type %in% unique(pred_cells_all$x)])
Idents(pbmcAll) <- pred_cells_all
cell <- DimPlot(pbmcAll, reduction="umap", cols = order_color$color[order_color$cell_type %in% unique(pred_cells_all)]) +
  labs(title = "SingleR - cells") + xlab("UMAP dimension 1") +
  ylab("UMAP dimension 2") + theme(legend.position="none")


pred_sctype <- factor(pred_sctype,
                      levels = order_color$cell_type[order_color$cell_type %in% unique(pred_sctype)])
Idents(pbmcAll) <- pred_sctype
sctype <- DimPlot(pbmcAll, reduction="umap", cols = order_color$color[order_color$cell_type %in% unique(pred_sctype)]) +
  labs(title = "scType") + xlab("UMAP dimension 1") +
  ylab("UMAP dimension 2") + theme(legend.position="none")



# Get legend

rand <- sample(order_color$cell_type, 6003, replace = TRUE)
rand <- factor(rand, levels = order_color$cell_type)
Idents(pbmcAll) <- rand

p_legend <- DimPlot(pbmcAll, reduction="umap", cols = order_color$color) +
  labs(title = "Legend") + xlab("UMAP dimension 1") +
  ylab("UMAP dimension 2")

legend <- cowplot::get_legend(p_legend)
grid.newpage()
legend <- grid.draw(legend)



all_figures <- ggarrange(true, legend, svm, cluster, cell, sctype,
                         ncol = 2, nrow = 3)
ggsave("figure_pbmc_monaco_all.png",all_figures, width = 10,height=11, dpi=300, units = "in")



# Misclassification

cells <- c("Cytotoxic", "Megakaryocyte", "Plasmacytoid")
all_cells <- c(order_color$cell_type[c(1:5, 9:13)])
color_cells <- c(order_color$color[c(1:5, 9:13)])

plot1 <- get_ggplot(true_all, pred_all, cells, all_cells,
                    labels_x_title = "only title", color = color_cells)
plot2 <- get_ggplot(true_all, pred_clusters_all, cells, all_cells,
                    labels_x_title = "no title and labels", color = color_cells)
plot3 <- get_ggplot(true_all, pred_cells_all, cells, all_cells,
                    labels_x_title = "only labels", color = color_cells)

plot <- ggarrange(
  plot1, plot2, plot3,
  font.label = list(size = 9),
  labels = c("SVM", "SingleR - clusters",
             "SingleR - cells"),
  hjust = 0,
  vjust = 1.2,
  common.legend = FALSE, nrow = 3,
  heights = c(0.3, 0.3, 0.3),
  align  = "v"
)

setwd("C:/Users/Usuario/Desktop/master/tfm/results/")
ggsave("misclassification_pbmc_monaco.png",plot, width = 10,height=11, dpi=300, units = "in")


# UMAP plot - Common cell types
setwd("C:/Users/Usuario/Desktop/master/tfm/PBMC_data/")
load("pbmc_processed_select.rda")

true_select <- read.csv("labels_pbmc_select.csv")
true_select_Monaco <- txt_recode(true_select$x, from = names(table(true_select)),
                                 to = c("B cells", "Monocytes", "Monocytes",
                                        "CD4+ T cells", "Dendritic cells",
                                        "NK cells"))

pred_select_Monaco <- read.csv("./results/SVM_linear_Monaco_scaled_Select.csv")

pred_clusters_select <- read.csv("./results/SingleR_clusters_Monaco_scaled_Select.csv")
pred_clusters_select$x[pred_clusters_select$x == "T cells"] <- "CD4+ T cells"

pred_cells_select <- read.csv("./results/SingleR_eachCell_Monaco_Select.csv")
pred_cells_select$x[pred_cells_select$x == "T cells"] <- "CD4+ T cells"


true_select <- factor(true_select_Monaco, levels = order_color$cell_type[1:5])
Idents(pbmcSelect) <- true_select
true <- UMAPPlot(pbmcSelect, reduction="umap", cols = order_color$color[1:5]
) +
  labs(title = "True annotation") + xlab("UMAP dimension 1") +
  ylab("UMAP dimension 2") + theme(legend.position="none")


pred_select <- factor(pred_select_Monaco$X0, levels = order_color$cell_type[1:5])
Idents(pbmcSelect) <- pred_select
svm <- DimPlot(pbmcSelect, reduction="umap", cols = order_color$color[1:5]) +
  labs(title = "SVM") + xlab("UMAP dimension 1") +
  ylab("UMAP dimension 2")+ theme(legend.position="none")


pred_clusters_select <- factor(pred_clusters_select$x,
                               levels = order_color$cell_type[1:5])
Idents(pbmcSelect) <- pred_clusters_select
cluster <- DimPlot(pbmcSelect, reduction="umap", cols = order_color$color[order_color$cell_type %in% unique(pred_clusters_all)]) +
  labs(title = "SingleR - clusters") + xlab("UMAP dimension 1") +
  ylab("UMAP dimension 2") + theme(legend.position="none")


pred_cells_select <- factor(pred_cells_select$x,
                            levels = order_color$cell_type[1:5])
Idents(pbmcSelect) <- pred_cells_select
cell <- DimPlot(pbmcSelect, reduction="umap", cols = order_color$color[1:5]) +
  labs(title = "SingleR - cells") + xlab("UMAP dimension 1") +
  ylab("UMAP dimension 2") + theme(legend.position="none")



all_figures <- ggarrange(true, svm, cluster, cell,
                         ncol = 4, nrow = 1)

ggsave("figure_pbmc_monaco_select.png",all_figures, width = 17,height=5, dpi=300, units = "in")


######

# General figure MCA - ImmGen #####


setwd("C:/Users/Usuario/Desktop/master/tfm/HCL_data")

true_all <- read.csv("labels_mca_all.csv")
true_all <- txt_recode(true_all$x, from = names(table(true_all)),
                       to = c("B cells", "B cell(Plasmocyte)", "T cells", "T cells",
                              "DC", "Erythroid cell", "Lymphocyte", "Macrophages",
                              "Mast cells", "Myoloid cell", "Neutrophils", "T cells"))

pred_all <- read.csv("./results/SVM_linear_ImmGen_scaled_All.csv")
pred_clusters_all <- read.csv("./results/SingleR_clusters_ImmGen_scaled_All.csv")
pred_clusters_all$x[pred_clusters_all$x == "Stem cell"] <- "Stem cells"
pred_cells_all <- read.csv("./results/SingleR_eachCell_ImmGen_scaled_All.csv")
pred_sctype <- read.csv("C:/Users/Usuario/Desktop/master/tfm/sctype/mca/pred_true.csv")
pred_sctype <- txt_recode(pred_sctype$pred, from = names(table(pred_sctype$pred)),
                          to = c("Basophils", "Erythroid cell", "Intermediate monocytes",
                                 "Macrophages", "B cells", "T cells", "NK cells",
                                 "Neutrophils", "B cells", "B cells"))


# Set order of cell types and color

only_reference <- sort(unique(
  c(pred_all$X0[!pred_all$X0 %in% names(table(true_all))],
    pred_clusters_all$x[!pred_clusters_all$x %in% names(table(true_all))],
    pred_cells_all$x[!pred_cells_all$x %in% names(table(true_all))])
))

only_sctype <- sort(unique(
  pred_sctype[!pred_sctype %in% names(table(true_all))]
))

only_sctype <- only_sctype[2]


order_color <- data.frame(
  cell_type = c("T cells",
                "B cells",
                "DC",
                "Neutrophils",
                "Macrophages",
                "Mast cells",
                "B cell(Plasmocyte)",
                "Erythroid cell",
                "Lymphocyte",
                "Myeloid cell",
                only_reference,
                only_sctype),
  color = cols25(23))


# UMAP - All cell types
load("mca_processed_all.rda")

true_all <- factor(true_all, levels = order_color$cell_type)
Idents(mcaAll) <- true_all
true <- UMAPPlot(mcaAll, reduction="umap", cols = order_color$color[order_color$cell_type %in% levels(true_all)]
) +
  labs(title = "True annotation") + xlab("UMAP dimension 1") +
  ylab("UMAP dimension 2") + theme(legend.position="none")


pred_all <- factor(pred_all$X0, levels = order_color$cell_type[order_color$cell_type %in% unique(pred_all$X0)])
Idents(mcaAll) <- pred_all
svm <- DimPlot(mcaAll, reduction="umap", cols = order_color$color[order_color$cell_type %in% unique(pred_all)]) +
  labs(title = "SVM") + xlab("UMAP dimension 1") +
  ylab("UMAP dimension 2") + theme(legend.position="none")


pred_clusters_all <- factor(pred_clusters_all$x,
                            levels = order_color$cell_type[order_color$cell_type %in% unique(pred_clusters_all$x)])
Idents(mcaAll) <- pred_clusters_all
cluster <- DimPlot(mcaAll, reduction="umap", cols = order_color$color[order_color$cell_type %in% unique(pred_clusters_all)]) +
  labs(title = "SingleR - clusters") + xlab("UMAP dimension 1") +
  ylab("UMAP dimension 2") + theme(legend.position="none")


pred_cells_all <- factor(pred_cells_all$x,
                         levels = order_color$cell_type[order_color$cell_type %in% unique(pred_cells_all$x)])
Idents(mcaAll) <- pred_cells_all
cell <- DimPlot(mcaAll, reduction="umap", cols = order_color$color[order_color$cell_type %in% unique(pred_cells_all)]) +
  labs(title = "SingleR - cells") + xlab("UMAP dimension 1") +
  ylab("UMAP dimension 2") + theme(legend.position="none")


pred_sctype <- factor(pred_sctype,
                      levels = order_color$cell_type[order_color$cell_type %in% unique(pred_sctype)])
Idents(mcaAll) <- pred_sctype
sctype <- DimPlot(mcaAll, reduction="umap", cols = order_color$color[order_color$cell_type %in% unique(pred_sctype)]) +
  labs(title = "scType") + xlab("UMAP dimension 1") +
  ylab("UMAP dimension 2") + theme(legend.position="none")



# Get legend

rand <- sample(order_color$cell_type, 7095, replace = TRUE)
rand <- factor(rand, levels = order_color$cell_type)
Idents(mcaAll) <- rand

p_legend <- DimPlot(mcaAll, reduction="umap", cols = order_color$color) +
  labs(title = "Legend") + xlab("UMAP dimension 1") +
  ylab("UMAP dimension 2")

legend <- cowplot::get_legend(p_legend)
grid.newpage()
legend <- grid.draw(legend)



all_figures <- ggarrange(true, legend, svm, cluster, cell, sctype,
                         ncol = 2, nrow = 3)

ggsave("figure_mca_immgen.png",all_figures, width = 10,height=11, dpi=300, units = "in")



# Misclassification
cells <- c("B cell(Plasmocyte)", "Erythroid cell", "Lymphocyte",
           "Myeloid cell")
all_cells <- c(order_color$cell_type[c(1:6, 11:22)])
color_cells <- c(order_color$color[c(1:6, 11:22)])

plot1 <- get_ggplot(true_all, pred_all, cells, all_cells, labels_x_title = "only title",
                    color = color_cells)
plot2 <- get_ggplot(true_all, pred_clusters_all, cells, all_cells,
                    labels_x_title = "no title and labels", color = color_cells)
plot3 <- get_ggplot(true_all, pred_cells_all, cells, all_cells,
                    labels_x_title = "only labels", color = color_cells)

plot <- ggarrange(
  plot1, plot2, plot3,
  font.label = list(size = 9),
  labels = c("SVM", "SingleR - clusters",
            "SingleR - cells"),
  hjust = 0,
  vjust = 1.2,
  common.legend = FALSE, nrow = 3,
  heights = c(0.3, 0.3, 0.3),
  align  = "v"
)

setwd("C:/Users/Usuario/Desktop/master/tfm/results/")
ggsave("misclassification_mca_immgen.png",plot, width = 10,height=11, dpi=300, units = "in")


# UMAP - Common cell types

load("mca_processed_select.rda")

true_select <- read.csv("labels_mca_select.csv")
true_select <- txt_recode(true_select$x, from = names(table(true_select)),
                          to = c("B cells", "T cells", "T cells",
                                 "DC", "Macrophages",
                                 "Mast cells", "Neutrophils", "T cells"))


pred_select <- read.csv("./results/SVM_linear_ImmGen_scaled_Select.csv")
pred_clusters_select <- read.csv("./results/SingleR_clusters_ImmGen_scaled_Select.csv")
pred_cells_select <- read.csv("./results/SingleR_eachCell_ImmGen_scaled_Select.csv")



true_select <- factor(true_select, levels = order_color$cell_type[1:6])
Idents(mcaSelect) <- true_select
true <- UMAPPlot(mcaSelect, reduction="umap", cols = order_color$color[1:6]
) +
  labs(title = "True annotation") + xlab("UMAP dimension 1") +
  ylab("UMAP dimension 2")  + theme(legend.position="none")


pred_select <- factor(pred_select$X0, levels = order_color$cell_type[1:6])
Idents(mcaSelect) <- pred_select
svm <- DimPlot(mcaSelect, reduction="umap", cols = order_color$color[1:6]) +
  labs(title = "SVM") + xlab("UMAP dimension 1") +
  ylab("UMAP dimension 2") + theme(legend.position="none")


pred_clusters_select <- factor(pred_clusters_select$x,
                               levels = order_color$cell_type[1:6])
Idents(mcaSelect) <- pred_clusters_select
cluster <- DimPlot(mcaSelect, reduction="umap", cols = order_color$color[order_color$cell_type %in% unique(pred_clusters_all)]) +
  labs(title = "SingleR - clusters") + xlab("UMAP dimension 1") +
  ylab("UMAP dimension 2") + theme(legend.position="none")


pred_cells_select <- factor(pred_cells_select$x,
                            levels = order_color$cell_type[1:6])
Idents(mcaSelect) <- pred_cells_select
cell <- DimPlot(mcaSelect, reduction="umap", cols = order_color$color[1:6]) +
  labs(title = "SingleR - cells") + xlab("UMAP dimension 1") +
  ylab("UMAP dimension 2") + theme(legend.position="none")



all_figures <- ggarrange(true, svm, cluster, cell,
                         ncol = 4, nrow = 1)

ggsave("figure_mca_immgen_select.png",all_figures, width = 17,height=5, dpi=300, units = "in")



######

# Figure misclassification MCA as train, PBMC as test ####
setwd("C:/Users/Usuario/Desktop/master/tfm/analysis_2sc")

true_pbmc <- read.csv("labels_pbmc.csv")

true_pbmc <- txt_recode(true_pbmc$x, from = names(table(true_pbmc)),
                        to = c("B cell", "CD14+", "CD16+", "Cd8+ T cell",
                               "Cytotoxic", "Dendritic cell", "Megakaryocyte",
                               "Natural killer", "Plasmacytoid"))

pred_SVM_MCA_train <- read.csv("./results/SVM_MCA_train_PBMC_test_notMNN.csv")
pred_SVM_MCA_train$X0[pred_SVM_MCA_train$X0 == "Myoloid cell"] <- "Myeloid cell"
pred_SVM_MCA_train_seurat <- read.csv("./results/SVM_MCA_train_PBMC_test_seurat.csv")
pred_SVM_MCA_train_seurat$X0[pred_SVM_MCA_train_seurat$X0 == "Myoloid cell"] <- "Myeloid cell"
pred_SVM_MCA_train_MNN <- read.csv("./results/SVM_MCA_train_PBMC_test_MNN.csv")
pred_SVM_MCA_train_MNN$X0[pred_SVM_MCA_train_MNN$X0 == "Myoloid cell"] <- "Myeloid cell"


cells <- c("CD14+", "CD16+", "Cytotoxic",
           "Megakaryocyte", "Natural killer", "Plasmacytoid")

all_cells <- unique(c(names(table(pred_SVM_MCA_train)),
                      names(table(pred_SVM_MCA_train_seurat)),
                      names(table(pred_SVM_MCA_train_MNN))))

plot1 <- get_ggplot(true_pbmc, pred_SVM_MCA_train$X0, cells,
                    all_cells, labels_x_title = "only title",
                    color = hue_pal()(9))
plot2 <- get_ggplot(true_pbmc, pred_SVM_MCA_train_seurat$X0, cells,
                    all_cells, labels_x_title = "no title and labels",
                    color = hue_pal()(9))
plot3 <- get_ggplot(true_pbmc, pred_SVM_MCA_train_MNN$X0, cells,
                    all_cells, labels_x_title = "only labels",
                    color = hue_pal()(9))


plot <- ggarrange(
  plot1, plot2, plot3,
  labels = c("A", "B", "C"),
  font.label = list(size = 12),
  common.legend = TRUE, legend = "right", nrow = 3,
  heights = c(0.34, 0.3, 0.48),
  align  = "v"
)

setwd("C:/Users/Usuario/Desktop/master/tfm/results/")
ggsave("misclassification_mcaTrain_pbmcTest.png",plot, width = 10,height=11, dpi=300, units = "in")


######

# Figure misclassification PBMC as train, MCA as test ####
setwd("C:/Users/Usuario/Desktop/master/tfm/analysis_2sc")

true_mca <- read.csv("labels_mca.csv")
true_mca <- txt_recode(true_mca$x, from = names(table(true_mca)),
                       to = c("B", "B cell(Plasmocyte)", "CD4+",
                              "CD4+", "Dendritic", "Erythroid cell",
                              "Lymphocyte", "Macrophage", "Mast cell",
                              "Myeloid cell", "Neutrophil", "CD4+"))

pred_SVM_PBMC_train <- read.csv("./results/SVM_PBMC_train_MCA_test_notMNN.csv")
pred_SVM_PBMC_train_seurat <- read.csv("./results/SVM_PBMC_train_MCA_test_seurat.csv")
pred_SVM_PBMC_train_MNN <- read.csv("./results/SVM_PBMC_train_MCA_test_MNN.csv")

cells <- c("B cell(Plasmocyte)", "Erythroid cell", "Lymphocyte",
           "Macrophage", "Mast cell", "Myeloid cell", "Neutrophil")
all_cells <- unique(c(names(table(pred_SVM_PBMC_train)),
                      names(table(pred_SVM_PBMC_train_MNN))))


plot1 <- get_ggplot(true_mca, pred_SVM_PBMC_train$X0, cells, all_cells, labels_x_title = "only title",
                    color = hue_pal()(9))
plot2 <- get_ggplot(true_mca, pred_SVM_PBMC_train_seurat$X0, cells,
                    all_cells, labels_x_title = "no title and labels",
                    color = hue_pal()(9))
plot3 <- get_ggplot(true_mca, pred_SVM_PBMC_train_MNN$X0, cells,
                    all_cells, labels_x_title = "only labels",
                    color = hue_pal()(9))

plot <- ggarrange(
  plot1, plot2, plot3,
  labels = c("A", "B", "C"),
  font.label = list(size = 12),
  common.legend = TRUE, legend = "right", nrow = 3,
  heights = c(0.34, 0.3, 0.48),
  align  = "v"
)

setwd("C:/Users/Usuario/Desktop/master/tfm/results/")
ggsave("misclassification_pbmcTrain_mcaTest.png",plot, width = 10,height=11, dpi=300, units = "in")



#####





#######################
# Computation time ####
#######################


# MCA - ImmGen ####


setwd("C:/Users/Usuario/Desktop/master/tfm/HCL_data/results/")

# All:

SVM_all <- read.csv("SVM_linear_Time_ImmGen_scaled_All.csv")
SingleR_cells_all <- read.table("SingleR_time_eachCell_ImmGen_scaled_All.csv", dec = ",", header = T)
SingleR_clusters_all <-read.table("SingleR_time_clusters_ImmGen_scaled_All.csv", dec = ",", header = T)
sctype <- read.table("time_sctype_mca.csv", dec = ",", header = T)

# Select:

SVM_select <- read.csv("SVM_linear_Time_ImmGen_scaled_Select.csv")
SingleR_cells_select <- read.table("SingleR_time_eachCell_ImmGen_scaled_Select.csv", dec = ",", header = T)
SingleR_clusters_select <- read.table("SingleR_time_clusters_ImmGen_scaled_Select.csv", dec = ",", header = T)


time_mca <- data.frame(time = c(SVM_all$Training.time, SingleR_clusters_all$x,
                                SingleR_cells_all$x, SVM_select$Training.time,
                                SingleR_clusters_select$x, SingleR_cells_select$x,
                                sctype$x),
                       setting = c(rep("All cell types", 3), rep("Common cell types", 3), "All cell types"),
                       method = c(rep(c("SVM",
                                        "SingleR - clusters", "SingleR - cells"), 2),
                                  "scType"))

time_mca$method <- factor(time_mca$method, levels = c("SVM",
                                                      "SingleR - clusters", "SingleR - cells", "scType"))

#####

# PBMC - Monaco ####

setwd("C:/Users/Usuario/Desktop/master/tfm/PBMC_data/results/")

# All:

SVM_all <- read.csv("SVM_linear_Time_Monaco_scaled_All.csv")
SingleR_cells_all <- read.table("SingleR_time_eachCell_Monaco_All.csv", dec = ",", header = T)
SingleR_clusters_all <-read.table("SingleR_time_clusters_Monaco_scaled_All.csv", dec = ",", header = T)
sctype <- read.table("time_sctype_pbmc.csv", dec = ",", header = T)


# Select:

SVM_select <- read.csv("SVM_linear_Time_Monaco_scaled_Select.csv")
SingleR_cells_select <- read.table("SingleR_time_eachCell_Monaco_Select.csv", dec = ",", header = T)
SingleR_clusters_select <- read.table("SingleR_time_clusters_Monaco_Select.csv", dec = ",", header = T)


time_pbmc <- data.frame(time = c(SVM_all$Training.time, SingleR_clusters_all$x,
                                 SingleR_cells_all$x, SVM_select$Training.time,
                                 SingleR_clusters_select$x, SingleR_cells_select$x,
                                 sctype$x),
                        setting = c(rep("All cell types", 3), rep("Common cell types", 3), "All cell types"),
                        method = c(rep(c("SVM",
                                         "SingleR - clusters", "SingleR - cells"), 2),
                                   "scType"))
time_pbmc$method <- factor(time_pbmc$method, levels = c("SVM",
                                                        "SingleR - clusters", "SingleR - cells", "scType"))

#####

# Plot with computation time for MCA-ImmGen and PBMC-Monaco ####

plot1 <- ggplot(time_mca, aes(method, time, fill = method)) +
  geom_col() +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x = element_blank(),
        plot.margin=unit(c(0.1,0.5,0,0), "cm")) +
  geom_text(aes(label= paste0( sprintf("%.2f", time), " s")),
            vjust=-1, color="black", size=2,
            position = position_dodge(.9)) +
  facet_wrap(~ setting, scales = "free_x", ncol = 4) +
  scale_fill_discrete(name = "Method") +
  xlab("Method") + ylab("Time (seconds)") +
  ylim(0, 180)

plot2 <- ggplot(time_pbmc, aes(method, time, fill = method)) +
  geom_col() +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x = element_blank()) +
  geom_text(aes(label= paste0( sprintf("%.2f", time), " s")),
            vjust=-1, color="black", size=2,
            position = position_dodge(.9)) +
  facet_wrap(~ setting, scales = "free_x", ncol = 4)  +
  scale_fill_discrete(name = "Method") +
  xlab("")  + ylab("Time (seconds)") +
  ylim(0, 180)


plot <- ggarrange(
  plot1, plot2,
  labels = c("A", "B"),
  font.label = list(size = 12),
  common.legend = TRUE, legend = "right", nrow = 2,
  heights = c(0.50, 0.5),
  align  = "v"
)

setwd("C:/Users/Usuario/Desktop/master/tfm/results/")
ggsave("computation_time.pdf",plot, width = 10,height=11, dpi=300, units = "in")

#####



# MCA ref ####

setwd("C:/Users/Usuario/Desktop/master/tfm/analysis_2sc/results/")

SVM_counts <- read.csv("SVM_Time_MCA_train_PBMC_test_notMNN.csv")
SVM_seurat <- read.csv("SVM_Time_MCA_train_PBMC_test_seurat.csv")
SVM_MNN <- read.csv("SVM_Time_MCA_train_PBMC_test_MNN.csv")
SVM_rej_counts <- read.csv("SVM_rej70_Time_MCA_train_PBMC_test_notMNN.csv")
SVM_rej_seurat <- read.csv("SVM_rej70_Time_MCA_train_PBMC_test_seurat.csv")
SVM_rej_MNN <- read.csv("SVM_rej70_Time_MCA_train_PBMC_test_MNN.csv")

time_mca <- data.frame(time = c(SVM_counts$Training.time,
                                SVM_seurat$Training.time,
                                SVM_MNN$Training.time,
                                SVM_rej_counts$Training.time,
                                SVM_rej_seurat$Training.time,
                                SVM_rej_MNN$Training.time),
                       processing = rep(c("Raw counts", "Seurat", "MNN"), 2),
                       method = c(rep("SVM", 3),
                                  rep("SVMrejection", 3))
)
time_mca$processing <- factor(time_mca$processing, levels = c("Raw counts", "Seurat", "MNN"))
time_mca$method <- factor(time_mca$method, levels = c("SVM", "SVMrejection"))

#####


# PBMC ref ####

SVM_counts <- read.csv("SVM_Time_PBMC_train_MCA_test_notMNN.csv")
SVM_seurat <- read.csv("SVM_Time_PBMC_train_MCA_test_seurat.csv")
SVM_MNN <- read.csv("SVM_Time_PBMC_train_MCA_test_MNN.csv")
SVM_rej_counts <- read.csv("SVM_rej70_Time_PBMC_train_MCA_test_notMNN.csv")
SVM_rej_seurat <- read.csv("SVM_rej70_Time_PBMC_train_MCA_test_seurat.csv")
SVM_rej_MNN <- read.csv("SVM_rej70_Time_PBMC_train_MCA_test_MNN.csv")

time_pbmc <- data.frame(time = c(SVM_counts$Training.time,
                                 SVM_seurat$Training.time,
                                 SVM_MNN$Training.time,
                                 SVM_rej_counts$Training.time,
                                 SVM_rej_seurat$Training.time,
                                 SVM_rej_MNN$Training.time),
                        processing = rep(c("Raw counts", "Seurat", "MNN"), 2),
                        method = c(rep("SVM", 3),
                                   rep("SVMrejection", 3))
)
time_pbmc$processing <- factor(time_pbmc$processing, levels = c("Raw counts", "Seurat", "MNN"))
time_pbmc$method <- factor(time_pbmc$method, levels = c("SVM", "SVMrejection"))

#####


# Plot with computation time for MCA-PBMCs ####

plot1 <- ggplot(time_mca, aes(method, time, fill = method)) +
  geom_col() +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x = element_blank(),
        plot.margin=unit(c(0.1,0.5,0,0), "cm")) +
  geom_text(aes(label= paste0( sprintf("%.2f", time), " s")),
            vjust=-1, color="black", size=2,
            position = position_dodge(.9)) +
  facet_wrap(~ processing, scales = "free_x", ncol = 4) +
  scale_fill_discrete(name = "Method") +
  xlab("Method") + ylab("Time (seconds)") +
  ylim(0, 200)

plot2 <- ggplot(time_pbmc, aes(method, time, fill = method)) +
  geom_col() +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x = element_blank()) +
  geom_text(aes(label= paste0( sprintf("%.2f", time), " s")),
            vjust=-1, color="black", size=2,
            position = position_dodge(.9)) +
  facet_wrap(~ processing, scales = "free_x", ncol = 4)  +
  scale_fill_discrete(name = "Method") +
  xlab("")  + ylab("Time (seconds)") +
  ylim(0, 200)


plot <- ggarrange(
  plot1, plot2,
  labels = c("A", "B"),
  font.label = list(size = 12),
  common.legend = TRUE, legend = "right", nrow = 2,
  heights = c(0.50, 0.5),
  align  = "v"
)

setwd("C:/Users/Usuario/Desktop/master/tfm/results/")
ggsave("computation_time_rej.pdf",plot, width = 10,height=11, dpi=300, units = "in")

#####


###########################################
# Assessment SVM with rejection option ####
###########################################


# 1) Reduced sample ####

setwd("C:/Users/Usuario/Desktop/master/tfm/HCL_data/")

pred_reducedSample <- read.csv("./results/SVM_ref_reduced_trainSample.csv")

setwd("C:/Users/Usuario/Desktop/master/tfm/HCL_data")

true_all <- read.csv("labels_mca_all.csv")

true_all <- txt_recode(true_all$x, from = names(table(true_all)),
                       to = c("B cells", "B cell(Plasmocyte)", "T cells", "T cells",
                              "DC", "Erythroid cell", "Lymphocyte", "Macrophages",
                              "Mast cells", "Myoloid cell", "Neutrophils", "T cells"))


SVM_reduced <- evaluation(true_all, pred_reducedSample$X0,
                          cells = c("T cells", "B cells", "DC",
                                    "Neutrophils", "Macrophages",
                                    "Mast cells", "B cell(Plasmocyte)",
                                    "Erythroid cell", "Lymphocyte",
                                    "Myoloid cell"),
                          n_correct = 6, n_cells = 10)
#####

# 3) SVM with rejection option using prefit option ####

true_all <- read.csv("labels_mca_all.csv")

true_all <- txt_recode(true_all$x, from = names(table(true_all)),
                       to = c("B cells", "Unknown", "T cells", "T cells",
                              "DC", "Unknown", "Unknown", "Macrophages",
                              "Mast cells", "Unknown", "Neutrophils", "T cells"))


pred_07 <- read.csv("./results/SVM_rej07.csv")
pred_03 <- read.csv("./results/SVM_rej03.csv")
pred_01 <- read.csv("./results/SVM_rej01.csv")




SVM_07 <- evaluation(true_all, pred_07$X0,
                     cells = c("T cells", "B cells", "DC",
                               "Neutrophils", "Macrophages",
                               "Mast cells", "Unknown"),
                     n_correct = 7, n_cells = 7)
table(pred_07$X0) / length(pred_07$X0)



SVM_03 <- evaluation(true_all, pred_03$X0,
                     cells = c("T cells", "B cells", "DC",
                               "Neutrophils", "Macrophages",
                               "Mast cells", "Unknown"),
                     n_correct = 7, n_cells = 7)

table(pred_03$X0)
table(pred_03$X0) / length(pred_03$X0)



SVM_01 <- evaluation(true_all, pred_01$X0,
                     cells = c("T cells", "B cells", "DC",
                               "Neutrophils", "Macrophages",
                               "Mast cells", "Unknown"),
                     n_correct = 7, n_cells = 7)

table(pred_01$X0)
table(pred_01$X0) / length(pred_01$X0)


#####


# 6) Random forest ####
pred_07 <- read.csv("./results/RF_07.csv")
pred_03 <- read.csv("./results/RF_03.csv")
pred_01 <- read.csv("./results/RF_01.csv")




RF_07 <- evaluation(true_all, pred_07$X0,
                    cells = c("T cells", "B cells", "DC",
                              "Neutrophils", "Macrophages",
                              "Mast cells", "Unknown"),
                    n_correct = 1, n_cells = 1)
table(pred_07$X0) / length(pred_07$X0)



RF_03 <- evaluation(true_all, pred_03$X0,
                    cells = c("T cells", "B cells", "DC",
                              "Neutrophils", "Macrophages",
                              "Mast cells", "Unknown"),
                    n_correct = 7, n_cells = 7)

table(pred_03$X0)
table(pred_03$X0) / length(pred_03$X0)



RF_01 <- evaluation(true_all, pred_01$X0,
                    cells = c("T cells", "B cells", "DC",
                              "Neutrophils", "Macrophages",
                              "Mast cells", "Unknown"),
                    n_correct = 7, n_cells = 7)

table(pred_01$X0)
table(pred_01$X0) / length(pred_01$X0)

#####
