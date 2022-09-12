
# Get SingleR annotations on JArribas dataset

# Packages ####
library(Seurat)
library(SeuratObject)
library(SingleR)
#####

# Import data ####
# Reference data
ref <- celldex::ImmGenData()

# Get processed test data
setwd("C:/Users/Usuario/Desktop/master/tfm/scripts_SingleR_merce/results")
dwIntegrated <- readRDS("dwIntegrated.20D.rds")

#####

# Function for SingleR annotations #####

SingleR_annotation <- function(test_data, ref_data, ref_labels,
                               setting){

  if(setting ==  "clusters"){

  pred_clust <- SingleR(test=as.SingleCellExperiment(test_data),
                        ref=ref_data,
                        clusters = test_data$seurat_clusters,
                        labels = ref_labels)
  new_cluster <- pred_clust$pruned.labels
  names(new_cluster) <- levels(test_data$seurat_clusters)
  test_data <- RenameIdents(test_data, new_cluster)
  return(Idents(test_data))
  }

  if(setting ==  "cells"){

    pred_cells <- SingleR(method = "single",
                          test=as.SingleCellExperiment(test_data),
                          ref=ref_data,
                          labels = ref_labels)
    return(pred_cells$pruned.labels)
  }
}

#####


# Obtaining annotations for each setting #####
# Setting 1a
SingleR_specific_cluster_SI_feature <- SingleR_annotation(dwIntegrated,
                                                          ref,
                                                          ref$label.fine,
                                                          setting = "clusters")

# Setting 1b
SingleR_specific_single_SI_feature <- SingleR_annotation(dwIntegrated,
                                                      ref,
                                                      ref$label.fine,
                                                      setting = "cells")


# Setting 2a
SingleR_main_cluster_SI_feature <- SingleR_annotation(dwIntegrated,
                                                      ref,
                                                      ref$label.main,
                                                      setting = "clusters")

# Setting 2b
SingleR_main_single_SI_feature <- SingleR_annotation(dwIntegrated,
                                                     ref,
                                                     ref$label.main,
                                                     setting = "cells")

# Setting 3a
DefaultAssay(dwIntegrated) <- "RNA"

SingleR_main_cluster_NO_feature <- SingleR_annotation(dwIntegrated,
                                                   ref,
                                                   ref$label.main,
                                                   setting = "clusters")

# Setting 3b
SingleR_main_single_NO_feature <- SingleR_annotation(dwIntegrated,
                                                     ref,
                                                     ref$label.main,
                                                     setting = "cells")
######


# Export annotations #####
setwd("C:/Users/Usuario/Desktop/master/tfm/scripts_SingleR_merce/results")
write.csv(SingleR_specific_cluster_SI_feature,
          "SingleR_specific_cluster_SI_feature.csv",
          row.names = F)
write.csv(SingleR_specific_single_SI_feature,
          "SingleR_specific_single_SI_feature.csv",
          row.names = F)
write.csv(SingleR_main_cluster_SI_feature,
          "SingleR_main_cluster_SI_feature.csv",
          row.names = F)
write.csv(SingleR_main_single_SI_feature,
          "SingleR_main_single_SI_feature.csv",
          row.names = F)
write.csv(SingleR_main_cluster_NO_feature,
          "SingleR_main_cluster_NO_feature.csv",
          row.names = F)
write.csv(SingleR_main_single_NO_feature,
          "SingleR_main_single_NO_feature.csv",
          row.names = F)
