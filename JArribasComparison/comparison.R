# Comparison of SingleR and SVM annotation on JArribas dataset

# Packages: ####
library(dplyr)
library(ggplot2)
library(networkD3)
library(webshot)
library(kableExtra)
#####

# Obtain SingleR annotations: ####

# Setting 1a
true1a <- read.csv("SingleR_specific_cluster_SI_feature.csv")

true1a[grep("^B cells", true1a)] <- "B cells"
true1a[grep("^DC", true1a)] <- "DC"
true1a[grep("^Epithelial", true1a)] <- "Epithelial"
true1a[grep("^Fibroblasts", true1a)] <- "Fibroblasts"
true1a[grep("^Macrophages", true1a)] <- "Macrophages"
true1a[grep("^NK cells", true1a)] <- "NK cells"
true1a[grep("^Neutrophils", true1a)] <- "Neutrophils"
true1a[grep("^T cells", true1a)] <- "T cells"
true1a[grep("^Tgd", true1a)] <- "Tgd"
true1a[grep("^ILC", true1a)] <- "ILC"
true1a[grep("^Stem", true1a)] <- "Stem"


# Setting 1b
true1b <- read.csv("SingleR_specific_single_SI_feature.csv")

# Setting 2a
true2a <- read.csv("SingleR_main_cluster_SI_feature.csv")

# Setting 2b
true2b <- read.csv("SingleR_main_single_SI_feature.csv")

# Setting 3a
true3a <- read.csv("SingleR_main_cluster_NO_feature.csv")

# Setting 3b
true3a <- read.csv("SingleR_main_single_NO_feature.csv")

#####

# Obtain SVM annotations ####

# Setting 1a - 1b
pred1a <- read.csv("SVM_specific_SI_feature.csv")

pred1a[grep("^B cells", pred1a)] <- "B cells"
pred1a[grep("^Basophils", pred1a)] <- "Basophils"
pred1a[grep("^DC", pred1a)] <- "DC"
pred1a[grep("^Endothelial", pred1a)] <- "Endothelial"
pred1a[grep("^Eosinophils", pred1a)] <- "Eosinophils"
pred1a[grep("^Epithelial", pred1a)] <- "Epithelial"
pred1a[grep("^Fibroblasts", pred1a)] <- "Fibroblasts"
pred1a[grep("^ILC", pred1a)] <- "ILC"
pred1a[grep("^Macrophages", pred1a)] <- "Macrophages"
pred1a[grep("^Mast", pred1a)] <- "Mast"
pred1a[grep("^Microglia", pred1a)] <- "Microglia"
pred1a[grep("^Monocytes", pred1a)] <- "Monocytes"
pred1a[grep("^NK cells", pred1a)] <- "NK cells"
pred1a[grep("^NKT", pred1a)] <- "NKT"
pred1a[grep("^Neutrophils", pred1a)] <- "Neutrophils"
pred1a[grep("^Stem cells", pred1a)] <- "Stem cells"
pred1a[grep("^Stromal", pred1a)] <- "Stromal"
pred1a[grep("^T cells", pred1a)] <- "T cells"
pred1a[grep("^Tgd", pred1a)] <- "Tgd"

pred1b <- pred1a


# Setting 2a - 2b
pred2a <- read.csv("SVM_main_SI_feature.csv")
pred2b <- read.csv("SVM_main_SI_feature.csv")

# Setting 3a - 3b
pred3a <- read.csv("SVM_main_NO_feature.csv")
pred3b <- read.csv("SVM_main_NO_feature.csv")
#####


# Barplot with % of global misclassification: ####

true_pred <- data.frame(true1a$x, pred1a$x,
                        true1b$x, pred1b$x,
                        true2a$x, pred2a$x,
                        true2b$x, pred2b$x,
                        true3a$x, pred3a$x,
                        true3b$x, pred3b$x)

mis_percentage <- rep(0, 6)
index <- seq(1, 11, by = 2)

for(i in 1:6){
  a <- index[i]
  mis_percentage[i] <- sum(true_pred[,a] != true_pred[,a + 1], na.rm = TRUE) /
    nrow(true_pred) * 100
}


mis_percentage_df <- data.frame(setting = c("Setting 1a", "Setting 1b",
                                            "Setting 2a", "Setting 2b",
                                            "Setting 3a", "Setting 3b"),
                                perc = mis_percentage)


p_global<-ggplot(data=mis_percentage_df, aes(x=setting, y=perc/100, fill = setting)) +
  geom_col() +
  theme_minimal() +
  theme(legend.position="none",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x = element_blank()) +
  scale_y_continuous(labels=scales::percent, limits = c(0,1)) +
  scale_fill_manual(name = "Setting", values = c("#1976D2", "#64B5F6",
                                                 "#689F38", "#9CCC65",
                                                 "#FFA000", "#FFCA28")) +
  geom_text(aes(label= paste0( sprintf("%.2f", perc), "%")), vjust=-1,
            color="black", size=3) +
  ylab("Percentatge of misclassification")
#####



# Barplot with % of misclassification per each cell type: ####
cell_types <- unique(c(true1a, true1b, true2a, true2b, true3a, true3b))
cell_types <- cell_types[!is.na(cell_types)]

mis_percentage_cellTypes <- list()

for(i in 1:6){
  a <- index[i]
  namesCellTypesTrue <- unique(true_pred[,a])
  namesCellTypesPred <- unique(true_pred[,a + 1])
  namesCellTypesBoth <- namesCellTypesTrue[namesCellTypesTrue %in% namesCellTypesPred]

  conf_matrix <- table(true_pred[,a], true_pred[,a + 1])
  n_correct_cellTypes <- sapply(namesCellTypesBoth, function(x) conf_matrix[x,x])
  n_cellTypes <- table(true_pred[,a])[sort(namesCellTypesTrue) %in% namesCellTypesBoth]

  n_correct_cellTypes <- n_correct_cellTypes[order(names(n_correct_cellTypes))]
  n_cellTypes <- n_cellTypes[order(names(n_cellTypes))]

  mis_percentage_cellTypes[[i]] <- (n_cellTypes - n_correct_cellTypes) / n_cellTypes * 100

  cells_not_predict <- namesCellTypesTrue[!namesCellTypesTrue %in% namesCellTypesPred]
  n_cells_not_predict <- rep(100, length(cells_not_predict))
  names(n_cells_not_predict) <- cells_not_predict

  mis_percentage_cellTypes[[i]] <- c(mis_percentage_cellTypes[[i]], n_cells_not_predict)

  cells_not_true <- cell_types[!cell_types %in% namesCellTypesTrue]
  n_cells_not_true <- rep(NA, length(cells_not_true))
  names(n_cells_not_true) <- cells_not_true

  mis_percentage_cellTypes[[i]] <- c(mis_percentage_cellTypes[[i]], n_cells_not_true)

  mis_percentage_cellTypes[[i]] <- mis_percentage_cellTypes[[i]][sort(names(mis_percentage_cellTypes[[i]]))]
}


n_cellTypes <- length(mis_percentage_cellTypes[[1]])
names_cellTypes <- names(mis_percentage_cellTypes[[1]])


aux <- matrix(unlist(mis_percentage_cellTypes), ncol = 6)
rownames(aux) <- names_cellTypes

n_na <- apply(aux, 1, function(x) sum(is.na(x)))
aux <- aux[n_na == 0, ]
names_aux <- rownames(aux)

mis_percentage_cellTypes_df <- data.frame(
  name = rep(names_aux, 6),
  percentage = c(aux),
  setting = rep(c("Setting 1a", "Setting 1b",
                  "Setting 2a", "Setting 2b",
                  "Setting 3a", "Setting 3b"), each = nrow(aux) ))



p_cells <- ggplot(mis_percentage_cellTypes_df, aes(setting, percentage/100, fill = setting)) +
  geom_col(position = position_dodge(.9)) +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x = element_blank(),
        legend.position="none") +
  scale_y_continuous(labels=scales::percent, limits = c(0,1)) +
  facet_wrap(~ name, scales = "free_x", ncol = 4) +
  geom_text(aes(label= paste0( sprintf("%.2f", percentage), "%")),
            vjust=-1, color="black", size=2,
            position = position_dodge(.9)) +
  scale_fill_manual(name = "Setting", values = c("#1976D2", "#64B5F6",
                                                 "#689F38", "#9CCC65",
                                                 "#FFA000", "#FFCA28")) +
  xlab("Cell type") +
  ylab("Percentatge of misclassification")

ggsave("C:/Users/Usuario/Desktop/master/tfm/results/figures_JArribas/p_global.pdf",
       p_global, width = 15,height=8, dpi=300, units = "in")
ggsave("C:/Users/Usuario/Desktop/master/tfm/results/figures_JArribas/p_cells.jpg",
       p_cells, width = 15,height=10, dpi=300, units = "in")

#####


# SANKEY DIAGRAM: #####

sankey <- function(true, pred){

  common_cell_types <- true %in% rownames(aux)

  true <- true[common_cell_types]
  pred <- pred[common_cell_types]
  pred <- paste(pred, " ")

  dt = data.frame(true, pred)

  names(dt) <- make.names(names(dt))
  dt <- dt %>%
    group_by(.dots=names(dt)) %>%
    summarise(count= n())

  nodes <- data.frame(name=c(as.character(dt$true),
                             as.character(dt$pred)) %>% unique())

  # With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
  dt$IDtrue=match(dt$true, nodes$name)-1
  dt$IDpred=match(dt$pred, nodes$name)-1

  # prepare colour scale
  ColourScal ='d3.scaleOrdinal() .
  range(["#FDE725FF","#B4DE2CFF","#6DCD59FF","#35B779FF",
  "#1F9E89FF","#26828EFF","#31688EFF","#3E4A89FF","#482878FF","#440154FF",
  "#FDE725DF","#B4DE2GFF","#6DCD58FF","#35B777FF",
  "#1F9E89FH","#26828EFB","#31688EMF","#3E4A87FF"])'

  dt$links$true <- sub(' .*', '',
                       dt$nodes[dt$links$true + 1, 'name'])

  # Make the Network
  p <- sankeyNetwork(Links = as.data.frame(dt), Nodes = nodes,
                     Source = "IDtrue", Target = "IDpred",
                     Value = "count", NodeID = "name",
                     sinksRight=TRUE, nodeWidth=30, fontSize=15,
                     nodePadding=20, width = 500, height = 800,
                     margin = list("left"=100), iterations = 0,
                     LinkGroup = 'true')


  return(p)
}

plot <- list()
for(i in 1:6){
  a <- index[i]
  plot[[i]] <- sankey(true_pred[,a], true_pred[,a + 1])
}


plot[[1]] %>%
  saveNetwork(file = 'C:/Users/Usuario/Desktop/master/tfm/results/figures_JArribas/setting1a.html')

plot[[2]] %>%
  saveNetwork(file = 'C:/Users/Usuario/Desktop/master/tfm/results/figures_JArribas/setting1b.html')

plot[[3]] %>%
  saveNetwork(file = 'C:/Users/Usuario/Desktop/master/tfm/results/figures_JArribas/setting2a.html')


plot[[4]] %>%
  saveNetwork(file = 'C:/Users/Usuario/Desktop/master/tfm/results/figures_JArribas/setting2b.html')

plot[[5]] %>%
  saveNetwork(file = 'C:/Users/Usuario/Desktop/master/tfm/results/figures_JArribas/setting3a.html')

plot[[6]] %>%
  saveNetwork(file = 'C:/Users/Usuario/Desktop/master/tfm/results/figures_JArribas/setting3b.html')



webshot('C:/Users/Usuario/Desktop/master/tfm/results/figures_JArribas/setting1a.html',
        'C:/Users/Usuario/Desktop/master/tfm/results/figures_JArribas/setting1a.png',
        vwidth = 500, vheight = 800)

webshot('C:/Users/Usuario/Desktop/master/tfm/results/figures_JArribas/setting1b.html',
        'C:/Users/Usuario/Desktop/master/tfm/results/figures_JArribas/setting1b.png',
        vwidth = 500, vheight = 800)

webshot('C:/Users/Usuario/Desktop/master/tfm/results/figures_JArribas/setting2a.html',
        'C:/Users/Usuario/Desktop/master/tfm/results/figures_JArribas/setting2a.png',
        vwidth = 500, vheight = 800)

webshot('C:/Users/Usuario/Desktop/master/tfm/results/figures_JArribas/setting2b.html',
        'C:/Users/Usuario/Desktop/master/tfm/results/figures_JArribas/setting2b.png',
        vwidth = 500, vheight = 800)

webshot('C:/Users/Usuario/Desktop/master/tfm/results/figures_JArribas/setting3a.html',
        'C:/Users/Usuario/Desktop/master/tfm/results/figures_JArribas/setting3a.png',
        vwidth = 500, vheight = 800)

webshot('C:/Users/Usuario/Desktop/master/tfm/results/figures_JArribas/setting3b.html',
        'C:/Users/Usuario/Desktop/master/tfm/results/figures_JArribas/setting3b.png',
        vwidth = 500, vheight = 800)
#####


# Scheme - table ####
setting = data.frame(
  Symbol=c("$\\square$", "$\\square$", "$\\square$",
           "$\\square$", "$\\square$", "$\\square$"),
  Color = c("#1976d2", "#64b5f6",
            "#689f38", "#9ccc65",
            "#ffa000", "#ffca28"),
  "Setting" = c("Setting 1a", "Setting 1b", "Setting 2a", "Setting 2b",
                "Setting 3a", "Setting 3b"),
  "Cell type label" = c(rep("Specific", 2), rep("Main", 4)),
  "Feature selection" = c(rep("Yes", 4), rep("No", 2)),
  "SingleR configuration" = c(rep(c("Per clusters", "Per cells"), 3))
)

colnames(setting) <- c("Symbol", "Color", " ", "Cell type label", "Feature selection",
                       "SingleR configuration")

ktable <- setting %>%
  mutate(Symbol = cell_spec(Symbol, color = Color, format = "latex",
                            escape = FALSE)) %>%
  select(-Color) %>%
  knitr::kable(escape = FALSE, booktabs = TRUE, format = "tex") %>%
  kable_styling(latex_options = c("striped", "scale_down")) %>%
  save_kable(file = "table_2.pdf", zoom = 1.5)
#####
