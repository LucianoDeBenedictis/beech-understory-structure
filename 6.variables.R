#' ---
#' title: "Structural and canopy variables"
#' author: "Luciano L.M. De Benedictis"
#' output: 
#'  pdf_document:
#'   latex_engine: xelatex
#'   keep_md: true
#' ---

# setup -------------------------------------------------------------------

library(tidyverse)
library(GGally)
library(klaR)
library(fclust)
library(party)

set.seed(92538)

# function for the upper triangle of pairs plot. Adapted from Solomon Kurz.

my_upper <- function(data, mapping, ...) {
  
  # get the x and y data to use the other code
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  r  <- unname(cor.test(x, y)$estimate)
  rt <- format(r, digits = 2)[1]
  tt <- as.character(rt)
  
  # plot the cor value
  ggally_text(
    label = tt, 
    mapping = aes(),
    size = 4)+
    theme_void()
}

# import data -------------------------------------------------------------

# DCP canopy indices
canopy <- read_csv("data/DCP.csv") |> 
  group_by(site) |> 
  summarise(LAImean = mean(LAI), LAIsd = sd(LAI)) |> 
  rename(LAI = LAImean)

# TLS structural indices
structure <- read_csv("data/TLS.csv") |> 
  mutate(links_norm = n_links_4_6/TopH,
         min_pc_norm = pc_min_sect/TopH,
         max_pc_norm = pc_max_sect/TopH,
         sd_perc_norm = sd_perc_clust/TopH) |> 
  dplyr::select(!c("Sampling_date", "n_links_4_6", "pc_min_sect", "pc_max_sect", "sd_perc_clust"))

# join together
variables <- structure |> 
  inner_join(canopy)

rm(structure, canopy)

# join categories
cat <- read_csv("data/classification.csv")
variables <- cat |> 
  inner_join(variables) |>
  mutate(category = as.factor(category))
rm(cat)

variables |> 
  dplyr::select(2:7, 12, 13) |> 
  print(n = 30)

#pairs plot
variables[3:13] |>
  ggpairs(upper = list(continuous =  my_upper))

# standardization
# FLM test is non-parametric, same results either way
# but standardization gives more meaningful coefficients and is required for other steps here

variables <- variables |> 
  mutate(across(where(is.numeric), scale))|>
  mutate(across(where(is.numeric), as.vector))

# random forest ------------------------------------------------------------------

forest <- cforest(category ~ ., data = variables[, -1], 
                  control = cforest_control(mtry = 3, ntree = 5000, mincriterion = 0)) #mtry default from randomForest

#variable importance from RF model
importance <- varimp(forest, conditional=T)

importance <- tibble(variable = names(importance), importance = importance) |> 
  mutate(variable = fct_reorder(variable, importance, .desc = T)) |> 
  arrange(variable)

#plot IV
importance |> 
  ggplot(aes(x = variable, y = importance, fill = variable))+
  geom_col(width = 0.5)+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank())+
  labs(x = NULL, y = "Importance")
ggsave("plots/importance.png", width = 190, height = 100, units = "mm", bg = 'white',
       scale = 0.8, dpi = 1000)
ggsave('plots/FIG1.eps', width = 190, height = 100, units = "mm", bg = 'white',
       scale = 0.8, dpi = 1000)

sel_var <- variables |> 
  dplyr::select(site, category, RH050, LAI)

saveRDS(sel_var, "selected_variables.rds")

# clustering --------------------------------------------------------------

# fuzzy k-means
fuzzy <- FKM(sel_var[-c(1, 2)], k = 3)

#add to df
sel_var <- sel_var |> 
  mutate(cluster_fuzzy = fuzzy$clus[,1],
         membership = fuzzy$clus[,2])

#assign the most plausible labels
contingency_fuzzy <- as.matrix(table(sel_var$category, sel_var$cluster_fuzzy))
assignment_fuzzy <- clue::solve_LSAP(t(contingency_fuzzy), maximum = T)

sel_var <- sel_var |> 
  mutate(cluster_fuzzy = assignment_fuzzy[cluster_fuzzy]) |> 
  mutate(cluster_fuzzy = factor(cluster_fuzzy, labels = levels(category)))

#check agreement

table(sel_var$category, sel_var$cluster_fuzzy)

#fuzzy Rand Index
RI.F(VC = sel_var$category, U = fuzzy$U)

sink("other results/external_validation.txt")
"Fuzzy Rand Index"
RI.F(VC = sel_var$category, U = fuzzy$U)
sink()

#perform PCA after clustering, only for 2D plotting

pca_results <- sel_var |> 
  dplyr::select(RH050, LAI) |> 
  prcomp(scale. = T, rank. = 2)

sel_var <- cbind(sel_var, pca_results$x)

#plot PCA

#adapted from ggbiplot source code

#variance explained
expl_var <- round(pca_results$sdev^2/sum(pca_results$sdev^2)*100, 2)

labels <- str_c(c("PC1 ", "PC2 "), 
                str_c("(", expl_var, "%)"))

#how much to shift variable labels
lbl_shift <- 0.1

#loadings and labels
loadings_tbl <- as_tibble(pca_results$rotation, rownames = "var") |> 
  mutate(xlab = PC1 + lbl_shift,
         ylab = PC2 + c(lbl_shift, -lbl_shift))

#default scaler from ggbiplot
scaler <- min(max(abs(sel_var$PC1))/max(abs(loadings_tbl$PC1)), 
              max(abs(sel_var$PC2))/max(abs(loadings_tbl$PC2)))*0.8

loadings_tbl <- loadings_tbl |> 
  mutate(across(PC1:ylab, ~ .x * scaler))

#clustering + PCA plot
sel_var |> 
  ggplot()+
  geom_point(aes(x = PC1, y = PC2, color = cluster_fuzzy,
                                 shape = category, size = membership)) +
  geom_segment(data = loadings_tbl, mapping = aes(x = 0, y = 0, xend = PC1,
                                                  yend = PC2), 
               arrow = grid::arrow(length = grid::unit(8, "points")), colour = "red", 
               linewidth = 0.5)+
  geom_text(data = loadings_tbl, aes(x = xlab, y = ylab, label = var), colour = "red")+
  scale_size_area(max_size = 3) +
  scale_color_brewer(palette = "Set1")+
  theme_minimal() +
  labs(color = "Cluster", shape = "Category", size = "Membership",
       x = labels[1], y = labels[2])

ggsave("plots/cluster_fuzzy_select.png", width = 190, height = 117,
       units = "mm", bg = 'white', scale = 1, dpi = 1000)
ggsave('plots/FIG2.eps', width = 190, height = 117,
       units = "mm", bg = 'white', scale = 1)

#summary statistics by cluster
sel_var |> 
  group_by(cluster_fuzzy) |> 
  summarise(across(3:5, mean))

# LDA ---------------------------------------------------------------------

# vector of classes
class <- variables[[2]]

# standardize variables
vars <- variables[3:13]

# variable selection, criterion "ability to separate"
step <- stepclass(x = vars, grouping = class, method = "lda", direction = "both",
                  criterion = "AS", improvement = 0.1)

wilks <- greedy.wilks(X = vars, grouping = class, niveau = 1) # keep going until all variables selected

write_csv(format(wilks$results, digits = 2), file ="other results/wilks_selection.csv")

# LDA with selected variables

selected <- vars[step[["model"]][["nr"]]]

lda <- lda(x = selected, grouping = class)

ggbiplot::ggbiplot(lda, groups = class, ellipse = T, varname.size = 4,
                   varname.adjust = 1.5,
                   ellipse.linewidth = NA, ellipse.alpha = 0.15, scale = 1)+
  labs(x = "LD1", y = "LD2", color = "Category", fill = "Category")+
  theme_minimal()+
  scale_color_brewer(palette = "Set1")+
  scale_fill_brewer(palette = "Set1")

ggsave("plots/LDA.png", width = 190, height = 117,
       units = "mm", bg = 'white', scale = 1.1, dpi = 1000)

drawparti(x = selected[[1]], y = selected[[2]], grouping = class, method = "lda",
          xlab = names(selected)[1], ylab = names(selected)[2], imageplot = F, col.mean = NULL)

png("plots/LDA_parti.png", bg = "white", width = 190, height = 117, units = "mm",
    pointsize = 12, res = 1000)
drawparti(x = selected[[1]], y = selected[[2]], grouping = class, method = "lda",
          xlab = names(selected)[1], ylab = names(selected)[2], imageplot = F, col.mean = NULL)
title("Partition plot")
dev.off()


# session info ------------------------------------------------------------

sessionInfo()