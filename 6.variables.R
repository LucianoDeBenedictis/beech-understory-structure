# setup -------------------------------------------------------------------

library(tidyverse)
library(ggfortify)
library(fclust)
library(party)

set.seed(403054)

cv <- function(x, na.rm = FALSE) {
  sd(x, na.rm = na.rm) / mean(x, na.rm = na.rm)
}


# import data -------------------------------------------------------------

#DCP canopy indices
canopy <- read_csv("data/DCP.csv") |> 
  group_by(site) |> 
  summarise(LAImean = mean(LAI), LAIcv = cv(LAI)) |> 
  rename(LAI = LAImean)

#TLS structural indices
structure <- read_csv("data/TLS.csv") |> 
  mutate(links_norm = n_links_4_6/TopH,
         min_pc_norm = pc_min_sect/TopH,
         max_pc_norm = pc_max_sect/TopH,
         sd_perc_norm = sd_perc_clust/TopH) |> 
  select(!c("Sampling_date", "n_links_4_6", "pc_min_sect", "pc_max_sect", "sd_perc_clust"))

#join together
variables <- structure |> 
  inner_join(canopy)

rm(structure, canopy)

# join categories
cat <- read_csv("data/classification.csv")
variables <- cat |> 
  inner_join(variables) |>
  mutate(category = as.factor(category))
rm(cat)

# variable selection ------------------------------------------------------------------

#random forest

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
ggsave('Plots/FIG1.eps', width = 190, height = 100, units = "mm", bg = 'white',
       scale = 0.8, dpi = 1000)

#select first 3 variables
selected <- importance |> 
  slice_head(n = 3) |> 
  pull(variable)

sel_var <- variables |> 
  select(site, category, all_of(selected))

#standardization
#FLM test is non-parametric, same results either way
#this just improves visualization of coefficients and is required for clustering

sel_var <- sel_var |>
  mutate(across(where(is.numeric), scale)) |>
  mutate(across(where(is.numeric), as.vector))

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
  select(RH050, LAI, RH075) |> 
  prcomp(scale. = T, rank. = 2)

sel_var <- cbind(sel_var, pca_results$x)

#plot PCA
pcaplot <- autoplot(pca_results, label = F, loadings = T, loadings.label = T,
         loadings.label.repel = T, scale = 0.05)

pcaplot$layers[[1]] <- NULL

#add clustering and categories
pcaplot+
  geom_point(data = sel_var, aes(x = PC1, y = PC2, color = cluster_fuzzy,
                                  shape = category, size = membership)) +
  scale_size_area(max_size = 3) +
  scale_color_brewer(palette = "Set1")+
  theme_minimal() +
  labs(color = "Cluster", shape = "Category", size = "Membership")

ggsave("plots/cluster_fuzzy_select.png", width = 190, height = 117,
       units = "mm", bg = 'white', scale = 1, dpi = 1000)
ggsave('Plots/FIG2.eps', width = 190, height = 117,
       units = "mm", bg = 'white', scale = 1)

#summary statistics by cluster
sel_var |> 
  group_by(cluster_fuzzy) |> 
  summarise(across(3:5, mean))
