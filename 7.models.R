# setup -------------------------------------------------------------------

library(tidyverse)
library(GET)
library(parallel)
library(patchwork)

source("0.functions.R")

RNGkind("L'Ecuyer-CMRG")
set.seed(24695)

# Initiate cluster
cl <- makeCluster(detectCores() - 1)
parallel::clusterSetRNGStream(cl = cl, iseed = 24695)

#load data
indexcurves <- readRDS("curves.rds")
vars <- readRDS("selected_variables.rds")


## FLM testing function ------------------------------------------------

flmtest <-  function (index, formula, nsim = 19999, cl=NULL){
  graph.flm(nsim = nsim,
            formula.full =  formula,
            formula.reduced = Y ~ 1,
            curve_sets = list(Y = index),
            factors = vars,
            cl = cl)
}

## plotting functions -------------------------------------------------------

#global envelope plot

plot_FLM <- function(x, title = NULL){
  plot <- plot(x)+
    xlab(expression(area~(italic(m^2))))
  subtitle <- gsub("Graphical functional GLM: ", "", plot$labels$title)
  plot[["layers"]][[1]][["aes_params"]]$fill <- rgb(188, 223, 235, maxColorValue = 255)
  plot[["layers"]][[1]][["aes_params"]]$alpha <- 1
  plot+
    theme_minimal()+
    labs(title = title,
         subtitle = subtitle)+
    theme(plot.subtitle = element_text(face = if(as.numeric(gsub("p [⁠[:punct:]] ", "", subtitle))>0.05)
      "plain" else "bold", size = 10),
      legend.position = "none",
      strip.text.x = element_text(size = 10))
}

#plot patchwork

patch <- function(x, ncol=2){
  wrap_plots(x, ncol = ncol, byrow = T, guides = "collect")+
    plot_layout(axis_titles = "collect")
}

# select indices from qdecomp output --------------------------------------

indices <- c("E_alpha", "E_gamma", "E_beta_mult", "redundancy_a", "U_gamma_star", "clustering")

#nicer labels for plotting
labels <- list(E_alpha = expression(E[alpha]),
               E_gamma = expression(E[gamma]),
               E_beta_mult = expression(E[beta]),
               redundancy_a = expression(R[alpha]^"*"),
               U_gamma_star = expression(U[gamma]^"*"),
               clustering = "Clustering")

indexcurves <- lapply(indexcurves, function (x) x[indices])


# functional linear models ------------------------------------------------

#main model
flm_FD <- lapply(indexcurves, function (x)
  lapply(x, function (y) flmtest(y, formula = formula(Y ~ RH050 + LAI), cl = cl))
)

#only  LAI
flm_LAI <- lapply(indexcurves, function (x)
  lapply(x, function (y) flmtest(y, formula = formula(Y ~ LAI), cl = cl))
)

#only RH050
flm_RH <- lapply(indexcurves, function (x)
  lapply(x, function (y) flmtest(y, formula = formula(Y ~ RH050), cl = cl))
)

#with quadratic terms
flm_quadratic <- lapply(indexcurves, function (x)
  lapply(x, function (y) flmtest(y, 
                                 formula = formula(Y ~ RH050 + LAI + I(RH050^2) + I(LAI^2)),
                                 cl = cl))
)

#save results
saveRDS(flm_FD, "FLM.rds")
saveRDS(flm_LAI, "FLM_LAI.rds")
saveRDS(flm_RH, "FLM_RH050.rds")
saveRDS(flm_quadratic, "FLM_quadratic.rds")

# #load results
# flm_FD <- readRDS("FLM.rds")
# flm_LAI <- readRDS("FLM_LAI.rds")
# flm_RH <- readRDS("FLM_RH050.rds")
# flm_quadratic <- readRDS("FLM_quadratic.rds")

# plotting results --------------------------------------------------------

## main model----

plots <- lapply(flm_FD, function (x)
  imap(x, ~ plot_FLM(.x, title = labels[[.y]]))
)

patch(plots$overall)
ggsave('plots/overall.png', width = 190, height = 117, units = "mm",
       bg = 'white', scale = 1.5, dpi = 1000)
ggsave('plots/FIG3.eps', width = 190, height = 117, units = "mm",
       bg = 'white', scale = 1.5, dpi = 1000)

patch(plots$aboveground)
ggsave('plots/above.png', width = 190, height = 117, units = "mm",
       bg = 'white', scale = 1.5, dpi = 1000)
ggsave('plots/FIG4.eps', width = 190, height = 117, units = "mm",
       bg = 'white', scale = 1.5, dpi = 1000)

patch(plots$clonal)
ggsave('plots/clonal.png', width = 190, height = 117, units = "mm",
       bg = 'white', scale = 1.5, dpi = 1000)
ggsave('plots/FIG5.eps', width = 190, height = 117, units = "mm",
       bg = 'white', scale = 1.5, dpi = 1000)

## single predictor models ----

### RH050 ----
plots_RH <- lapply(flm_RH, function (x)
  imap(x, ~ plot_FLM(.x, title = labels[[.y]]))
)

patch(plots_RH$overall)
ggsave('plots/RHoverall.png', width = 190, height = 117, units = "mm",
       bg = 'white', scale = 1.5, dpi = 1000)

patch(plots_RH$aboveground)
ggsave('plots/RHabove.png', width = 190, height = 117, units = "mm",
       bg = 'white', scale = 1.5, dpi = 1000)

patch(plots_RH$clonal)
ggsave('plots/RHclonal.png', width = 190, height = 117, units = "mm",
       bg = 'white', scale = 1.5, dpi = 1000)


### LAI ----
plots_LAI <- lapply(flm_LAI, function (x)
  imap(x, ~ plot_FLM(.x, title = labels[[.y]]))
)

patch(plots_LAI$overall)
ggsave('plots/LAIoverall.png', width = 190, height = 117, units = "mm",
       bg = 'white', scale = 1.5, dpi = 1000)

patch(plots_LAI$aboveground)
ggsave('plots/LAIabove.png', width = 190, height = 117, units = "mm",
       bg = 'white', scale = 1.5, dpi = 1000)

patch(plots_LAI$clonal)
ggsave('plots/LAIclonal.png', width = 190, height = 117, units = "mm",
       bg = 'white', scale = 1.5, dpi = 1000)

## quadratic model ----
plots_quad <- lapply(flm_quadratic, function (x)
  imap(x, ~ plot_FLM(.x, title = labels[[.y]]))
)

patch(plots_quad$overall)
ggsave('plots/quad_overall.png', width = 190, height = 117, units = "mm",
       bg = 'white', scale = 1.5, dpi = 1000)

patch(plots_quad$aboveground)
ggsave('plots/quad_above.png', width = 190, height = 117, units = "mm",
       bg = 'white', scale = 1.5, dpi = 1000)

patch(plots_quad$clonal)
ggsave('plots/quad_clonal.png', width = 190, height = 117, units = "mm",
       bg = 'white', scale = 1.5, dpi = 1000)


# functional axes models --------------------------------------------------

#load data
curves_axes <- readRDS("curves_axes.rds")

#select indices
curves_axes <- lapply(curves_axes, function (x) x[indices])

#model
flm_axes <- lapply(curves_axes, function (x)
  lapply(x, function (y) flmtest(y, formula = formula(Y ~ RH050 + LAI), cl = cl))
)

#save results
saveRDS(flm_axes, "FLM_axes.rds")

# #load results
# flm_axes <- readRDS("FLM_axes.rds")

#plot
plots_axes <- lapply(flm_axes, function (x)
  imap(x, ~ plot_FLM(.x, title = labels[[.y]]))
)

patch(plots_axes$above1)
ggsave('plots/aboveRC1.png', width = 190, height = 117, units = "mm",
       bg = 'white', scale = 1.5, dpi = 1000)

patch(plots_axes$above2)
ggsave('plots/aboveRC2.png', width = 190, height = 117, units = "mm",
       bg = 'white', scale = 1.5, dpi = 1000)

patch(plots_axes$clonal1)
ggsave('plots/clonalRC1.png', width = 190, height = 117, units = "mm",
       bg = 'white', scale = 1.5, dpi = 1000)

patch(plots_axes$clonal2)
ggsave('plots/clonalRC2.png', width = 190, height = 117, units = "mm",
       bg = 'white', scale = 1.5, dpi = 1000)


# supplementary indices ---------------------------------------------------

#load data
supp <- readRDS("curves.rds")

#select indices

indices_supp <- c("alpha", "gamma", "beta_add", "E_alpha_sp", "E_gamma_sp", "E_beta_sp")

#nicer labels for plotting
labels_supp <- list(alpha = expression(Q[alpha]),
                   gamma = expression(Q[gamma]),
                   beta_add = expression(Q[beta]),
                   E_alpha_sp = expression(E[alpha~tax]),
                   E_gamma_sp = expression(E[gamma~tax]),
                   E_beta_sp = expression(E[beta~tax]))

supp <- lapply(supp, function (x) x[indices_supp])

#models

flm_supp <- lapply(supp, function (x)
  lapply(x, function (y) flmtest(y, formula = formula(Y ~ RH050 + LAI), cl = cl))
)

#plot
plots_supp <- lapply(flm_supp, function (x)
  imap(x, ~ plot_FLM(.x, title = labels_supp[[.y]]))
)

patch(plots_supp$overall)
ggsave('plots/overall_supp.png', width = 190, height = 90, units = "mm",
       bg = 'white', scale = 1.5, dpi = 1000)

patch(plots_supp$above)
ggsave('plots/above_supp.png', width = 190, height = 90, units = "mm",
       bg = 'white', scale = 1.5, dpi = 1000)

patch(plots_supp$clonal)
ggsave('plots/clonal_supp.png', width = 190, height = 90, units = "mm",
       bg = 'white', scale = 1.5, dpi = 1000)
