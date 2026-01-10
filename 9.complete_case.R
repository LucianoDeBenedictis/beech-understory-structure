#' ---
#' title: "Complete-case analysis"
#' author: "Luciano L.M. De Benedictis"
#' output: 
#'  pdf_document:
#'   latex_engine: xelatex
#'   keep_md: true
#' ---

# setup -------------------------------------------------------------------

library(tidyverse)
library(GET)
library(vegan)
library(parallel)
library(patchwork)

resampled <- readRDS("transect_resampled.rds")
speciestraits <- readRDS("traits_noimp.rds")
source("0.functions.R")

indices <- vector(mode = "list")
indexcurves <- vector(mode = "list")

#update boxcox to support NAs
boxcox <- function(x, check = F) {
  
  x1 <- na.omit(x)
  
  if (any(x1 == 0)) x1 <- x1 + 1
  
  boxcox_result <- MASS::boxcox(x1 ~ 1, plotit = F)
  lambda <- boxcox_result$x[which.max(boxcox_result$y)]
  
  message("Box-Cox transform with lambda = ", lambda)
  
  if (check == T) return(lambda)
  else{
    if (lambda == 0) {
      x1 <- log(x1)
    } else {
      x1 <- (x1 ^ lambda - 1) / lambda
    }
    x[!is.na(x)] <- x1
    return(x)
  }
}

# keep complete traits ----------------------------------------------------

remove_overl <- speciestraits |> 
  filter(if_any(leaf_area_mm2:spread, is.na)) |> 
  pull(species)

remove_above <- speciestraits |> 
  filter(if_any(leaf_area_mm2:ssd_combined_mg_mm3, is.na)) |> 
  pull(species)

remove_clo <- speciestraits |> 
  filter(if_any(BBRsize:spread, is.na)) |> 
  pull(species)

resampled_overl <- resampled |> 
  lapply(function(x) lapply(x, function(y) filter(y, !species %in% remove_overl)))

resampled_above <- resampled |> 
  lapply(function(x) lapply(x, function(y) filter(y, !species %in% remove_above)))

resampled_clo <- resampled |> 
  lapply(function(x) lapply(x, function(y) filter(y, !species %in% remove_clo)))

rm(resampled)

# transform traits --------------------------------------------------------

speciestraits |> 
  select(-n) |> 
  summarize(across(where(is.numeric),  ~ boxcox(., check = T))) |> 
  pivot_longer(everything(), names_to = "trait", values_to = "lambda") |> 
  arrange(lambda, trait)

#' This is a dry run check. The traits will be transformed as such:
#' 
#' - ~ 1/sqrt N mass, LMA, offspring
#' - ~ log spread, SSD, leaf area, diaspore mass
#' - ~ sqrt height
#' - ~ identity BBR size
#' - ~ square persistence

data <- speciestraits |> 
  select(species, where(is.numeric)) |> 
  select(-n) |> 
  mutate(across(where(is.numeric),  ~ boxcox(.))) |> 
  as.data.frame()

# check trait distributions

data |> 
  pivot_longer(-species, names_to = "trait", values_to = "value") |> 
  histdensity(value)+
  facet_wrap(~ trait, scales = "free")+
  theme_minimal()

# overall -----------------------------------------------------------------

## calculate distances ----------------------------------------------------

distances <- data |> 
  filter(!species %in% remove_overl) |> 
  select(-1) |> 
  vegdist(method = "mahalanobis", diag = T, upper = T) |> 
  as.matrix()

#scaling to 0-1

distances <- distances / max(distances)

#double checking the squared euclidean property

ade4::is.euclid(as.dist(sqrt(distances)))

#check distance distribution

hist(distances)

#pairs with distance 0

as.data.frame(as.table(distances)) |> 
  filter(Freq == 0 & Var1 != Var2)

## calculate indices -------------------------------------------------------

species <- data |> 
  filter(!species %in% remove_overl) |> 
  pull(species)

raotrait <- lapply(resampled_overl, function (x)
  lapply(x, function (y) qdecomp(y, distances, species)) |> 
    bind_rows(.id = "step") |> 
    mutate(step = as.numeric(step))
) |> 
  bind_rows(.id = "site")

#scale steps in units of area

raotrait <- raotrait |> 
  mutate(step = step / 100)

#add to list
indices[[1]] <- raotrait
names(indices)[1] <- "overall"


## create curve set --------------------------------------------------------

curvesets <- vector(mode = "list", length = sum(map_lgl(raotrait, is.numeric))-1)

for (i in seq_along(curvesets)){
  column <- colnames(raotrait)[i+2]
  curves <- raotrait |> 
    select(site, step, {{column}}) |> 
    pivot_wider(names_from = site, values_from = 3) |> 
    filter(if_all(where(is.numeric), ~ !is.na(.)))
  curvesets[[i]] <- curve_set(obs = as.data.frame(curves[-1]), r = curves[[1]])
  names(curvesets)[i] <- column
}

indexcurves[[1]] <- curvesets
names(indexcurves)[1] <- "overall"


# aboveground -------------------------------------------------------------


## calculate distances ----------------------------------------------------

distances <- data |> 
  filter(!species %in% remove_above) |> 
  select(leaf_area_mm2:ssd_combined_mg_mm3) |> 
  vegdist(method = "mahalanobis", diag = T, upper = T) |> 
  as.matrix()

#scaling to 0-1

distances <- distances / max(distances)

#double checking the squared euclidean property

ade4::is.euclid(as.dist(sqrt(distances)))

#check distance distribution

hist(distances)

#pairs with distance 0

as.data.frame(as.table(distances)) |> 
  filter(Freq == 0 & Var1 != Var2)

## calculate indices -------------------------------------------------------

species <- data |> 
  filter(!species %in% remove_above) |> 
  pull(species)

raotrait <- lapply(resampled_above, function (x)
  lapply(x, function (y) qdecomp(y, distances, species)) |> 
    bind_rows(.id = "step") |> 
    mutate(step = as.numeric(step))
) |> 
  bind_rows(.id = "site")

#scale steps in units of area

raotrait <- raotrait |> 
  mutate(step = step / 100)

#add to list

indices[[2]] <- raotrait
names(indices)[2] <- "aboveground"


## create curve set --------------------------------------------------------

curvesets <- vector(mode = "list", length = sum(map_lgl(raotrait, is.numeric))-1)

for (i in seq_along(curvesets)){
  column <- colnames(raotrait)[i+2]
  curves <- raotrait |> 
    select(site, step, {{column}}) |> 
    pivot_wider(names_from = site, values_from = 3) |> 
    filter(if_all(where(is.numeric), ~ !is.na(.)))
  curvesets[[i]] <- curve_set(obs = as.data.frame(curves[-1]), r = curves[[1]])
  names(curvesets)[i] <- column
}

indexcurves[[2]] <- curvesets
names(indexcurves)[2] <- "aboveground"

# clonal -------------------------------------------------------------

## calculate distances ----------------------------------------------------
distances <- data |> 
  filter(!species %in% remove_clo) |> 
  select(BBRsize:spread) |> 
  vegdist(method = "mahalanobis", diag = T, upper = T) |> 
  as.matrix()

#scaling to 0-1

distances <- distances / max(distances)

#double checking the squared euclidean property

ade4::is.euclid(as.dist(sqrt(distances)))

#check distance distribution

hist(distances)

#pairs with distance 0

as.data.frame(as.table(distances)) |> 
  filter(Freq == 0 & Var1 != Var2)

## calculate indices -------------------------------------------------------

species <- data |> 
  filter(!species %in% remove_clo) |> 
  pull(species)

raotrait <- lapply(resampled_clo, function (x)
  lapply(x, function (y) qdecomp(y, distances, species)) |> 
    bind_rows(.id = "step") |> 
    mutate(step = as.numeric(step))
) |> 
  bind_rows(.id = "site")

#scale steps in units of area

raotrait <- raotrait |> 
  mutate(step = step / 100)

#add to list

indices[[3]] <- raotrait
names(indices)[3] <- "clonal"


## create curve set --------------------------------------------------------

curvesets <- vector(mode = "list", length = sum(map_lgl(raotrait, is.numeric))-1)

for (i in seq_along(curvesets)){
  column <- colnames(raotrait)[i+2]
  curves <- raotrait |> 
    select(site, step, {{column}}) |> 
    pivot_wider(names_from = site, values_from = 3) |> 
    filter(if_all(where(is.numeric), ~ !is.na(.)))
  curvesets[[i]] <- curve_set(obs = as.data.frame(curves[-1]), r = curves[[1]])
  names(curvesets)[i] <- column
}

indexcurves[[3]] <- curvesets
names(indexcurves)[3] <- "clonal"


rm(resampled_above, resampled_clo, resampled_overl)

# models ------------------------------------------------------------------

RNGkind("L'Ecuyer-CMRG")
set.seed(24695)

# Initiate cluster
cl <- makeCluster(detectCores())
parallel::clusterSetRNGStream(cl = cl, iseed = 24695)

#load data

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
    scale_x_continuous( # convert x-axis from area to length
      labels = function(x) x * 10,
      breaks = c(0.01, seq(from = 0.5, to = 2.5, by = 0.5)),
      name = expression("length of sampling units" ~ (italic(m)))
    )
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

## select indices from qdecomp output --------------------------------------

indices <- c("E_alpha", "E_gamma", "E_beta_mult", "redundancy_a", "U_gamma_star", "clustering")

#nicer labels for plotting
labels <- list(E_alpha = expression(E[alpha]),
               E_gamma = expression(E[gamma]),
               E_beta_mult = expression(E[beta]),
               redundancy_a = expression(R[alpha]^"*"),
               U_gamma_star = expression(U[gamma]^"*"),
               clustering = "Clustering")

indexcurves <- lapply(indexcurves, function (x) x[indices])

## fit models ------------------------------------------------

flm_FD <- lapply(indexcurves, function (x)
  lapply(x, function (y) flmtest(y, formula = formula(Y ~ RH050 + LAI), cl = cl))
)

## plotting results -------------------------------------------------------

plots <- lapply(flm_FD, function (x)
  imap(x, ~ plot_FLM(.x, title = labels[[.y]]))
)

#+ fig.height=6
patch(plots$overall)
ggsave('plots/overall_noimp.png', width = 190, height = 117, units = "mm",
       bg = 'white', scale = 1.5, dpi = 1000)

patch(plots$aboveground)
ggsave('plots/above_noimp.png', width = 190, height = 117, units = "mm",
       bg = 'white', scale = 1.5, dpi = 1000)

patch(plots$clonal)
ggsave('plots/clonal_noimp.png', width = 190, height = 117, units = "mm",
       bg = 'white', scale = 1.5, dpi = 1000)

# session info ------------------------------------------------------------

sessionInfo()