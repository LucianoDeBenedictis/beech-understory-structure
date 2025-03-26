# setup -------------------------------------------------------------------

library(tidyverse)
library(GET)
library(vegan)
#library(ggfortify)
library(psych)

resampled <- readRDS("transect_resampled.rds")
speciestraits <- readRDS("traits_imputation.rds")
source("0.functions.R")

indicesaxes <- vector(mode = "list")
curvesaxes <- vector(mode = "list")

# transform traits --------------------------------------------------------

data <- speciestraits |> 
  select(species, where(is.numeric)) |> 
  select(-n) |> 
  mutate(across(where(is.numeric),  ~ boxcox(.))) |> 
  as.data.frame()


# check dimensions --------------------------------------------------------

PCA_above <- data |> 
  select(2:7) |>
  psych::principal(nfactors = 2, rotate = "varimax")

PCA_above

png('plots/PCA_above.png', bg = "white", width = 7, height = 7, units = "in",
    res = 1000)
biplot(PCA_above, main = NULL)
dev.off()

data <- cbind(data, PCA_above$scores) |> 
  rename(above1 = RC1, above2 = RC2)

PCA_clonal <- data |> 
  select(8:11) |> 
  psych::principal(nfactors = 2, rotate = "varimax")

PCA_clonal

png('plots/PCA_clonal.png', bg = "white", width = 7, height = 7, units = "in",
    res = 1000)
biplot(PCA_clonal, main = NULL)
dev.off()

data <- cbind(data, PCA_clonal$scores) |> 
  rename(clonal1 = RC1, clonal2 = RC2)


# aboveground--------------------------------------------------------------

## aboveground RC1 -----
### calculate distances ----------------------------------------------------

distances <- data |> 
  select(above1) |> 
  dist(method = "euclidean", diag = T, upper = T) |> 
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

### calculate indices -----------------------------------------------------

species <- data$species

raotrait <- lapply(resampled, function (x)
  lapply(x, function (y) qdecomp(y, distances, species)) |> 
    bind_rows(.id = "step") |> 
    mutate(step = as.numeric(step))
) |> 
  bind_rows(.id = "site")

#scale steps in units of area

raotrait <- raotrait |> 
  mutate(step = step / 100)

#add to list
indicesaxes[[1]] <- raotrait
names(indicesaxes)[1] <- "above1"

### create curve set ------------------------------------------------------

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

curvesaxes[[1]] <- curvesets
names(curvesaxes)[1] <- "above1"

## aboveground RC2 -----
### calculate distances ----------------------------------------------------

distances <- data |> 
  select(above2) |> 
  dist(method = "euclidean", diag = T, upper = T) |> 
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

### calculate indices -----------------------------------------------------

species <- data$species

raotrait <- lapply(resampled, function (x)
  lapply(x, function (y) qdecomp(y, distances, species)) |> 
    bind_rows(.id = "step") |> 
    mutate(step = as.numeric(step))
) |> 
  bind_rows(.id = "site")

#scale steps in units of area

raotrait <- raotrait |> 
  mutate(step = step / 100)

#add to list
indicesaxes[[2]] <- raotrait
names(indicesaxes)[2] <- "above2"

### create curve set ------------------------------------------------------

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

curvesaxes[[2]] <- curvesets
names(curvesaxes)[2] <- "above2"

# clonal --------------------------------------------------------------

## clonal RC1 -----
### calculate distances ----------------------------------------------------

distances <- data |> 
  select(clonal1) |> 
  dist(method = "euclidean", diag = T, upper = T) |> 
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

### calculate indices -----------------------------------------------------

species <- data$species

raotrait <- lapply(resampled, function (x)
  lapply(x, function (y) qdecomp(y, distances, species)) |> 
    bind_rows(.id = "step") |> 
    mutate(step = as.numeric(step))
) |> 
  bind_rows(.id = "site")

#scale steps in units of area

raotrait <- raotrait |> 
  mutate(step = step / 100)

#add to list
indicesaxes[[3]] <- raotrait
names(indicesaxes)[3] <- "clonal1"

### create curve set ------------------------------------------------------

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

curvesaxes[[3]] <- curvesets
names(curvesaxes)[3] <- "clonal1"

## clonal RC2 -----
### calculate distances ----------------------------------------------------

distances <- data |> 
  select(clonal2) |> 
  dist(method = "euclidean", diag = T, upper = T) |> 
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

### calculate indices -----------------------------------------------------

species <- data$species

raotrait <- lapply(resampled, function (x)
  lapply(x, function (y) qdecomp(y, distances, species)) |> 
    bind_rows(.id = "step") |> 
    mutate(step = as.numeric(step))
) |> 
  bind_rows(.id = "site")

#scale steps in units of area

raotrait <- raotrait |> 
  mutate(step = step / 100)

#add to list
indicesaxes[[4]] <- raotrait
names(indicesaxes)[4] <- "clonal2"

### create curve set ------------------------------------------------------

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

curvesaxes[[4]] <- curvesets
names(curvesaxes)[4] <- "clonal2"


# save results ------------------------------------------------------------

saveRDS(indicesaxes, "indices_axes.rds")
saveRDS(curvesaxes, "curves_axes.rds")
