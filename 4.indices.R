# setup -------------------------------------------------------------------

library(tidyverse)
library(GET)
library(vegan)

resampled <- readRDS("transect_resampled.rds")
speciestraits <- readRDS("traits_imputation.rds")
source("0.functions.R")

indices <- vector(mode = "list")
indexcurves <- vector(mode = "list")

# transform traits --------------------------------------------------------

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
#ggsave("plots/traits_transformed.png", bg = "white")

# overall -----------------------------------------------------------------

## calculate distances ----------------------------------------------------

distances <- vegdist(data[-1], method = "mahalanobis", diag = T, upper = T) |> 
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

distances <- vegdist(data[2:7], method = "mahalanobis", diag = T, upper = T) |> 
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

distances <- vegdist(data[8:11], method = "mahalanobis", diag = T, upper = T) |> 
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


# plotting indices --------------------------------------------------------
categories <- read_csv("data/classification.csv")

indices$overall |> 
  left_join(categories) |> 
  pivot_longer(3:21, names_to = "index", values_to = "value") |> 
  ggplot(aes(x = step, y = value, group = site, colour = category))+
  geom_line()+
  facet_wrap(~index, scales = "free_y")

# save results ------------------------------------------------------------

saveRDS(indices, "indices.rds")
saveRDS(indexcurves, "curves.rds")
