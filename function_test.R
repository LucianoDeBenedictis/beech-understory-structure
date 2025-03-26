library(tidyverse)

source("0.functions.R")

# checking Rao with de Bello et al. (2010) ----
## figure 1 ----

distances <- matrix(c(0,1,2,1,0,2,2,2,0), nrow = 3, ncol = 3,
                    dimnames = list(letters[1:3],letters[1:3]))

distances

props <- c(0.25,0.25,0.5)

data <- data.frame(species = letters[1:3],
                   site1 = c(1, 1, 2),
                   site2 = c(28, 1, 1),
                   site3 = c(1, 1, 2))

#checking raoq
props <- apply(data[2:4], 2, function (x) x/sum(x))

raos <- apply(props, 2, function (x) raoq(props=x, dist=distances))

raos
mean(raos)

#checking qdecomp

#transform to qdecomp long format
data <- data |> 
  pivot_longer(site1:site3, values_to = "abund", names_to = "site", names_prefix = "site") |> 
  uncount(abund) |> 
  rename("X" = "site")

decomp <- qdecomp(data = data, distances = distances, species = unique(data$species))
decomp

## figure 2 case 1 ----
forms <- rep(letters[1:4], each = 2)
distances <- outer(forms, forms, function(x, y) ifelse(x == y, 0, 1))

data <- data.frame(species = letters[1:8],
                   X = rep(1:2, each = 4))

decomp <- qdecomp(data, distances, species = letters[1:8])
decomp

## figure 2 case 2 ----
forms <- rep(letters[1:6], times = rep(c(2, 1, 1), 2))
distances <- outer(forms, forms, function(x, y) ifelse(x == y, 0, 1))

#same data as case 1

decomp <- qdecomp(data, distances, species = letters[1:8])
decomp

# simple resampling test ----
steps <- c(1,2,4)
df <- tibble(x = 1:20, species = letters[1:20])

resample_test <- resample(df, steps, 20)
