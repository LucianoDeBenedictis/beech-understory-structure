#' ---
#' title: "Gather soil and topography data"
#' author: "Luciano L.M. De Benedictis"
#' output: 
#'  pdf_document:
#'   latex_engine: xelatex
#'   keep_md: true
#' ---


# setup -------------------------------------------------------------------

library(tidyverse)
library(soilDB)
library(aqp)
library(terra)

points <- read_csv("map/coords.csv") |> 
  select(1:3) |> 
  rename(id = site, lat = Latitude, lon = Longitude)


# soil --------------------------------------------------------------------

#' Downloading SoilGrids data takes time, uncomment to run and save data.

# soil <- fetchSoilGrids(points, verbose = T)
# 
# saveRDS(soil, "soilgrids.rds")

soil <- readRDS("soilgrids.rds")

summary <- horizons(soil) |> 
  select(label, id, nitrogenmean, phh2omean) |> 
  filter(label %in% c("0-5", "5-15", "15-30")) |> 
  group_by(id) |> 
  summarise(nitrogen = mean(nitrogenmean), ph = mean(phh2omean))

sites <- read_csv("data/classification.csv") |> 
  left_join(summary, by = join_by(site == id)) |> 
  drop_na()

sites |> 
  pivot_longer(c(nitrogen, ph), names_to = "variable", values_to = "value") |> 
  ggplot(aes(x = category, y = value, color = category))+
  #geom_boxplot()+
  geom_point(position = 'jitterdodge')+
  facet_wrap(~variable, scales = 'free_y')


# topography --------------------------------------------------------------

dtm_casentino <- rast("map/DTM/w48570_s10/w48570_s10.tif")
plot(dtm_casentino)

v <- vect(points, crs = "WGS84")
plot(v)
WGS84 <- "+init=EPSG:4326"
dtm_casentino<- terra::project(dtm_casentino, WGS84) #reprojecting to epsg:4326#
plot(dtm_casentino)
points(v)

slope <- terrain(dtm_casentino, v = "slope", unit = "degrees")
aspect <- terrain(dtm_casentino, v = "aspect", unit = "degrees")

altitude <- extract(dtm_casentino, v, ID=F, bind=T) |> 
  values() |> 
  rename(site = 1, altitude = 2)

slope <- extract(slope, v, ID=F, bind=T) |> 
  values() |> 
  rename(site = 1, slope = 2)

aspect <- extract(aspect, v, ID=F, bind=T) |> 
  values() |> 
  rename(site = 1, aspect = 2)

sites <- sites |> 
  left_join(altitude) |> 
  left_join(slope) |> 
  left_join(aspect)

print(sites, n = 50)

write_csv(sites, "other results/sitetype.csv")

# summary -----------------------------------------------------------------

library(circular)

summary <- sites |> 
  mutate(aspect = as.circular(aspect, units = "degrees", type = "angles",
                              template = "none", modulo = "2pi", zero = 0,
                              rotation = "clock")) |> 
  group_by(category) |> 
  summarise(across(c(nitrogen:slope), list(mean = mean, sd = sd)),
            aspect_mean = mean.circular(aspect),
            aspect_sd = sd.circular(aspect))

glimpse(summary)

write_csv(summary, "other results/sitesummary.csv")

# session info ------------------------------------------------------------

sessionInfo()