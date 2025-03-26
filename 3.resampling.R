# setup -------------------------------------------------------------------

dat <- readRDS("transect_data.rds")
source("0.functions.R")

# create vector of steps ----------------------------------------------------
#geometric progression, fixed ratio between steps

n <- 20 #desired steps
max <- 250 #maximum step

z <- n - 1
steps <- vector(length = 0)

while(length(steps) < n){
r <- max ^ (1 / z) #step ratio
steps <- unique(round(r^(seq(0,z)))) #keep unique after rounding
z <- z + 1 #add one step if less than desired
}

# resample ----------------------------------------------------------------

resampled <- lapply(dat, function (x) resample(x, steps, 1000))

saveRDS(resampled, file = "transect_resampled.rds")
