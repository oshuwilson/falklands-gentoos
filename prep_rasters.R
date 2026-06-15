#-------------------------------------------------------------------------------
# Create Rasters for Habitat Mapping
# Joshua Wilson
#-------------------------------------------------------------------------------

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Placement/Writeup/")

library(tidyverse)
library(terra)
library(sf)
library(tidyterra)

# read in depth file
depth <- rast("data/rasters/ENV_Depths_FIzones2.tif")
depth[depth == 0] <- NA
plot(depth)

# create slope raster from depth raster
slope <- terrain(depth, v = "slope", neighbors = 4, unit = "degrees")
plot(slope)

# read in template 0.05 degree raster
temp <- rast("data/rasters/CB1213Gdive_densities.tif")
temp <- rast(ext(depth %>% project("EPSG:4326")), res = res(temp), crs = crs(temp))

#-------------------------------------------------------------------------------
# Format Rasters
#-------------------------------------------------------------------------------

# resample depth and slope to target resolution
depth <- project(depth, crs(temp)) %>%
  resample(temp)
slope <- project(slope, crs(temp)) %>%
  resample(temp)

# read in Falklands shapefile
fk <- read_sf("data/falklands_shapefile/All_Falkls.shp") %>%
  vect()

# create coordinates for colonies of interest in decimal degrees
bullroad <- data.frame(x = -59.386580, y = -52.309625, name = "BR") %>%
  vect(geom = c("x", "y"), crs = "epsg:4326")
cowbay <- data.frame(x=-57.869319,y=-51.428663) %>%
  vect(geom = c("x", "y"), crs = "epsg:4326")
steeple <- data.frame(x=-61.214336,y=-51.035945) %>%
  vect(geom = c("x", "y"), crs = "epsg:4326")

# mask out depth raster
depth <- mask(depth, fk, inverse = T)
plot(depth)

# mask out slope
slope <- mask(slope, fk, inverse = T)
plot(slope)

# create dist2coast raster
fk_rast <- rasterize(fk, depth, field = 1)
fk_rast[is.na(fk_rast)] <- 0
dist2coast <- gridDist(fk_rast, target = 1) %>%
  mask(fk, inverse = T)
plot(dist2coast) 

# create cost raster to allow dist2colony computation
costrast <- depth
costrast[!is.na(costrast)] <- 1
costrast[is.na(costrast)] <- 2
plot(costrast)

# get cell numbers that contain each colony
br_no <- extract(costrast, bullroad, cells = T) %>%
  pull(cell)
cb_no <- extract(costrast, cowbay, cells = T) %>%
  pull(cell)
sj_no <- extract(costrast, steeple, cells = T) %>%
  pull(cell)

# change value of cells that contain colonies to 3
br_cost <- costrast
br_cost[br_no] <- 3

cb_cost <- costrast
cb_cost[cb_no] <- 3

sj_cost <- costrast
sj_cost[sj_no] <- 3

# make land values NA
br_cost[br_cost == 2] <- NA
cb_cost[cb_cost == 2] <- NA
sj_cost[sj_cost == 2] <- NA

# compute distance to each colony using the land mask
br_dist <- costDist(br_cost, target = 3)
plot(br_dist)

cb_dist <- costDist(cb_cost, target = 3)
plot(cb_dist)

sj_dist <- costDist(sj_cost, target = 3)
plot(sj_dist)

# export all rasters
writeRaster(depth, "output/habmod/rasters/depth.tif", overwrite = T)
writeRaster(slope, "output/habmod/rasters/slope.tif", overwrite = T)
writeRaster(dist2coast, "output/habmod/rasters/dist2coast.tif", overwrite = T)
writeRaster(br_dist, "output/habmod/rasters/dist2BR.tif", overwrite = T)
writeRaster(cb_dist, "output/habmod/rasters/dist2CB.tif", overwrite = T)
writeRaster(sj_dist, "output/habmod/rasters/dist2SJ.tif", overwrite = T)


#-------------------------------------------------------------------------------
# Dynamic Predictors
#-------------------------------------------------------------------------------

# get unique dates of tracks
br <- read.csv("data/tracks/bull_roads.csv") %>%
  filter(breed_stage == "brood_guard")
cb <- read.csv("data/tracks/cow_bay.csv") %>%
  filter(breed_stage == "brood_guard")
sj <- read.csv("data/tracks/steeple_jason.csv")

br_dates <- unique(as_date(br$date_gmt, format = "%d/%m/%Y"))
cb_dates <- unique(as_date(cb$date_gmt, format = "%d/%m/%Y"))
sj_dates <- unique(as_date(sj$date_gmt, format = "%d/%m/%Y"))
all_dates <- unique(c(br_dates, cb_dates, sj_dates))
c(8, 10, 6, 13, 18)

# read in SST and Chlorophyll
sst12 <- rast("E:/Satellite_Data/daily/sst/sst_2012.nc")
sst13 <- rast("E:/Satellite_Data/daily/sst/sst_2013.nc")
sst14 <- rast("E:/Satellite_Data/daily/sst/sst_2014.nc")
chl12 <- rast("E:/Satellite_Data/daily/chl/resampled/chl_2012_resampled.nc")
chl13 <- rast("E:/Satellite_Data/daily/chl/resampled/chl_2013_resampled.nc")
chl14 <- rast("E:/Satellite_Data/daily/chl/resampled/chl_2014_resampled.nc")

# limit rasters to dates
sst12 <- sst12[[as_date(time(sst12)) %in% all_dates]]
sst13 <- sst13[[as_date(time(sst13)) %in% all_dates]]
sst14 <- sst14[[as_date(time(sst14)) %in% all_dates]]
chl12 <- chl12[[as_date(time(chl12)) %in% all_dates]]
chl13 <- chl13[[as_date(time(chl13)) %in% all_dates]]
chl14 <- chl14[[as_date(time(chl14)) %in% all_dates]]

# crop to falklands
e <- ext(depth)
sst12 <- crop(sst12, e)
sst13 <- crop(sst13, e)
sst14 <- crop(sst14, e)
chl12 <- crop(chl12, e)
chl13 <- crop(chl13, e)
chl14 <- crop(chl14, e)

# isolate dates for each colony/season
br_13_dates <- br_dates[year(round_date(br_dates, "year")) == 2013]
br_14_dates <- br_dates[year(round_date(br_dates, "year")) == 2014]

cb_13_dates <- cb_dates[year(round_date(cb_dates, "year")) == 2013]
cb_14_dates <- cb_dates[year(round_date(cb_dates, "year")) == 2014]

# collate all sst and chl rasters for each colony/season
br_sst_13 <- c(sst12[[as_date(time(sst12)) %in% br_13_dates]])
br_sst_14 <- c(sst13[[as_date(time(sst13)) %in% br_14_dates]])

cb_sst_13 <- c(sst12[[as_date(time(sst12)) %in% cb_13_dates]])
cb_sst_14 <- c(sst13[[as_date(time(sst13)) %in% cb_14_dates]],
                sst14[[as_date(time(sst14)) %in% cb_14_dates]])

sj_sst_13 <- c(sst12[[as_date(time(sst12)) %in% sj_dates]])

br_chl_13 <- c(chl12[[as_date(time(chl12)) %in% br_13_dates]])
br_chl_14 <- c(chl13[[as_date(time(chl13)) %in% br_14_dates]])

cb_chl_13 <- c(chl12[[as_date(time(chl12)) %in% cb_13_dates]])
cb_chl_14 <- c(chl13[[as_date(time(chl13)) %in% cb_14_dates]],
                chl14[[as_date(time(chl14)) %in% cb_14_dates]])

sj_chl_13 <- c(chl12[[as_date(time(chl12)) %in% sj_dates]])

# compute mean values for each season
br_sst_13 <- app(br_sst_13, mean, na.rm = T)
br_sst_14 <- app(br_sst_14, mean, na.rm = T)

cb_sst_13 <- app(cb_sst_13, mean, na.rm = T)
cb_sst_14 <- app(cb_sst_14, mean, na.rm = T)

sj_sst_13 <- app(sj_sst_13, mean, na.rm = T)

br_chl_13 <- app(br_chl_13, mean, na.rm = T)
br_chl_14 <- app(br_chl_14, mean, na.rm = T)

cb_chl_13 <- app(cb_chl_13, mean, na.rm = T)
cb_chl_14 <- app(cb_chl_14, mean, na.rm = T)

sj_chl_13 <- app(sj_chl_13, mean, na.rm = T)

# resample to match habitat rasters
br_sst_13 <- resample(br_sst_13, depth, method = "bilinear")
br_sst_14 <- resample(br_sst_14, depth, method = "bilinear")
cb_sst_13 <- resample(cb_sst_13, depth, method = "bilinear")
cb_sst_14 <- resample(cb_sst_14, depth, method = "bilinear")
sj_sst_13 <- resample(sj_sst_13, depth, method = "bilinear")
br_chl_13 <- resample(br_chl_13, depth, method = "bilinear")
br_chl_14 <- resample(br_chl_14, depth, method = "bilinear")
cb_chl_13 <- resample(cb_chl_13, depth, method = "bilinear")
cb_chl_14 <- resample(cb_chl_14, depth, method = "bilinear")
sj_chl_13 <- resample(sj_chl_13, depth, method = "bilinear")


# export
writeRaster(br_sst_13, "output/habmod/rasters/br_sst_13.tif", overwrite = T)
writeRaster(br_sst_14, "output/habmod/rasters/br_sst_14.tif", overwrite = T)
writeRaster(cb_sst_13, "output/habmod/rasters/cb_sst_13.tif", overwrite = T)
writeRaster(cb_sst_14, "output/habmod/rasters/cb_sst_14.tif", overwrite = T)
writeRaster(sj_sst_13, "output/habmod/rasters/sj_sst_13.tif", overwrite = T)
writeRaster(br_chl_13, "output/habmod/rasters/br_chl_13.tif", overwrite = T)
writeRaster(br_chl_14, "output/habmod/rasters/br_chl_14.tif", overwrite = T)
writeRaster(cb_chl_13, "output/habmod/rasters/cb_chl_13.tif", overwrite = T)
writeRaster(cb_chl_14, "output/habmod/rasters/cb_chl_14.tif", overwrite = T)
writeRaster(sj_chl_13, "output/habmod/rasters/sj_chl_13.tif", overwrite = T)


#-------------------------------------------------------------------------------
# Extract Values to Dive Intensity Rasters
#-------------------------------------------------------------------------------

# read in dist2cost and colony distance rasters
dist2coast <- rast("output/habmod/rasters/dist2coast.tif")
br_dist <- rast("output/habmod/rasters/dist2BR.tif")
cb_dist <- rast("output/habmod/rasters/dist2CB.tif")
sj_dist <- rast("output/habmod/rasters/dist2SJ.tif")

# read in dive intensity rasters
cb1213 <- rast("data/rasters/CB1213Gdive_densities.tif")
cb1314 <- rast("data/rasters/CB1314Gdive_densities.tif")

br1213 <- rast("data/rasters/BR1213Gdive_densities.tif")
br1314 <- rast("data/rasters/BR1314Gdive_densities.tif")
br1314[is.na(br1314)] <- 0

sj <- rast("data/rasters/SJ1213Idive_densities.tif")
     

# resample rasters to match habitat rasters
cb1213 <- project(cb1213, crs(temp)) %>%
  resample(temp)
br1213 <- project(br1213, crs(temp)) %>%
  resample(temp)
sj <- project(sj, crs(temp)) %>%
  resample(temp)

# create a point within each non-NA cell of the dive rasters
cb13_pts <- as.data.frame(cb1213, xy = T) %>%
  filter(!is.na(CB1213Gdive_densities)) %>%
  vect(geom = c("x", "y"), crs = crs(cb))

cb14_pts <- as.data.frame(cb1314, xy = T) %>%
  filter(!is.na(CB1314Gdive_densities)) %>%
  vect(geom = c("x", "y"), crs = crs(cb))

br13_pts <- as.data.frame(br1213, xy = T) %>%
  filter(!is.na(BR1213Gdive_densities)) %>%
  vect(geom = c("x", "y"), crs = crs(br))

br14_pts <- as.data.frame(br1314, xy = T) %>%
  filter(!is.na(BR1314Gdive_densities)) %>%
  vect(geom = c("x", "y"), crs = crs(br))

sj_pts <- as.data.frame(sj, xy = T) %>%
  filter(!is.na(SJ1213Idive_densities)) %>%
  vect(geom = c("x", "y"), crs = crs(sj))

# extract values from habitat rasters to points
cb13_pts$depth <- extract(depth, cb13_pts, ID = F)
cb13_pts$slope <- extract(slope, cb13_pts, ID = F)
cb13_pts$dist2coast <- extract(dist2coast, cb13_pts, ID = F)
cb13_pts$dist2CB <- extract(cb_dist, cb13_pts, ID = F)
cb13_pts$sst <- extract(cb_sst_13, cb13_pts, ID = F)
cb13_pts$chl <- extract(cb_chl_13, cb13_pts, ID = F)

cb14_pts$depth <- extract(depth, cb14_pts, ID = F)
cb14_pts$slope <- extract(slope, cb14_pts, ID = F)
cb14_pts$dist2coast <- extract(dist2coast, cb14_pts, ID = F)
cb14_pts$dist2CB <- extract(cb_dist, cb14_pts, ID = F)
cb14_pts$sst <- extract(cb_sst_14, cb14_pts, ID = F)
cb14_pts$chl <- extract(cb_chl_14, cb14_pts, ID = F)

br13_pts$depth <- extract(depth, br13_pts, ID = F)
br13_pts$slope <- extract(slope, br13_pts, ID = F)
br13_pts$dist2coast <- extract(dist2coast, br13_pts, ID = F)
br13_pts$dist2BR <- extract(br_dist, br13_pts, ID = F)
br13_pts$sst <- extract(br_sst_13, br13_pts, ID = F)
br13_pts$chl <- extract(br_chl_13, br13_pts, ID = F)

br14_pts$depth <- extract(depth, br14_pts, ID = F)
br14_pts$slope <- extract(slope, br14_pts, ID = F)
br14_pts$dist2coast <- extract(dist2coast, br14_pts, ID = F)
br14_pts$dist2BR <- extract(br_dist, br14_pts, ID = F)
br14_pts$sst <- extract(br_sst_14, br14_pts, ID = F)
br14_pts$chl <- extract(br_chl_14, br14_pts, ID = F)

sj_pts$depth <- extract(depth, sj_pts, ID = F)
sj_pts$slope <- extract(slope, sj_pts, ID = F)
sj_pts$dist2coast <- extract(dist2coast, sj_pts, ID = F)
sj_pts$dist2SJ <- extract(sj_dist, sj_pts, ID = F)
sj_pts$sst <- extract(sj_sst_13, sj_pts, ID = F)
sj_pts$chl <- extract(sj_chl_13, sj_pts, ID = F)

# convert points to dataframes
cb13_pts <- as.data.frame(cb13_pts, geom = "xy")
cb14_pts <- as.data.frame(cb14_pts, geom = "xy")
br13_pts <- as.data.frame(br13_pts, geom = "xy")
br14_pts <- as.data.frame(br14_pts, geom = "xy")
sj_pts <- as.data.frame(sj_pts, geom = "xy")

# rename columns
cb13_pts <- cb13_pts %>%
  rename(dive_intensity = CB1213Gdive_densities,
         dist2col = dist2CB)
cb14_pts <- cb14_pts %>%
  rename(dive_intensity = CB1314Gdive_densities,
         dist2col = dist2CB)
br13_pts <- br13_pts %>%
  rename(dive_intensity = BR1213Gdive_densities,
         dist2col = dist2BR)
br14_pts <- br14_pts %>%
  rename(dive_intensity = BR1314Gdive_densities,
         dist2col = dist2BR)
sj_pts <- sj_pts %>%
  rename(dive_intensity = SJ1213Idive_densities,
         dist2col = dist2SJ)

# export
saveRDS(cb13_pts, "output/habmod/extraction/cowbay13.rds")
saveRDS(cb14_pts, "output/habmod/extraction/cowbay14.rds")
saveRDS(br13_pts, "output/habmod/extraction/bullroad13.rds")
saveRDS(br14_pts, "output/habmod/extraction/bullroad14.rds")
saveRDS(sj_pts, "output/habmod/extraction/steeplejason.rds")