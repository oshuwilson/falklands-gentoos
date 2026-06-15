#-------------------------------------------------------------------------------
# Plot GPS Tracks over Dynamic Variables to Explain Exclusion
# Joshua Wilson
#-------------------------------------------------------------------------------

rm(list=ls())
setwd("E:/")

library(tidyverse)
library(terra)
library(tidyterra)
library(sf)

# read in SST and Chlorophyll
sst <- rast("Satellite_Data/monthly/sst/sst.nc")
chl12 <- rast("Satellite_Data/daily/chl/resampled/chl_2012_resampled.nc")
chl13 <- rast("Satellite_Data/daily/chl/resampled/chl_2013_resampled.nc")

# limit sst to December 2012 and December 2013
sst12 <- sst[[time(sst) == as_date("2012-12-01")]]
sst13 <- sst[[time(sst) == as_date("2013-12-01")]]
rm(sst)

# limit chl to December of the respective years
chl12 <- chl12[[month(time(chl12)) == 12]]
chl13 <- chl13[[month(time(chl13)) == 12]]

# falklands extent
e <- ext(-62, -56, -53, -50)

# crop to falklands
sst12 <- crop(sst12, e)
sst13 <- crop(sst13, e)

chl12 <- crop(chl12, e)
chl13 <- crop(chl13, e)

# compute average chlorophyll over December
chl12 <- app(chl12, mean, na.rm = T)
chl13 <- app(chl13, mean, na.rm = T)

# read in Falklands shapefile
fk <- read_sf("~/OneDrive - University of Southampton/Documents/Placement/Writeup/data/falklands_shapefile/All_Falkls.shp") %>%
  vect()

# create bounding boxes for each study site
br <- ext(-59.5, -59, -52.5, -52.2)
cb <- ext(-58, -57.2, -51.5, -51)
sj <- ext(-61.8, -61, -51.2, -50.8)

# convert to spatvectors
br <- vect(br)
cb <- vect(cb)
sj <- vect(sj)

# add CRS
crs(br) = "EPSG:4326"
crs(cb) = "EPSG:4326"
crs(sj) = "EPSG:4326"

# plot overlaid on SST/CHL
p1 <- ggplot() +
  geom_spatraster(data = chl13) +
  geom_spatvector(data = fk, col = NA, fill = "grey70") +
  geom_spatvector(data = br, fill = NA, col = "white") +
  geom_spatvector(data = cb, fill = NA, col = "white") +
  geom_spatvector(data = sj, fill = NA, col = "white") +
  scale_fill_viridis_c(na.value = "grey70", trans = "log", 
                       name = "Chlorophyll Concentration\n2013 (mg/m3)") +
  theme_minimal() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(panel.border = element_rect(colour = "black", fill = NA))
p1

p2 <- ggplot() +
  geom_spatraster(data = chl12) +
  geom_spatvector(data = fk, col = NA, fill = "grey70") +
  geom_spatvector(data = br, fill = NA, col = "white") +
  geom_spatvector(data = cb, fill = NA, col = "white") +
  geom_spatvector(data = sj, fill = NA, col = "white") +
  scale_fill_viridis_c(na.value = "grey70", trans = "log", 
                       name = "Chlorophyll Concentration\n2012 (mg/m3)") +
  theme_minimal() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(panel.border = element_rect(colour = "black", fill = NA))
p2

p3 <- ggplot() +
  geom_spatraster(data = sst13) +
  geom_spatvector(data = fk, col = NA, fill = "grey70") +
  geom_spatvector(data = br, fill = NA, col = "white") +
  geom_spatvector(data = cb, fill = NA, col = "white") +
  geom_spatvector(data = sj, fill = NA, col = "white") +
  scale_fill_viridis_c(na.value = "grey70", name = "Sea Surface Temperature\n2013 (°C)") +
  theme_minimal() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(panel.border = element_rect(colour = "black", fill = NA))
p3

p4 <- ggplot() +
  geom_spatraster(data = sst12) +
  geom_spatvector(data = fk, col = NA, fill = "grey70") +
  geom_spatvector(data = br, fill = NA, col = "white") +
  geom_spatvector(data = cb, fill = NA, col = "white") +
  geom_spatvector(data = sj, fill = NA, col = "white") +
  scale_fill_viridis_c(na.value = "grey70", name = "Sea Surface Temperature\n2012 (°C)") +
  theme_minimal() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(panel.border = element_rect(colour = "black", fill = NA))
p4

# extract values to boxes
brsd <- data.frame(
  sst12 = extract(sst12, br, fun = sd, na.rm = T, ID = F) %>% pull(1),
  sst13 = extract(sst13, br, fun = sd, na.rm = T, ID = F) %>% pull(1),
  chl12 = extract(chl12, br, fun = sd, na.rm = T, ID = F) %>% pull(1),
  chl13 = extract(chl13, br, fun = sd, na.rm = T, ID = F) %>% pull(1))
brsd

brmean <- data.frame(
  sst12 = extract(sst12, br, fun = mean, na.rm = T, ID = F) %>% pull(1),
  sst13 = extract(sst13, br, fun = mean, na.rm = T, ID = F) %>% pull(1),
  chl12 = extract(chl12, br, fun = mean, na.rm = T, ID = F) %>% pull(1),
  chl13 = extract(chl13, br, fun = mean, na.rm = T, ID = F) %>% pull(1))
brmean

cbsd <- data.frame(
  sst12 = extract(sst12, cb, fun = sd, na.rm = T, ID = F) %>% pull(1),
  sst13 = extract(sst13, cb, fun = sd, na.rm = T, ID = F) %>% pull(1),
  chl12 = extract(chl12, cb, fun = sd, na.rm = T, ID = F) %>% pull(1),
  chl13 = extract(chl13, cb, fun = sd, na.rm = T, ID = F) %>% pull(1))
cbsd

cbmean <- data.frame(
  sst12 = extract(sst12, cb, fun = mean, na.rm = T, ID = F) %>% pull(1),
  sst13 = extract(sst13, cb, fun = mean, na.rm = T, ID = F) %>% pull(1),
  chl12 = extract(chl12, cb, fun = mean, na.rm = T, ID = F) %>% pull(1),
  chl13 = extract(chl13, cb, fun = mean, na.rm = T, ID = F) %>% pull(1))
cbmean

sjsd <- data.frame(
  sst12 = extract(sst12, sj, fun = sd, na.rm = T, ID = F) %>% pull(1),
  chl12 = extract(chl12, sj, fun = sd, na.rm = T, ID = F) %>% pull(1))
sjsd

sjmean <- data.frame(
  sst12 = extract(sst12, sj, fun = mean, na.rm = T, ID = F) %>% pull(1),
  chl12 = extract(chl12, sj, fun = mean, na.rm = T, ID = F) %>% pull(1))
sjmean

# extract values for entire rasters
global <- data.frame(
  sst12sd = global(sst12, sd, na.rm = T) %>% pull(1),
  sst13sd = global(sst13, sd, na.rm = T) %>% pull(1),
  chl12sd = global(chl12, sd, na.rm = T) %>% pull(1),
  chl13sd = global(chl13, sd, na.rm = T) %>% pull(1))
global

globalmean <- data.frame(
  sst12mean = global(sst12, mean, na.rm = T) %>% pull(1),
  sst13mean = global(sst13, mean, na.rm = T) %>% pull(1),
  chl12mean = global(chl12, mean, na.rm = T) %>% pull(1),
  chl13mean = global(chl13, mean, na.rm = T) %>% pull(1))
globalmean

# number of cells in each box
ncell(crop(sst12, br))
ncell(crop(sst12, cb))
ncell(crop(sst12, sj))

# plot all together
env_vars <- cowplot::plot_grid(p2, p1, p4, p3, ncol = 2)
env_vars + ggview::canvas(width = 14, height = 10)

# export
ggsave("~/OneDrive - University of Southampton/Documents/Placement/Writeup/output/imagery/supplementary/dynamic_vars.png",
       env_vars, width = 14, height = 10)


# extract values for rasters
sst12_vals <- values(sst12, na.rm = T)
sst13_vals <- values(sst13, na.rm = T)
chl12_vals <- values(chl12, na.rm = T)
chl13_vals <- values(chl13, na.rm = T)

# test for interannual differences 
wilcox.test(sst12_vals, sst13_vals)
wilcox.test(chl12_vals, chl13_vals)

# repeat for each bounding box
sst12_br_vals <- values(crop(sst12, br), na.rm = T)
sst13_br_vals <- values(crop(sst13, br), na.rm = T)
chl12_br_vals <- values(crop(chl12, br), na.rm = T)
chl13_br_vals <- values(crop(chl13, br), na.rm = T)
wilcox.test(sst12_br_vals, sst13_br_vals)
wilcox.test(chl12_br_vals, chl13_br_vals)

sst12_cb_vals <- values(crop(sst12, cb), na.rm = T)
sst13_cb_vals <- values(crop(sst13, cb), na.rm = T)
chl12_cb_vals <- values(crop(chl12, cb), na.rm = T)
chl13_cb_vals <- values(crop(chl13, cb), na.rm = T)
wilcox.test(sst12_cb_vals, sst13_cb_vals)
wilcox.test(chl12_cb_vals, chl13_cb_vals)
