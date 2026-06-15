#-------------------------------------------------------------------------------
# Visualise the Outputs from Boosted Regression Trees
# Joshua Wilson
#-------------------------------------------------------------------------------

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Placement/Writeup/")

# load libraries 
{
  library(tidyverse)
  library(terra) 
  library(tidyterra)
  library(sf)
}

#-------------------------------------------------------------------------------
# 1. Plot Raster Predictions
#-------------------------------------------------------------------------------

# read in Falklands shapefile
fk <- read_sf("data/falklands_shapefile/All_Falkls.shp") %>%
  vect()

# read in Falklands bathymetry
bathy <- readRDS("data/falklands_bathy.RDS") %>%
  crop(ext(-62, -56, -53, -50))

# read in zero models
zero_br <- rast("output/habmod/predictions/brt_zero_pred_bullroad.tif")
zero_cb <- rast("output/habmod/predictions/brt_zero_pred_cowbay.tif")
zero_sj <- rast("output/habmod/predictions/brt_zero_pred_steeplejason.tif")

# combine zero models into one raster
zero <- merge(zero_br, zero_cb, zero_sj)

# if 0 < 0.2, floor to 0
zero <- ifel(zero < 0.2, 0, zero)

# make NA values 0
zero[is.na(zero)] <- 0
plot(zero)

# read in intensity models
intensity_br <- rast("output/habmod/predictions/brt_intensity_pred_bullroad.tif")
intensity_cb <- rast("output/habmod/predictions/brt_intensity_pred_cowbay.tif")
intensity_sj <- rast("output/habmod/predictions/brt_intensity_pred_steeplejason.tif")

# combine intensity models into one raster
intensity <- merge(intensity_br, intensity_cb, intensity_sj)

# make NA values 0 
intensity[is.na(intensity)] <- 0

# revalue negative values to 0
intensity[intensity < 0] <- 0
plot(intensity)

# combine intensity and zero model
combined <- zero * intensity
plot(combined)

# # read in combined models
# combined_br <- rast("output/habmod/predictions/brt_final_pred_bullroad.tif")
# combined_cb <- rast("output/habmod/predictions/brt_final_pred_cowbay.tif")
# combined_sj <- rast("output/habmod/predictions/brt_final_pred_steeplejason.tif")
# 
# # combine combined models into one raster
# combined <- merge(combined_br, combined_cb, combined_sj)

# make NA values 0 
combined[is.na(combined)] <- 0
plot(combined)

# revalue low occurrence values to NA
zero[zero < 0.2] <- NA
plot(zero)

# extents for the three regions
br_ext <- ext(-60, -58.5, -53, -52)
cb_ext <- ext(-59, -56.5, -52, -50.5)
sj_ext <- ext(-62, -60, -51.5, -50.5)

# crop combined rasters to the three regions
bullroad <- crop(combined, br_ext)
cowbay <- crop(combined, cb_ext)
steeple <- crop(combined, sj_ext)
plot(bullroad)
plot(cowbay)
plot(steeple)

# scale regional predictions to 0-1
bullroad <- (bullroad - min(values(bullroad), na.rm = T)) / (max(values(bullroad), na.rm = T) - min(values(bullroad), na.rm = T))
cowbay <- (cowbay - min(values(cowbay), na.rm = T)) / (max(values(cowbay), na.rm = T) - min(values(cowbay), na.rm = T))
steeple <- (steeple - min(values(steeple), na.rm = T)) / (max(values(steeple), na.rm = T) - min(values(steeple), na.rm = T))

# plot - whole of falklands
p1 <- ggplot() +
  geom_spatraster(data = bullroad) +
  geom_spatraster(data = cowbay) +
  geom_spatraster(data = steeple) +
  geom_spatvector(data = bathy, col = "grey70") +
  geom_spatvector(data = fk, fill = "grey70", col = "grey70") +
  scale_fill_viridis_c(trans = "sqrt", na.value = "transparent", 
                       name = "Predicted Dive\nIntensity") +
  theme_bw() +
  scale_y_continuous(limits = c(-53, -50.5), expand = c(0,0)) +
  scale_x_continuous(limits = c(-62, -56.5), expand = c(0,0)) +
  theme(panel.background = element_rect(fill = "#440154FF"),
        panel.grid = element_blank())

# plot - bull roads
p2 <- ggplot() +
  geom_spatraster(data = bullroad) +
  geom_spatvector(data = bathy, col = "grey70") +
  geom_spatvector(data = fk, fill = "grey70", col = "grey70") +
  scale_fill_viridis_c(trans = "sqrt", na.value = "transparent", guide = "none") +
  theme_bw() +
  scale_y_continuous(limits = c(-53, -52), expand = c(0,0)) +
  scale_x_continuous(limits = c(-60, -58.5), expand = c(0,0)) +
  theme(panel.background = element_rect(fill = "#440154FF"),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()) 

# plot - cow bay
p3 <- ggplot() +
  geom_spatraster(data = cowbay) +
  geom_spatvector(data = bathy, col = "grey70") +
  geom_spatvector(data = fk, fill = "grey70", col = "grey70") +  
  scale_fill_viridis_c(trans = "sqrt", na.value = "transparent", guide = "none") +
  theme_bw() +
  scale_y_continuous(limits = c(-52, -50.6), expand = c(0,0)) +
  scale_x_continuous(limits = c(-58.8, -56.8), expand = c(0,0)) +
  theme(panel.background = element_rect(fill = "#440154FF"),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())

# plot - steeple
p4 <- ggplot() +
  geom_spatraster(data = steeple) +
  geom_spatvector(data = bathy, col = "grey70") +
  geom_spatvector(data = fk, fill = "grey70", col = "grey70") +
  scale_fill_viridis_c(trans = "sqrt", na.value = "transparent", guide = "none") +
  theme_bw() +
  scale_y_continuous(limits = c(-51.5, -50.6), expand = c(0,0)) +
  scale_x_continuous(limits = c(-61.8, -60.4), expand = c(0,0),
                     breaks = seq(-62, -60, 0.5)) +
  theme(panel.background = element_rect(fill = "#440154FF"),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())

# get legend from p1
legend <- cowplot::get_legend(p1)
plot(legend)

# remove legend from p1
p1 <- p1 + theme(legend.position = "none")

# export plots for compilation
p1 + ggview::canvas(width = 8, height= 6)
ggsave("output/imagery/brt_raster/fktest2.png", p1,
       width = 180, height = 135, units = "mm", dpi = 300)
ggsave("output/imagery/brt_raster/brtest2.png", p2,
       width = 130, height = 100, units = "mm", dpi = 300)
ggsave("output/imagery/brt_raster/cbtest2.png", p3,
       width = 8, height = 6, units = "in", dpi = 300)
ggsave("output/imagery/brt_raster/sjtest2.png", p4,
       width = 120, height = 90, units = "mm", dpi = 300)

ggsave("output/imagery/brt_raster/legend.png", legend,
       width = 2, height = 3, units = "in", dpi = 300)


#-------------------------------------------------------------------------------
# Repeat for intensity model only
#-------------------------------------------------------------------------------

# crop intensity rasters to the three regions
bullroad <- crop(intensity, br_ext)
cowbay <- crop(intensity, cb_ext)
steeple <- crop(intensity, sj_ext)

# scale regional predictions to 0-1
bullroad <- (bullroad - min(values(bullroad), na.rm = T)) / (max(values(bullroad), na.rm = T) - min(values(bullroad), na.rm = T))
cowbay <- (cowbay - min(values(cowbay), na.rm = T)) / (max(values(cowbay), na.rm = T) - min(values(cowbay), na.rm = T))
steeple <- (steeple - min(values(steeple), na.rm = T)) / (max(values(steeple), na.rm = T) - min(values(steeple), na.rm = T))

# plot - whole of falklands
p1 <- ggplot() +
  geom_spatraster(data = bullroad) +
  geom_spatraster(data = cowbay) +
  geom_spatraster(data = steeple) +
  geom_spatvector(data = bathy, col = "grey70") +
  geom_spatvector(data = fk, fill = "grey70", col = "grey70") +
  scale_fill_viridis_c(trans = "sqrt", na.value = "transparent", 
                       name = "Predicted Dive Intensity") +
  theme_bw() +
  scale_y_continuous(limits = c(-53, -50.5), expand = c(0,0)) +
  scale_x_continuous(limits = c(-62, -56.5), expand = c(0,0)) +
  theme(panel.background = element_rect(fill = "#440154FF"),
        panel.grid = element_blank())
p1 + ggview::canvas(width = 12, height = 10)

# export plotp1 + ggview::canvas(width = 8, height= 6)
ggsave("output/imagery/supplementary/intensity_raster.png", p1,
       width = 12, height = 10, dpi = 300)


#-------------------------------------------------------------------------------
# Repeat for zero model only
#-------------------------------------------------------------------------------

# crop zero rasters to the three regions
bullroad <- crop(zero, br_ext)
cowbay <- crop(zero, cb_ext)
steeple <- crop(zero, sj_ext)

# scale regional predictions to 0-1
bullroad <- (bullroad - min(values(bullroad), na.rm = T)) / (max(values(bullroad), na.rm = T) - min(values(bullroad), na.rm = T))
cowbay <- (cowbay - min(values(cowbay), na.rm = T)) / (max(values(cowbay), na.rm = T) - min(values(cowbay), na.rm = T))
steeple <- (steeple - min(values(steeple), na.rm = T)) / (max(values(steeple), na.rm = T) - min(values(steeple), na.rm = T))

# plot - whole of falklands
p1 <- ggplot() +
  geom_spatraster(data = bullroad) +
  geom_spatraster(data = cowbay) +
  geom_spatraster(data = steeple) +
  geom_spatvector(data = bathy, col = "grey70") +
  geom_spatvector(data = fk, fill = "grey70", col = "grey70") +
  scale_fill_viridis_c(na.value = "transparent", 
                       name = "Probability of Presence") +
  theme_bw() +
  scale_y_continuous(limits = c(-53, -50.5), expand = c(0,0)) +
  scale_x_continuous(limits = c(-62, -56.5), expand = c(0,0)) +
  theme(panel.background = element_rect(fill = "#440154FF"),
        panel.grid = element_blank())
p1 + ggview::canvas(width = 12, height = 10)

# export plotp1 + ggview::canvas(width = 8, height= 6)
ggsave("output/imagery/supplementary/zero_raster.png", p1,
       width = 12, height = 10, dpi = 300)
