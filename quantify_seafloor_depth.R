#-------------------------------------------------------------------------------
# Seafloor depth characterisation
#-------------------------------------------------------------------------------

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Placement/Writeup/")

library(tidyverse)
library(terra)
library(tidyterra)
library(nlme) 
library(emmeans)

# load data
br <- read.csv("data/tracks/bull_roads.csv") %>%
  filter(breed_stage == "brood_guard")
cb <- read.csv("data/tracks/cow_bay.csv") %>%
  filter(breed_stage == "brood_guard")
sj <- read.csv("data/tracks/steeple_jason.csv")

# combine together
tracks <- bind_rows(br %>% mutate(colony = "bull_roads"),
                        cb %>% mutate(colony = "cow_bay"),
                        sj %>% mutate(colony = "steeple_jason"))

# make spatvector
trax <- vect(tracks, geom = c("longitude", "latitude"), crs = "EPSG:4326")
plot(trax, pch = ".")

# read in depth raster
depth <- rast("data/rasters/depth.tif")

# compute slope
slope <- terrain(depth, v = "slope", unit = "degrees")

# extract depth values at each location
trax$depth <- extract(depth, trax, ID = F)

# extract slope values at each location
trax$slope <- extract(slope, trax, ID = F)

# remove positive depths
trax <- trax %>%
  filter(depth < -1)

# mean depth values by colony
trax %>%
  as.data.frame() %>%
  mutate(date_gmt = as_date(date_gmt, format = "%d/%m/%Y")) %>%
  mutate(season = round_date(date_gmt, "year")) %>%
  group_by(colony) %>%
  summarise(mean_depth = mean(depth, na.rm = T),
            sd_depth = sd(depth, na.rm = T),
            max_depth = min(depth, na.rm = T),
            mean_slope = mean(slope, na.rm = T),
            sd_slope = sd(slope, na.rm = T),
            max_slope = max(slope, na.rm = T),
            n = n())
plot(slope %>% crop(ext(trax)) %>% clamp(-100, 100))
plot(trax, pch = ".", add = T)

# convert to dataframe
tracks <- as.data.frame(trax)

# create the variance structure object
vindent <- varIdent(form=~1|colony)

# fit primary model
model1 <-lme(depth~colony,random=~1|bird_id, method ="REML",data=tracks,
             weights=vindent)

# summary of model
summary(model1)
emmeans(model1, pairwise ~ colony)

# fit slope model
model2 <-lme(slope~colony,random=~1|bird_id, method ="REML",data=tracks,
             weights=vindent)
summary(model2)
emmeans(model2, pairwise ~ colony)

tracks %>%
  group_by(colony) %>%
  summarise(min = min(slope, na.rm = T),
            max = max(slope, na.rm = T),
            mean = mean(slope, na.rm = T),
            sd = sd(slope, na.rm = T),
            uq = quantile(slope, 0.75, na.rm = T),
            lq = quantile(slope, 0.25, na.rm = T),
            n = n())
