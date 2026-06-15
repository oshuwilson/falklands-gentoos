#-------------------------------------------------------------------------------
# Mixed effects models to analyse gentoo penguin foraging trip data
# Joshua Wilson (adapted from Jonathan Handley)
#-------------------------------------------------------------------------------

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Placement/Writeup/")

#load libraries
library(nlme) 
library(ape)
library(tidyverse)
library(MuMIn)
library(emmeans)

#-------------------------------------------------------------------------------
# 1. Data Cleaning
#-------------------------------------------------------------------------------

# read in data
data <- readRDS("data/PercentBenthicDives.RDS")

# isolate colonies and seasons of interest
gentoos2 <- data %>%
  filter(Colony == "CB" & Season %in% c("2012_13", "2013_14") & Stage == "Guard" |
           Colony == "BR" & Season == "2012_13" & Stage == "Guard" |
           Colony == "BR" & Season == "2013_14" & Stage == "Incubation" & DepID != "BR13_2013-14_GentooFI" |
           Colony == "SJ" & Season == "2012_13" & Stage == "Incubation")

# summary stats
ds <- gentoos2 %>%
  group_by(Colony, Season) %>%
  summarise(
    mean = mean(PercentBenthic),
    sd = sd(PercentBenthic),
    min = min(PercentBenthic),
    max = max(PercentBenthic)
  ) 
ds

# export summary stats
saveRDS(ds, "output/summaries/prop_benthic.RDS")

#-------------------------------------------------------------------------------
# 2. Mixed Effects Model Fitting
#-------------------------------------------------------------------------------

# create the variance structure object
vindent <- varIdent(form=~1|Colony)

# fixed effects model
model1 <- lme(PercentBenthic ~ Colony, random = ~1|DepID, method = "REML", 
              weights = vindent, data = gentoos2)
summary(model1)
emmeans(model1, pairwise ~ Colony)

# build model with and without fixed effect (both ML method)
model1.ml <-lme(PercentBenthic~Colony,random=~1|DepID, method ="ML",data=gentoos2,
                weights=vindent, control = lmeControl(msMaxIter = 1000, msMaxEval = 1000))
model1.nofixed <-lme(PercentBenthic~1,random=~1|DepID, method ="ML",data=gentoos2,
                     weights=vindent, control = lmeControl(msMaxIter = 1000, msMaxEval = 1000))

# compare 
anova(model1.ml, model1.nofixed) 

# randomly mutate season for steeple jason to allow interaction
gentoos2 <- gentoos2 %>%
  mutate(Season = as.character(Season)) %>%
  mutate(Season = ifelse(Colony == "SJ", 
                         sample(c("2012_13", "2013_14"), size = n(), 
                                replace = TRUE), 
                         Season))

# build model with interaction
model1.season <- lme(PercentBenthic~Colony*Season,random=~1|DepID, method ="ML",data=gentoos2,
                 weights=vindent, control = lmeControl(msMaxIter = 1000, msMaxEval = 1000))
summary(model1.season)

emmeans.season <- emmeans(model1.season, pairwise ~ Season | Colony)
emmeans.season
emmeans(model1, pairwise ~ Colony)
