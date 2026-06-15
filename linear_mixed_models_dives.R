#-------------------------------------------------------------------------------
# Mixed effects models to analyse gentoo penguin diving and time budget data
# Joshua Wilson (adapted from Norman Ratcliffe & Jonathan Handley)
#-------------------------------------------------------------------------------

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Placement/Writeup/")

#load libraries
library(nlme) #for mixed effects model fitting
library(ape)
library(tidyverse)
library(emmeans)

#read data in csv format, which has following columns: Gender (m or f), Individual (sequential integer for each bird)
#trip (sequential integer for each trip across ALL birds, not WITHIN birds, otherwise models where trip is not nested within
#individual will regard the first trip made by each bird as repeat samples of the same trip, for example),
#botttim (maximim dive depth in m) and Mass (body mass in g)
gentoos <- read.csv("data/FI_Gentoo_Dives2_forLME.csv", header=T)

#-------------------------------------------------------------------------------
# 1. Data Cleaning
#-------------------------------------------------------------------------------

# rename "incubating" to "incubation" within stage
gentoos <- gentoos %>%
  mutate(Stage = ifelse(Stage == "Incubating", "Incubation", Stage))

# convert to factors
gentoos <- gentoos %>%
  mutate(Stage = as.factor(Stage),
         Colony = as.factor(Colony),
         Season = as.factor(Season),
         DayNight = as.factor(DayNight))

# check for missing values
which(is.na(gentoos$botttim))

# subset data by colony, season, and stage (only those that have tracking data)
CB1213 <- subset(gentoos, Colony=="CB" & Season == "2012_13" & Stage == "Guard")
CB1314 <- subset(gentoos, Colony=="CB" & Season == "2013_14" & Stage == "Guard")
BR1213 <- subset(gentoos, Colony=="BR" & Season == "2012_13" & Stage == "Guard")
BR1314 <- subset(gentoos, Colony=="BR" & Season == "2013_14" & Stage == "Incubation" & DepID != "BR13_TDR_2013.14_GentooFI") 
SJ1213 <- subset(gentoos, Colony=="SJ" & Season == "2012_13" & Stage == "Incubation")

# collate into one data frame
gentoos2 <- rbind(CB1213,CB1314,BR1213,BR1314,SJ1213)

# drop unused factor levels
gentoos2 <- droplevels(gentoos2)

# need to correct DiveNo so that it is unique across individual trips
# create a list of all trips of interest
datalist <- unique(gentoos2$TripAll)

# loop over each trip
for(this_trip in datalist){
  
  # isolate data for this trip
  data <- gentoos2 %>%
    filter(TripAll == this_trip) 
  
  # create a new column for dive numbers per trip
  data <- data %>%
    mutate(DiveNo_PerTrip = row_number())
  
  # join to other data
  if(this_trip == datalist[1]){
    final_data <- data
  } else {
    final_data <- bind_rows(final_data, data)
  }
}

# rename dataset to continue
gentoos <- final_data

# clean up
rm(list = setdiff(ls(), "gentoos"))

# create the variance structure object
vindent <- varIdent(form=~1|Colony)

# create the autocorrelation structure object
csT <- corAR1(form=~DiveNo_PerTrip|TripAll)


#-------------------------------------------------------------------------------
# 2. Bottom Time (bottime)
#-------------------------------------------------------------------------------

# check for NAs
sum(is.na(gentoos$botttim)) # lots of NAs - remove
gentoos2 <- gentoos %>%
  filter(!is.na(botttim))

# randomly mutate season for SJ to allow colony * season interaction
gentoos3 <- gentoos2 %>%
  mutate(Season = as.character(Season)) %>%
  mutate(Season = ifelse(Colony == "SJ", 
                         sample(c("2012_13", "2013_14"), size = n(), 
                                replace = TRUE), 
                         Season))

# summarise by colony
ds <- gentoos2 %>%
  group_by(Colony, Season) %>%
  summarise(mean = mean(botttim),
            sd = sd(botttim),
            min = min(botttim),
            max = max(botttim))
ds

ggplot(ds, aes(x = Colony, y = mean, group = Season, col = Season)) +
  geom_point(position=position_dodge(width=0.5), size=3) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(width=0.5)) +
  labs(x = "Colony", y = "Bottom Time (s)", title = "Bottom Time by Colony and Season") +
  theme_classic() +
  theme(legend.position = "top")

# export summary stats
saveRDS(ds, "output/summaries/bottime.rds")

# fit base model 
# model1 <- lme(botttim~Colony,random=~1|TripAll, method ="REML",data=gentoos2, 
#               correlation=csT, weights=vindent,
#               control=lmeControl(maxIter=10, msMaxEval=50))

# save model
#saveRDS(model1, "output/nlme/bottime_model.rds")
model1 <- readRDS("output/nlme/bottime_model.rds")

# estimated means
emmeans1b <- emmeans(model1, pairwise ~ Colony)
emmeans1b

# # refit model C with ML for null comparison - SLOW
# modelC.ml <- lme(botttim~Colony,random=~1|TripAll, method ="ML",data=gentoos2, correlation=csT, weights=vindent,
#                  control=lmeControl(maxIter=10, msMaxEval=50))
# 
# # null model without fixed effects - SLOW
# model.nofixed <- lme(botttim~1,random=~1|TripAll, method ="ML",data=gentoos2, correlation=csT, weights=vindent,
#                      control=lmeControl(maxIter=10, msMaxEval=50))
# 
# # compare models
# anova(modelC.ml,model.nofixed) 

# fit model with interactive season effect
# model1.season <- lme(botttim~Colony*Season,random=~1|TripAll, method ="REML",data=gentoos3, 
#                      correlation=csT, weights=vindent,
#                      control=lmeControl(maxIter=10, msMaxEval=50))


# save model
# saveRDS(model1.season, "output/nlme/bottime_model_season.rds")
model1.season <- readRDS("output/nlme/bottime_model_season.rds")

# estimated means
emmeans1.season <- emmeans(model1.season, pairwise ~ Season | Colony)
emmeans1.season$emmeans
emmeans1.season$contrasts


#-------------------------------------------------------------------------------
# 3. Maximum Depth (maxdep)
#-------------------------------------------------------------------------------

# clean up 
rm(list = setdiff(ls(), c("gentoos", "vindent", "csT")))

# check for NAs
sum(is.na(gentoos$maxdep)) # lots of NAs - remove
gentoos2 <- gentoos %>%
  filter(!is.na(maxdep))

# randomly mutate season for SJ to allow colony * season interaction
gentoos3 <- gentoos2 %>%
  mutate(Season = as.character(Season)) %>%
  mutate(Season = ifelse(Colony == "SJ", 
                         sample(c("2012_13", "2013_14"), size = n(), 
                                replace = TRUE), 
                         Season))

# summarise by colony
ds <- gentoos2 %>%
  group_by(Colony, Season) %>%
  summarise(mean = mean(maxdep),
            sd = sd(maxdep),
            min = min(maxdep),
            max = max(maxdep),
            se = sd/sqrt(n()))
ds

ggplot(ds, aes(x = Colony, y = mean, group = Season, col = Season)) +
  geom_point(position=position_dodge(width=0.5), size=3) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(width=0.5)) +
  labs(x = "Colony", y = "Maximum Depth (m)") +
  theme_classic() +
  theme(legend.position = "top")

# export summary stats
saveRDS(ds, "output/summaries/maxdep.rds")

# fit base model 
# model1 <- lme(maxdep~Colony,random=~1|TripAll, method ="REML",data=gentoos2,
#               correlation=csT, weights=vindent,
#               control=lmeControl(maxIter=10, msMaxEval=50))

# save model
#saveRDS(model1, "output/nlme/maxdep_model.rds")
model1 <- readRDS("output/nlme/maxdep_model.rds")

# estimated means
emmeans1b <- emmeans(model1, pairwise ~ Colony)
emmeans1b

# # refit model C with ML for null comparison - SLOW
# modelC.ml <- lme(maxdep~Colony,random=~1|TripAll, method ="ML",data=gentoos2, correlation=csT, weights=vindent,
#                  control=lmeControl(maxIter=10, msMaxEval=50))
# 
# # null model without fixed effects - SLOW
# model.nofixed <- lme(maxdep~1,random=~1|TripAll, method ="ML",data=gentoos2, correlation=csT, weights=vindent,
#                      control=lmeControl(maxIter=10, msMaxEval=50))
# 
# # compare models
# anova(modelC.ml,model.nofixed) 

# fit model with interactive season effect
# model1.season <- lme(maxdep~Colony*Season,random=~1|TripAll, method ="REML",data=gentoos3,
#                      correlation=csT, weights=vindent,
#                      control=lmeControl(maxIter=10, msMaxEval=50))
# 
# 
# # save model
# saveRDS(model1.season, "output/nlme/maxdep_model_season.rds")
model1.season <- readRDS("output/nlme/maxdep_model_season.rds")

# estimated means
emmeans1.season <- emmeans(model1.season, pairwise ~ Season|Colony)
emmeans1.season$emmeans
emmeans1.season$contrasts


#-------------------------------------------------------------------------------
# 4. Dive Duration (divetime)
#-------------------------------------------------------------------------------

# clean up 
rm(list = setdiff(ls(), c("gentoos", "vindent", "csT")))

# check for NAs
sum(is.na(gentoos$divetim)) # lots of NAs - remove
gentoos2 <- gentoos %>%
  filter(!is.na(divetim))

# randomly mutate season for SJ to allow colony * season interaction
gentoos3 <- gentoos2 %>%
  mutate(Season = as.character(Season)) %>%
  mutate(Season = ifelse(Colony == "SJ", 
                         sample(c("2012_13", "2013_14"), size = n(), 
                                replace = TRUE), 
                         Season))

# summarise by colony
ds <- gentoos2 %>%
  group_by(Colony, Season) %>%
  summarise(mean = mean(divetim),
            sd = sd(divetim), 
            min = min(divetim),
            max = max(divetim))
ds

ggplot(ds, aes(x = Colony, y = mean, group = Season, col = Season)) +
  geom_point(position=position_dodge(width=0.5), size=3) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(width=0.5)) +
  labs(x = "Colony", y = "Maximum Depth (m)") +
  theme_classic() +
  theme(legend.position = "top")

# export summary stats
saveRDS(ds, "output/summaries/divetime.rds")

# fit base model 
# model1 <- lme(divetim~Colony,random=~1|TripAll, method ="REML",data=gentoos2, 
#               correlation=csT, weights=vindent,
#               control=lmeControl(maxIter=10, msMaxEval=50))

# save model
#saveRDS(model1, "output/nlme/divetime_model.rds")
model1 <- readRDS("output/nlme/divetime_model.rds")

# estimated means
emmeans1b <- emmeans(model1, pairwise ~ Colony)
emmeans1b

# # refit model C with ML for null comparison - SLOW
# modelC.ml <- lme(divetim~Colony,random=~1|TripAll, method ="ML",data=gentoos2, correlation=csT, weights=vindent,
#                  control=lmeControl(maxIter=10, msMaxEval=50))
# 
# # null model without fixed effects - SLOW
# model.nofixed <- lme(divetim~1,random=~1|TripAll, method ="ML",data=gentoos2, correlation=csT, weights=vindent,
#                      control=lmeControl(maxIter=10, msMaxEval=50))
# 
# # compare models
# anova(modelC.ml,model.nofixed) 

# fit model with interactive season effect
# model1.season <- lme(divetim~Colony*Season,random=~1|TripAll, method ="REML",data=gentoos3,
#                      correlation=csT, weights=vindent,
#                      control=lmeControl(maxIter=10, msMaxEval=50))
# 
# 
# # save model
# saveRDS(model1.season, "output/nlme/divetime_model_season.rds")
model1.season <- readRDS("output/nlme/divetime_model_season.rds")

# estimated means
emmeans1.season <- emmeans(model1.season, pairwise ~ Season|Colony)
emmeans1.season
