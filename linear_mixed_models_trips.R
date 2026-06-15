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

# read in GPS data
GPS <- read.csv("data/GentooAllGPSTrack.csv")

# round some columns to 2 decimal places
GPS <- GPS %>%
  mutate(TravDist_km = round(TravDist_km, 2),
         MaxDist_km = round(MaxDist_km, 2),
         Duration_hr = round(Duration_hr, 2),
         AvSpeed_kmh = round(AvSpeed_kmh, 2)) 

# rename columns 
GPS <- GPS %>%
  rename(Colony = Colony.L,
         Season = Season.L,
         Stage = Stage.L,
         TravDist = TravDist_km,
         MaxDist = MaxDist_km,
         Duration = Duration_hr,
         AvSpeed = AvSpeed_kmh)

# isolate colonies and seasons of interest
gentoos2 <- GPS %>%
  filter(Colony == "CB" & Season %in% c("2012_13", "2013_14") & Stage == "Guard" |
           Colony == "BR" & Season == "2012_13" & Stage == "Guard" |
           Colony == "BR" & Season == "2013_14" & Stage == "Incubation" & Dep.ID != "BR13_2013-14_GentooFI" |
           Colony == "SJ" & Season == "2012_13" & Stage == "Incubation")

# create the variance structure object
vindent <- varIdent(form=~1|Colony)

# randomly mutate season for SJ to allow colony * season interaction
gentoos3 <- gentoos2 %>%
  mutate(Season = ifelse(Colony == "SJ", 
                         sample(c("2012_13", "2013_14"), size = n(), 
                                replace = TRUE), 
                         Season))

#-------------------------------------------------------------------------------
# 2. Average speed (AvSpeed) 
#-------------------------------------------------------------------------------

# fit primary model
model1 <-lme(AvSpeed~Colony,random=~1|Bird, method ="REML",data=gentoos2,
             weights=vindent)

# summary of model
summary(model1)

# build model with and without fixed effect (use ML method not REML)
model1.ml <-lme(AvSpeed~Colony,random=~1|Bird, method ="ML",data=gentoos2,
                weights=vindent, control = lmeControl(msMaxIter = 1000, msMaxEval = 1000))
model1.nofixed <-lme(AvSpeed~1,random=~1|Bird, method ="ML",data=gentoos2,
                     weights=vindent, control = lmeControl(msMaxIter = 1000, msMaxEval = 1000))

# compare
anova(model1.ml,model1.nofixed) 

# fit model for cow bay and bull roads using season as interactive effect
model1.season <- lme(AvSpeed~Colony*Season,random=~1|Bird, method ="REML",
                     data=gentoos3,
                     weights=vindent)

# compare means for colonies across seasons
emmeans1a <- emmeans(model1.season, pairwise ~ Season | Colony)
emmeans1b <- emmeans(model1, pairwise ~ Colony)
emmeans1a
emmeans1b

# statistical differences
emmeans1_p <- emmeans1a$contrasts

# identify groups with no differences
diffs1 <- emmeans1_p %>%
  as.data.frame() %>%
  #filter(p.value > 0.05) %>%
  select(contrast, p.value, t.ratio) %>%
  mutate(colony_season1 = str_split(contrast, " - ", simplify = TRUE)[, 1],
         colony_season2 = str_split(contrast, " - ", simplify = TRUE)[, 2]) %>%
  select(colony_season1, colony_season2, p.value, t.ratio) %>%
  mutate(p.value = round(p.value, 4)) 
diffs1

# create summary table
emmeans1a <- emmeans1a$emmeans %>%
  as.data.frame() %>%
  filter(Colony != "SJ")
emmeans1b <- emmeans1b$emmeans %>%
  as.data.frame() %>%
  filter(Colony == "SJ") %>%
  mutate(Season = as.factor("2012_13"))
summary_table <- rbind(emmeans1a, emmeans1b) %>%
  rename(mean = emmean) %>%
  as.data.frame()

summary_table <- gentoos2 %>% 
  group_by(Colony, Season) %>%
  summarise(mean = mean(AvSpeed),
            sd = sd(AvSpeed),
            min = min(AvSpeed),
            max = max(AvSpeed))
summary_table

ggplot(summary_table %>% mutate(group = paste(Colony, Season)), 
       aes(x = group, y = mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
  labs(x = "Colony", y = "Average Speed (km/h)", title = "Average Speed by Colony and Season") +
  theme_minimal() +
  theme(legend.title = element_blank())

# export summary stats
saveRDS(summary_table, "output/summaries/avspeed.rds")


#-------------------------------------------------------------------------------
# 3. Maximum distance (MaxDist)
#-------------------------------------------------------------------------------

# fit primary model
model1 <-lme(MaxDist~Colony,random=~1|Bird, method ="REML",data=gentoos2,
             weights=vindent)

# summary of model
summary(model1)

# build model with and without fixed effect (use ML method not REML)
model1.ml <-lme(MaxDist~Colony,random=~1|Bird, method ="ML",data=gentoos2,
                weights=vindent, control = lmeControl(msMaxIter = 1000, msMaxEval = 1000))
model1.nofixed <-lme(MaxDist~1,random=~1|Bird, method ="ML",data=gentoos2,
                     weights=vindent, control = lmeControl(msMaxIter = 1000, msMaxEval = 1000))

# compare
anova(model1.ml,model1.nofixed) 

# fit model for cow bay and bull roads using season as interactive effect
model1.season <- lme(MaxDist~Colony*Season,random=~1|Bird, method ="REML",
                     data=gentoos3,
                     weights=vindent)

# compare means for colonies across seasons
emmeans1a <- emmeans(model1.season, pairwise ~ Season | Colony)
emmeans1b <- emmeans(model1, pairwise ~ Colony)
emmeans1a
emmeans1b

# statistical differences
emmeans1_p <- emmeans1a$contrasts

# identify groups with no differences
diffs1 <- emmeans1_p %>%
  as.data.frame() %>%
  #filter(p.value > 0.05) %>%
  select(contrast, p.value, t.ratio) %>%
  mutate(colony_season1 = str_split(contrast, " - ", simplify = TRUE)[, 1],
         colony_season2 = str_split(contrast, " - ", simplify = TRUE)[, 2]) %>%
  select(colony_season1, colony_season2, p.value, t.ratio) %>%
  mutate(p.value = round(p.value, 4)) 
diffs1

# create summary table
emmeans1a <- emmeans1a$emmeans %>%
  as.data.frame() %>%
  filter(Colony != "SJ")
emmeans1b <- emmeans1b$emmeans %>%
  as.data.frame() %>%
  filter(Colony == "SJ") %>%
  mutate(Season = as.factor("2012_13"))
summary_table <- rbind(emmeans1a, emmeans1b) %>%
  rename(mean = emmean) %>%
  as.data.frame()

summary_table <- gentoos2 %>% 
  group_by(Colony, Season) %>%
  summarise(mean = mean(MaxDist),
            sd = sd(MaxDist),
            min = min(MaxDist),
            max = max(MaxDist))
summary_table

ggplot(summary_table %>% mutate(group = paste(Colony, Season)), 
       aes(x = group, y = mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
  labs(x = "Colony", y = "Maximum Distance (km)", title = "Maximum Distance by Colony and Season") +
  theme_minimal() +
  theme(legend.title = element_blank())

# export summary stats
saveRDS(summary_table, "output/summaries/maxdist.rds")


#-------------------------------------------------------------------------------
# 4. Trip duration (Duration)
#-------------------------------------------------------------------------------

# fit primary model
model1 <-lme(Duration~Colony,random=~1|Bird, method ="REML",data=gentoos2,
             weights=vindent)

# summary of model
summary(model1)

# build model with and without fixed effect (use ML method not REML)
model1.ml <-lme(Duration~Colony,random=~1|Bird, method ="ML",data=gentoos2,
                weights=vindent, control = lmeControl(msMaxIter = 1000, msMaxEval = 1000))
model1.nofixed <-lme(Duration~1,random=~1|Bird, method ="ML",data=gentoos2,
                     weights=vindent, control = lmeControl(msMaxIter = 1000, msMaxEval = 1000))

# compare
anova(model1.ml,model1.nofixed)

# fit model for cow bay and bull roads using season as interactive effect
model1.season <- lme(Duration~Colony*Season,random=~1|Bird, method ="REML",
                     data=gentoos3,
                     weights=vindent)

# compare means for colonies across seasons
emmeans1a <- emmeans(model1.season, pairwise ~ Season | Colony)
emmeans1b <- emmeans(model1, pairwise ~ Colony)
emmeans1a
emmeans1b

# statistical differences
emmeans1_p <- emmeans1a$contrasts

# identify groups with no differences
diffs1 <- emmeans1_p %>%
  as.data.frame() %>%
  #filter(p.value > 0.05) %>%
  select(contrast, p.value, t.ratio) %>%
  mutate(colony_season1 = str_split(contrast, " - ", simplify = TRUE)[, 1],
         colony_season2 = str_split(contrast, " - ", simplify = TRUE)[, 2]) %>%
  select(colony_season1, colony_season2, p.value, t.ratio) %>%
  mutate(p.value = round(p.value, 4)) 
diffs1

# create summary table
emmeans1a <- emmeans1a$emmeans %>%
  as.data.frame() %>%
  filter(Colony != "SJ")
emmeans1b <- emmeans1b$emmeans %>%
  as.data.frame() %>%
  filter(Colony == "SJ") %>%
  mutate(Season = as.factor("2012_13"))
summary_table <- rbind(emmeans1a, emmeans1b) %>%
  rename(mean = emmean) %>%
  as.data.frame()

summary_table <- gentoos2 %>% 
  group_by(Colony, Season) %>%
  summarise(mean = mean(Duration),
            sd = sd(Duration),
            min = min(Duration),
            max = max(Duration))
summary_table

ggplot(summary_table %>% mutate(group = paste(Colony, Season)), 
       aes(x = group, y = mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
  labs(x = "Colony", y = "Trip Duration (hours)", title = "Trip Duration by Colony and Season") +
  theme_minimal() +
  theme(legend.title = element_blank())

# export summary stats
saveRDS(summary_table, "output/summaries/duration.rds")


#-------------------------------------------------------------------------------
# 5. Foraging trip distance (TravDist)
#-------------------------------------------------------------------------------

# fit primary model
model1 <-lme(TravDist~Colony,random=~1|Bird, method ="REML",data=gentoos2,
             weights=vindent)

# summary of model
summary(model1)

# build model with and without fixed effect (use ML method not REML)
model1.ml <-lme(TravDist~Colony,random=~1|Bird, method ="ML",data=gentoos2,
                weights=vindent, control = lmeControl(msMaxIter = 1000, msMaxEval = 1000))
model1.nofixed <-lme(TravDist~1,random=~1|Bird, method ="ML",data=gentoos2,
                     weights=vindent, control = lmeControl(msMaxIter = 1000, msMaxEval = 1000))

# compare
anova(model1.ml,model1.nofixed) 

# fit model for cow bay and bull roads using season as interactive effect
model1.season <- lme(TravDist~Colony*Season,random=~1|Bird, method ="REML",
                     data=gentoos3,
                     weights=vindent)

# compare means for colonies across seasons
emmeans1a <- emmeans(model1.season, pairwise ~ Season | Colony)
emmeans1b <- emmeans(model1, pairwise ~ Colony)
emmeans1a
emmeans1b

# statistical differences
emmeans1_p <- emmeans1a$contrasts

# identify groups with no differences
diffs1 <- emmeans1_p %>%
  as.data.frame() %>%
  #filter(p.value > 0.05) %>%
  select(contrast, p.value, t.ratio) %>%
  mutate(colony_season1 = str_split(contrast, " - ", simplify = TRUE)[, 1],
         colony_season2 = str_split(contrast, " - ", simplify = TRUE)[, 2]) %>%
  select(colony_season1, colony_season2, p.value, t.ratio) %>%
  mutate(p.value = round(p.value, 4)) 
diffs1

# create summary table
emmeans1a <- emmeans1a$emmeans %>%
  as.data.frame() %>%
  filter(Colony != "SJ")
emmeans1b <- emmeans1b$emmeans %>%
  as.data.frame() %>%
  filter(Colony == "SJ") %>%
  mutate(Season = as.factor("2012_13"))
summary_table <- rbind(emmeans1a, emmeans1b) %>%
  rename(mean = emmean) %>%
  as.data.frame()

summary_table <- gentoos2 %>% 
  group_by(Colony, Season) %>%
  summarise(mean = mean(TravDist),
            sd = sd(TravDist),
            min = min(TravDist),
            max = max(TravDist))
summary_table

ggplot(summary_table %>% mutate(group = paste(Colony, Season)), 
       aes(x = group, y = mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
  labs(x = "Colony", y = "Foraging Trip Distance (km)", title = "Foraging Trip Distance by Colony and Season") +
  theme_minimal() +
  theme(legend.title = element_blank())

# export summary stats
saveRDS(summary_table, "output/summaries/travdist.rds")
