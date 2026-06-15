#-------------------------------------------------------------------------------
# Bootstrap Tests
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


#-------------------------------------------------------------------------------
# 2. Bootstrapping 
#-------------------------------------------------------------------------------

# list all bull roads and cow bay individuals
BR_individuals <- gentoos2 %>%
  filter(Colony == "BR") %>%
  pull(Dep.ID) %>%
  unique()

CB_individuals <- gentoos2 %>%
  filter(Colony == "CB") %>%
  pull(Dep.ID) %>%
  unique()

# this iteration
for(i in 1:1000){
  
  if (i %% 100 == 0) {
    print(paste("Iteration:", i))
  }
  
  # set seed
  set.seed(i)
  
  # pick a random 8 individuals from Cow Bay/Bull Roads
  CB_sample <- sample(CB_individuals, 8)
  BR_sample <- sample(BR_individuals, 8)
  
  # filter data to these individuals
  boot_data <- gentoos2 %>%
    filter(Colony %in% c("CB", "BR") & Dep.ID %in% c(CB_sample, BR_sample))
  
  # get mean and sd of trip characteristics
  boot_summary <- boot_data %>%
    group_by(Colony) %>%
    summarise(mean_dist = mean(TravDist),
              sd_dist = sd(TravDist),
              mean_dur = mean(Duration),
              sd_dur = sd(Duration),
              mean_speed = mean(AvSpeed),
              sd_speed = sd(AvSpeed),
              mean_max = mean(MaxDist),
              sd_max = sd(MaxDist))
  
  # add iteration
  boot_summary <- boot_summary %>%
    mutate(iteration = i)
  
  # bind together
  if (i == 1) {
    boot_results <- boot_summary
  } else {
    boot_results <- bind_rows(boot_results, boot_summary)
  }
  
}

# get overall means for gentoos2 data
overall_means <- gentoos2 %>%
  filter(Colony != "SJ") %>%
  group_by(Colony) %>%
  summarise(mean_dist = mean(TravDist),
            mean_dur = mean(Duration),
            mean_speed = mean(AvSpeed),
            mean_max = mean(MaxDist),
            sd_dist = sd(TravDist),
            sd_dur = sd(Duration),
            sd_speed = sd(AvSpeed),
            sd_max = sd(MaxDist),
            se_dist = sd(TravDist)/sqrt(n()),
            se_dur = sd(Duration)/sqrt(n()),
            se_speed = sd(AvSpeed)/sqrt(n()),
            se_max = sd(MaxDist)/sqrt(n()))

# steeple jason means
SJ_means <- gentoos2 %>%
  filter(Colony == "SJ") %>%
  group_by(Colony) %>%
  summarise(mean_dist = mean(TravDist),
            mean_dur = mean(Duration),
            mean_speed = mean(AvSpeed),
            mean_max = mean(MaxDist),
            sd_dist = sd(TravDist),
            sd_dur = sd(Duration),
            sd_speed = sd(AvSpeed),
            sd_max = sd(MaxDist),
            se_dist = sd(TravDist)/sqrt(n()),
            se_dur = sd(Duration)/sqrt(n()),
            se_speed = sd(AvSpeed)/sqrt(n()),
            se_max = sd(MaxDist)/sqrt(n())) %>%
  select(-Colony)

# pivot boot_results longer
boot_results_long <- boot_results %>%
  pivot_longer(cols = c(mean_dist, mean_dur, mean_speed, mean_max),
               names_to = "metric",
               values_to = "mean") %>%
  mutate(metric = recode(metric,
                         mean_dist = "Trip Distance (km)",
                         mean_dur = "Trip Duration (hr)",
                         mean_speed = "Average Speed (km/h)",
                         mean_max = "Max Distance (km)"),
         Colony = recode(Colony,
                         CB = "Cow Bay",
                         BR = "Bull Roads"))

# pivot overall means
overall_means_long <- overall_means %>% 
  pivot_longer(cols = c(mean_dist, mean_dur, mean_speed, mean_max),
               names_to = "metric",
               values_to = "mean") %>%
  mutate(metric = recode(metric,
                         mean_dist = "Trip Distance (km)",
                         mean_dur = "Trip Duration (hr)",
                         mean_speed = "Average Speed (km/h)",
                         mean_max = "Max Distance (km)"),
         Colony = recode(Colony,
                         CB = "Cow Bay",
                         BR = "Bull Roads"))

# pivot overall standard error
overall_means_se_long <- overall_means_long %>% 
  pivot_longer(cols = c(se_dist, se_dur, se_speed, se_max),
               names_to = "se_metric",
               values_to = "se") %>%
  mutate(se_metric = recode(se_metric,
                         se_dist = "Trip Distance (km)",
                         se_dur = "Trip Duration (hr)",
                         se_speed = "Average Speed (km/h)",
                         se_max = "Max Distance (km)"),
         Colony = recode(Colony,
                         CB = "Cow Bay",
                         BR = "Bull Roads")) %>%
  mutate(lower = mean - 2*se,
         upper = mean + 2*se) %>%
  filter(metric == se_metric)

# pivot overall standard deviation
overall_means_sd_long <- overall_means_long %>% 
  pivot_longer(cols = c(sd_dist, sd_dur, sd_speed, sd_max),
               names_to = "sd_metric",
               values_to = "sd") %>%
  mutate(sd_metric = recode(sd_metric,
                            sd_dist = "Trip Distance (km)",
                            sd_dur = "Trip Duration (hr)",
                            sd_speed = "Average Speed (km/h)",
                            sd_max = "Max Distance (km)"),
         Colony = recode(Colony,
                         CB = "Cow Bay",
                         BR = "Bull Roads")) %>%
  mutate(lower = mean - sd,
         upper = mean + sd) %>%
  filter(metric == sd_metric)

# pivot SJ means
SJ_means_long <- SJ_means %>%
  pivot_longer(cols = c(mean_dist, mean_dur, mean_speed, mean_max),
               names_to = "metric",
               values_to = "mean") %>%
  mutate(metric = recode(metric,
                         mean_dist = "Trip Distance (km)",
                         mean_dur = "Trip Duration (hr)",
                         mean_speed = "Average Speed (km/h)",
                         mean_max = "Max Distance (km)"))



# plot
deviation_plot <- ggplot(boot_results_long, aes(x = mean, fill = Colony)) +
  geom_histogram(alpha = 0.5, position = "identity") +
  geom_vline(data = overall_means_long, aes(xintercept = mean, col = Colony), 
             linetype = "solid", linewidth = 1) +
  geom_vline(data = overall_means_se_long, aes(xintercept = lower, col = Colony), 
             linetype = "dashed", linewidth = .5) +
  geom_vline(data = overall_means_se_long, aes(xintercept = upper, col = Colony), 
             linetype = "dashed", linewidth = .5) +
  geom_vline(data = overall_means_sd_long, aes(xintercept = lower, col = Colony), 
             linetype = "dotted", linewidth = .75) +
  geom_vline(data = overall_means_sd_long, aes(xintercept = upper, col = Colony), 
             linetype = "dotted", linewidth = .75) +
  facet_wrap(~metric * Colony, scales = "free", ncol = 2) +
  xlab("") +
  ylab("Count") +
  theme_minimal() +
  scale_y_continuous(expand = c(0,0)) + 
  scale_fill_manual(values = c("Cow Bay" = "steelblue4", "Bull Roads" = "orange")) +
  scale_color_manual(values = c("Cow Bay" = "steelblue4", "Bull Roads" = "orange")) +
  theme(legend.title = element_blank(),
        legend.position = "none") 
deviation_plot + ggview::canvas(7, 10)

# export
ggsave("output/nlme/bootstrap/deviation_plot.png", width = 7, height = 10)


# for each metric, what proportion of bootstrapped means lie beyond the standard error
boot_results_long <- boot_results_long %>%
  left_join(overall_means_se_long %>% select(Colony, metric, lower, upper), by = c("Colony", "metric")) %>%
  rename(lower_se = lower, upper_se = upper) %>%
  left_join(overall_means_sd_long %>% select(Colony, metric, lower, upper), by = c("Colony", "metric")) %>%
  rename(lower_sd = lower, upper_sd = upper)

proportions <- boot_results_long %>%
  group_by(Colony, metric) %>%
  summarise(proportion_beyond_se = mean(mean < lower_se | mean > upper_se),
            proportion_beyond_sd = mean(mean < lower_sd | mean > upper_sd)) %>%
  mutate(proportion_within_se = 1 - proportion_beyond_se,
         proportion_within_sd = 1 - proportion_beyond_sd)


# add proportions to plot - left align
deviation_plot2 <- deviation_plot +
  geom_text(data = proportions, aes(x = Inf, y = Inf, label = paste0("SE: ", round(proportion_within_se*100, 1), "%\nSD: ", round(proportion_within_sd*100, 1), "%")), 
            hjust = 1.7, vjust = 1.1, size = 3, fontface = "italic")
deviation_plot2 + ggview::canvas(7, 10)

# export
ggsave("output/nlme/bootstrap/deviation_plot2.png", width = 7, height = 10)

mean(proportions$proportion_within_sd)
mean(proportions$proportion_within_se)
