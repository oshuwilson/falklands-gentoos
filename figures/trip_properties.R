#-------------------------------------------------------------------------------
# Visualise Outputs from Linear Mixed Effects Models
# Joshua Wilson
#-------------------------------------------------------------------------------

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Placement/Writeup")

library(tidyverse)

# list all dive/trip summary data
files <- list.files("output/summaries/", full.names = TRUE)
files

# for each file
for(file in files){
  
  # get name from filename
  name <- str_remove(basename(file), ".rds")
  
  # read in file
  df <- readRDS(file)
  
  # if se is a capitalised column name, make lowercase
  if("SE" %in% colnames(df)){
    df <- df %>%
      rename(se = SE)
  }
  
  # append variable name
  df <- df %>%
    mutate(variable = name)
  
  # combine to other files
  if(file == files[1]){
    data <- df
  } else {
    data <- rbind(data, df)
  }
}

# reorder variables so it runs from trip properties to dive properties
data$variable <- factor(data$variable, levels = c("avspeed", "maxdist", "travdist", "duration",
                                                  "bottime", "divetime", "maxdep", "prop_benthic"))

# recode variables
data <- data %>%
  mutate(variable = recode(variable,
                           avspeed = "Average Speed (km/h)",
                           maxdist = "Maximum Distance from Colony (km)",
                           travdist = "Total Travel Distance (km)",
                           duration = "Trip Duration (hr)",
                           bottime = "Bottom Time (s)",
                           divetime = "Dive Duration (s)",
                           maxdep = "Maximum Depth (m)",
                           prop_benthic = "Proportion of Benthic Dives (%)"))

# recode colony
data <- data %>%
  mutate(Colony = recode(Colony,
                         BR = "Bull Roads",
                         CB = "Cow Bay",
                         SJ = "Steeple Jason Neck"))

# make group a factor
# data$group <- factor(data$group)

# make upper and lower bounds
data <- data %>%
  mutate(lower = mean - sd,
         upper = mean + sd)

# if lower bound exceeds 0, set to 0
data <- data %>%
  mutate(lower = ifelse(lower < 0, 0, lower))

# change season names to 2012 and 2013
data <- data %>%
  mutate(Season = recode(Season,
                         "2012_13" = "2012",
                         "2013_14" = "2013"))

# plot data
p1 <- ggplot(data, aes(x = Colony, y = mean, group = Season, col = Season)) +
  geom_point(position = position_dodge(width = .5)) +
  geom_linerange(aes(ymin = lower, ymax = upper),
                 position = position_dodge(width = .5)) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)), limits = c(0,NA)) +
  facet_wrap(~ variable, scales = "free_y", nrow = 2) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = c("darkred", "steelblue4"))
p1 + ggview::canvas(width = 11, height = 6)

# export
ggsave("output/imagery/trip_and_dive_properties.png", 
       plot = p1, width = 11, height = 6, units = "in", dpi = 300)

# export EPS
ggsave("output/imagery/submission figs/trip_and_dive_properties.eps", 
       plot = p1, width = 11, height = 6, units = "in", dpi = 300)
