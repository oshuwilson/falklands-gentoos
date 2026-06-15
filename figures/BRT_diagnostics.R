#-------------------------------------------------------------------------------
# Visualise the Habitat Relationships from Boosted Regression Trees
# Joshua Wilson
#-------------------------------------------------------------------------------

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Placement/Writeup/")

# load libraries 
{
  library(tidyverse)
  library(cowplot)
}


#-------------------------------------------------------------------------------
# Variable Importance
#-------------------------------------------------------------------------------

# read in variable importance scores for each colony and model
br_zero <- readRDS("output/habmod/supp/brt_zero_varimp_bullroad.rds")
br_int <- readRDS("output/habmod/supp/brt_intensity_varimp_bullroad.rds")

cb_zero <- readRDS("output/habmod/supp/brt_zero_varimp_cowbay.rds")
cb_int <- readRDS("output/habmod/supp/brt_intensity_varimp_cowbay.rds")

sj_zero <- readRDS("output/habmod/supp/brt_zero_varimp_steeplejason.rds")
sj_int <- readRDS("output/habmod/supp/brt_intensity_varimp_steeplejason.rds")

# add colony and model labels
br_zero <- br_zero %>%
  mutate(colony = "Bull Roads", model = "Presence-Absence")
br_int <- br_int %>%
  mutate(colony = "Bull Roads", model = "Dive Intensity")

cb_zero <- cb_zero %>%
  mutate(colony = "Cow Bay", model = "Presence-Absence")
cb_int <- cb_int %>%
  mutate(colony = "Cow Bay", model = "Dive Intensity")

sj_zero <- sj_zero %>%
  mutate(colony = "Steeple Jason Neck", model = "Presence-Absence")
sj_int <- sj_int %>%
  mutate(colony = "Steeple Jason Neck", model = "Dive Intensity")

# merge together
varimp <- rbind(br_zero, br_int, cb_zero, cb_int, sj_zero, sj_int)

# recode predictors
varimp <- varimp %>%
  mutate(Variable = recode(Variable, 
                           "dist2col" = "Distance to Colony",
                           "dist2coast" = "Distance to Coast",
                           "depth" = "Depth",
                           "slope" = "Slope"))

# set model so presence-absence comes first
varimp$model <- factor(varimp$model, levels = c("Presence-Absence", "Dive Intensity"))

# plot
p1 <- ggplot(varimp, aes(x = reorder(Variable, -Importance), y = Importance)) +
  geom_bar(stat = "identity") +
  facet_grid(model ~ colony) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Variable Importance") +
  xlab("") +
  scale_y_continuous(limits = c(0,1)) +
  theme(panel.border = element_rect(colour = "black", fill=NA))
p1 + ggview::canvas(width = 10, height = 8)

# export
ggsave("output/imagery/brt_varimp.png", p1, width = 10, height = 8, dpi = 300)


#-------------------------------------------------------------------------------
# Partial Dependence Plots
#-------------------------------------------------------------------------------

# clean up
rm(list=ls())

# read in pdp data for each colony and model
br_zero_pdp <- readRDS("output/habmod/supp/brt_zero_pdp_bullroad.rds")
br_int_pdp <- readRDS("output/habmod/supp/brt_intensity_pdp_bullroad.rds")

cb_zero_pdp <- readRDS("output/habmod/supp/brt_zero_pdp_cowbay.rds")
cb_int_pdp <- readRDS("output/habmod/supp/brt_intensity_pdp_cowbay.rds")

sj_zero_pdp <- readRDS("output/habmod/supp/brt_zero_pdp_steeplejason.rds")
sj_int_pdp <- readRDS("output/habmod/supp/brt_intensity_pdp_steeplejason.rds")

# append colony and model labels
br_zero_pdp <- br_zero_pdp %>%
  mutate(colony = "Bull Roads", model = "Presence-Absence")
br_int_pdp <- br_int_pdp %>%
  mutate(colony = "Bull Roads", model = "Dive Intensity")

cb_zero_pdp <- cb_zero_pdp %>%
  mutate(colony = "Cow Bay", model = "Presence-Absence")
cb_int_pdp <- cb_int_pdp %>%
  mutate(colony = "Cow Bay", model = "Dive Intensity")

sj_zero_pdp <- sj_zero_pdp %>%
  mutate(colony = "Steeple Jason Neck", model = "Presence-Absence")
sj_int_pdp <- sj_int_pdp %>%
  mutate(colony = "Steeple Jason Neck", model = "Dive Intensity")

# combine together
pdp <- rbind(br_zero_pdp, br_int_pdp, cb_zero_pdp, cb_int_pdp, sj_zero_pdp, sj_int_pdp)

# remove unnecessary predictors
pdp <- pdp %>%
  filter(var %in% c("depth", "dist2coast", "dist2col", "slope"))



# reorder predictors to order of importance (colony, coast, depth, slope)
pdp$var <- factor(pdp$var, levels = c("dist2col", "dist2coast", "depth", "slope"))

# recode predictors
pdp <- pdp %>%
  mutate(var = recode(var, 
                      "dist2col" = "Distance to Colony",
                      "dist2coast" = "Distance to Coast",
                      "depth" = "Depth",
                      "slope" = "Slope"))

# invert depth
pdp <- pdp %>%
  mutate(x = ifelse(var == "Depth", -x, x))

# split into presence-absence and dive intensity models
pdp_pa <- pdp %>% filter(model == "Presence-Absence")
pdp_int <- pdp %>% filter(model == "Dive Intensity")

# remove extreme values where lines flatline from presence-absence model
pdp_pa <- pdp_pa %>%
  filter(var != "Distance to Colony" | x < 70000) %>%
  filter(var != "Distance to Coast" | x < 65000) %>%
  filter(var != "Depth" | x < 250 & x >= 0) %>%
  filter(var != "Slope" | x < 1.4)

# divide distance metrics to kilometres
pdp_pa <- pdp_pa %>%
  mutate(x = ifelse(var %in% c("Distance to Colony", "Distance to Coast"), x / 1000, x))

# convert slope from radians to degrees
pdp_pa <- pdp_pa %>%
  mutate(x = ifelse(var == "Slope", x * (180 / pi), x))

# recode variables to include units
pdp_pa <- pdp_pa %>%
  mutate(var = recode(var,
                      "Distance to Colony" = "Distance to Colony (km)",
                      "Distance to Coast" = "Distance to Coast (km)",
                      "Depth" = "Depth (m)",
                      "Slope" = "Slope (°)"))

# plot presence-absence results
p2 <- ggplot(pdp_pa, aes(x = x, y = yhat)) +
  geom_line(aes(col = colony), lwd = 1) +
  facet_wrap(~var, scales = "free_x", nrow = 4) + 
  theme_minimal() +
  ylab("Partial Effect") +
  xlab("") + 
  ggtitle("Presence-Absence Model") +
  theme(panel.background = element_rect(fill = NA, colour = "black"),
        plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("#DB941A", "#007E8F", "#DB5943"), name = "") 
p2  

# limit dive intensity pdp to same x-axis limits as presence-absence model
pdp_int <- pdp_int %>%
  filter(var != "Distance to Colony" | x < 40000) %>%
  filter(var != "Distance to Coast" | x < 45000) %>%
  filter(var != "Depth" | x < 150 & x >= 0) %>%
  filter(var != "Slope" | x < 1.5)

# revalue values below 0 to 0
pdp_int <- pdp_int %>%
  mutate(yhat = ifelse(yhat < 0, 0, yhat))

# divide distance metrics to kilometres
pdp_int <- pdp_int %>%
  mutate(x = ifelse(var %in% c("Distance to Colony", "Distance to Coast"), x / 1000, x))

# convert slope from radians to degrees
pdp_int <- pdp_int %>%
  mutate(x = ifelse(var == "Slope", x * (180 / pi), x))

# recode variables to include units
pdp_int <- pdp_int %>%
  mutate(var = recode(var,
                      "Distance to Colony" = "Distance to Colony (km)",
                      "Distance to Coast" = "Distance to Coast (km)",
                      "Depth" = "Depth (m)",
                      "Slope" = "Slope (°)"))


# plot dive intensity results
p3 <- ggplot(pdp_int, aes(x = x, y = yhat)) +
  geom_line(aes(col = colony), lwd = 1) +
  facet_wrap(~var, scales = "free_x", nrow = 4) +
  theme_minimal() +
  ylab("Partial Effect") +
  xlab("") +
  ggtitle("Dive Intensity Model") +
  theme(panel.background = element_rect(fill = NA, colour = "black"),
        plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("#DB941A", "#007E8F", "#DB5943"), name = "") 
p3

# get legend 
legend <- get_legend(p2)
plot(legend)

# remove legend from plots
p2 <- p2 + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "none") + ylab("")

# plot together
plots <- plot_grid(p2, p3)

# add legend
grid <- plot_grid(plots, legend, ncol = 2, rel_widths = c(1, 0.2))
grid + ggview::canvas(width = 10, height = 12)

# export
ggsave("output/imagery/brt_pdp.png", grid, width = 10, height = 12, dpi = 300)


# export as EPS
ggsave("output/imagery/submission figs/brt_pdp.eps", grid, width = 10, height = 12, dpi = 300)
