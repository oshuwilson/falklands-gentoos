#-------------------------------------------------------------------------------
# Boosted Regression Trees for each colony
# Joshua Wilson (adapted from Jonathan Handley)
#-------------------------------------------------------------------------------

# rerun dive intensity models with scaled dive intensity for colonies 1 (BR)

rm(list=ls())
setwd("~/OneDrive - University of Southampton/Documents/Placement/Writeup/")

# load libraries 
library(tidyverse)
library(terra) 
library(sf) 
library(geosphere)
library(tidymodels)
library(tidysdm)
library(bonsai)
library(lightgbm)
library(themis)
library(DALEXtra)

# list all files to read in 
files <- list.files(path = "output/habmod/extraction/", pattern = "*.rds", full.names = T, recursive = T)

# choose while file to work with this time
for(file in files[2:3]){
  print(file)
  data <- readRDS(file)
  
  # get study name from file
  study <- tools::file_path_sans_ext(basename(file))
  study
  
  
  #-------------------------------------------------------------------------------
  # Data Preparation
  #-------------------------------------------------------------------------------
  
  # scale intensity column 
  data <- data %>%
    mutate(dive_intensity = dive_intensity / max(dive_intensity)) %>%
    ungroup()
  
  # remove cells over land
  data <- data %>%
    filter(!is.na(depth))
  
  # use hurdle model approach to deal with zeros in the data
  # convert all non-zero presence values to 1
  data <- data %>%
    mutate(pa = ifelse(dive_intensity == 0, 0, 1))
  
  # convert pa to ordered factor
  data <- data %>% mutate(pa = as.factor(pa))
  data$pa <- ordered(data$pa, levels = c("1", "0"))
  
  # subset final data to only include columns of interest for modelling
  df <- data %>%
    select(pa, depth, dist2col, dist2coast, slope)
  
  
  #-------------------------------------------------------------------------------
  # Fit Boosted Regression Trees - Binary 0-1 Model
  #-------------------------------------------------------------------------------
  
  # define BRT
  brt_mod <- boost_tree() %>%
    set_mode("classification") %>%
    set_engine("lightgbm", num_threads = 1 #use lightgbm package
    ) %>%
    set_args(trees = tune(),
             tree_depth = tune(), 
             learn_rate = tune(), 
             min_n = 20) 
  
  # create workflow
  brt_wf <- workflow() %>%
    add_model(brt_mod)
  
  # define hyperparameter values to vary over 
  learn.rate <- c(0.005, 0.01)
  tree.depth <- 1:3
  trees <- seq(500, 2000, 500)
  grid <- expand_grid(learn_rate = learn.rate, tree_depth = tree.depth, trees = trees)
  
  # create cross-validation folds
  folds <- vfold_cv(data = df, v = 10)
  
  # define formula for modelling
  rec <- recipe(pa ~ ., data = df) %>%
    step_downsample(pa)
  
  # update workflow
  brt_wf <- brt_wf %>%
    add_recipe(rec)
  
  # run models with tuning
  tun <- tune_grid(brt_wf,
                   resamples = folds,
                   grid = grid,
                   metrics = sdm_metric_set(),
                   control = control_grid(verbose=T)) 
  
  # get metric scores for each tuning value
  metrics <- collect_metrics(tun, summarize = F)
  
  # extract best model
  best <- show_best(tun, metric = "boyce_cont") %>%
    filter(n == 10)
  
  # get validation metrics for best model
  best_metrics <- collect_metrics(tun, summarize = T) %>%
    filter(trees == best$trees[1],
           tree_depth == best$tree_depth[1],
           learn_rate == best$learn_rate[1])
  
  # set up model to use best hyperparameters
  best_mod <- boost_tree() %>%
    set_engine(engine = "lightgbm") %>%
    set_mode("classification") %>%
    set_args(min_n = 20, trees = best$trees[1], 
             tree_depth = best$tree_depth[1], learn_rate = best$learn_rate[1])
  
  # update workflow
  best_wf <- brt_wf %>%
    update_model(best_mod)
  
  # run best model on all data
  best_fit <- best_wf %>%
    fit(df)
  
  # save model outputs
  saveRDS(best_fit, 
          paste0("output/habmod/models/brt_zero_model_", study, ".rds"))
  saveRDS(best_metrics, 
          paste0("output/habmod/metrics/brt_zero_metrics_", study, ".rds"))
  
  
  #-------------------------------------------------------------------------------
  # Plot Model Predictions and Supplementary Info
  #-------------------------------------------------------------------------------
  
  # read in environmental rasters
  depth <- rast("output/habmod/rasters/depth.tif")
  slope <- rast("output/habmod/rasters/slope.tif")
  dist2coast <- rast("output/habmod/rasters/dist2coast.tif")
  if(study == "bullroad"){
    dist2col <- rast("output/habmod/rasters/dist2BR.tif")
  } else if(study == "cowbay"){
    dist2col <- rast("output/habmod/rasters/dist2CB.tif")
  } else if(study == "steeplejason"){
    dist2col <- rast("output/habmod/rasters/dist2SJ.tif")
  }
  
  # limit to extent of this region
  if(study == "bullroad"){
    ext_rast <- rast("data/rasters/BR1314Gdive_densities.tif")
  }
  if(study == "cowbay"){
    ext_rast <- rast("data/rasters/CB1314Gdive_densities.tif")
  }
  if(study == "steeplejason"){
    ext_rast <- rast("data/rasters/SJ1213Idive_densities.tif")
  }
  e <- ext(ext_rast)
  
  # stack predictors
  preds <- c(depth, slope, dist2coast, dist2col)
  
  # rename predictors
  names(preds) <- c("depth", "slope", "dist2coast", "dist2col")
  
  # crop predictors
  preds <- crop(preds, e)
  
  # predict final model
  pred <- predict_raster(best_fit, preds, type = "prob")
  pred <- pred[[names(pred) == ".pred_1"]]
  plot(pred)
  
  # export prediction
  writeRaster(pred, 
              paste0("output/habmod/predictions/brt_zero_pred_", study, ".tif"),
              overwrite = T)
  
  # get variable importance
  vip <- vip::vi(best_fit)
  vip
  
  # export variable importance scores
  saveRDS(vip,
          paste0("output/habmod/supp/brt_zero_varimp_", study, ".rds"))
  
  #get explainer
  explainer <- explain_tidymodels(model = best_fit, 
                                  data = dplyr::select(df, -pa),
                                  y = as.integer(df$pa),
                                  verbose = T)
  
  #compute partial dependence
  pdps <- model_profile(explainer, 
                        variables = names(df)[!names(df) %in% c("pa")],
                        N = 500)
  
  #extract pdp predictive values
  pdp_ovr <- as_tibble(pdps$agr_profiles) %>%
    rename(x = `_x_`, yhat = `_yhat_`, var = `_vname_`) %>%
    dplyr::select(var, x, yhat) %>%
    mutate(yhat = 1-yhat)
  
  # plot
  ggplot(pdp_ovr, aes(x = x, y = yhat)) +
    geom_line() +
    facet_wrap(~var, scales = "free_x") +
    theme_bw()
  
  # export PDP values
  saveRDS(pdp_ovr, 
          paste0("output/habmod/supp/brt_zero_pdp_", study, ".rds"))
  
  
  #-------------------------------------------------------------------------------
  # Fit Boosted Regression Trees - Dive Intensity Model
  #-------------------------------------------------------------------------------
  
  # remove absences from dataframe
  df2 <- data %>%
    filter(pa == 1)
  
  # select only columns of interest
  df2 <- df2 %>%
    select(dive_intensity, depth, slope, dist2coast, dist2col)
  
  # define BRT - regression this time
  brt_mod <- boost_tree() %>%
    set_mode("regression") %>%
    set_engine("lightgbm", num_threads = 1 #use lightgbm package
    ) %>%
    set_args(trees = tune(),
             tree_depth = tune(), 
             learn_rate = tune(), 
             min_n = 20) 
  
  # create workflow
  brt_wf <- workflow() %>%
    add_model(brt_mod)
  
  # define hyperparameter values to vary over 
  learn.rate <- c(0.005, 0.01)
  tree.depth <- 1:3
  trees <- seq(500, 2000, 500)
  grid <- expand_grid(learn_rate = learn.rate, tree_depth = tree.depth, trees = trees)
  
  # create cross-validation folds
  folds <- vfold_cv(data = df2, v = 10)
  
  # define formula for modelling
  rec <- recipe(dive_intensity ~ ., data = df2) 
  
  # update workflow
  brt_wf <- brt_wf %>%
    add_recipe(rec)
  
  # run models with tuning
  tun <- tune_grid(brt_wf,
                   resamples = folds,
                   grid = grid,
                   metrics = metric_set(rmse, rsq, mae),
                   control = control_grid(verbose=T)) 
  
  # get metric scores for each tuning value
  metrics <- collect_metrics(tun, summarize = F)
  
  # extract best model
  best <- show_best(tun, metric = "mae") %>%
    filter(n == 10)
  
  # get validation metrics for best model
  best_metrics <- collect_metrics(tun, summarize = T) %>%
    filter(trees == best$trees[1],
           tree_depth == best$tree_depth[1],
           learn_rate == best$learn_rate[1])
  
  # set up model to use best hyperparameters
  best_mod <- boost_tree() %>%
    set_engine(engine = "lightgbm") %>%
    set_mode("regression") %>%
    set_args(min_n = 20, trees = best$trees[1], 
             tree_depth = best$tree_depth[1], learn_rate = best$learn_rate[1])
  
  # update workflow
  best_wf <- brt_wf %>%
    update_model(best_mod)
  
  # run best model on all data
  best_fit <- best_wf %>%
    fit(df2)
  
  # save model outputs
  saveRDS(best_fit, 
          paste0("output/habmod/models/brt_intensity_model_", study, ".rds"))
  saveRDS(best_metrics, 
          paste0("output/habmod/metrics/brt_intensity_metrics_", study, ".rds"))
  
  
  #-------------------------------------------------------------------------------
  # Plot Model Predictions and Supplementary Info
  #-------------------------------------------------------------------------------
  
  # read in environmental rasters
  depth <- rast("output/habmod/rasters/depth.tif")
  slope <- rast("output/habmod/rasters/slope.tif")
  dist2coast <- rast("output/habmod/rasters/dist2coast.tif")
  if(study == "bullroad"){
    dist2col <- rast("output/habmod/rasters/dist2BR.tif")
  } else if(study == "cowbay"){
    dist2col <- rast("output/habmod/rasters/dist2CB.tif")
  } else if(study == "steeplejason"){
    dist2col <- rast("output/habmod/rasters/dist2SJ.tif")
  }
  
  # limit to extent of this region
  if(study == "bullroad"){
    ext_rast <- rast("data/rasters/BR1314Gdive_densities.tif")
  }
  if(study == "cowbay"){
    ext_rast <- rast("data/rasters/CB1314Gdive_densities.tif")
  }
  if(study == "steeplejason"){
    ext_rast <- rast("data/rasters/SJ1213Idive_densities.tif")
  }
  e <- ext(ext_rast)
  
  # stack predictors
  preds <- c(depth, slope, dist2coast, dist2col)
  
  # rename predictors
  names(preds) <- c("depth", "slope", "dist2coast", "dist2col")
  
  # crop predictors
  preds <- crop(preds, e)
  
  # predict final model
  pred <- predict_raster(best_fit, preds)
  plot(pred)
  
  # export prediction
  writeRaster(pred, 
              paste0("output/habmod/predictions/brt_intensity_pred_", study, ".tif"),
              overwrite = T)
  
  # get variable importance
  vip <- vip::vi(best_fit)
  vip
  
  # export variable importance scores
  saveRDS(vip,
          paste0("output/habmod/supp/brt_intensity_varimp_", study, ".rds"))
  
  #get explainer
  explainer <- explain_tidymodels(model = best_fit, 
                                  data = dplyr::select(df2, -dive_intensity),
                                  y = df2$dive_intensity,
                                  verbose = T)
  
  #compute partial dependence
  pdps <- model_profile(explainer, 
                        variables = names(df2)[!names(df2) %in% c("dive_intensity")],
                        N = 500)
  
  #extract pdp predictive values
  pdp_ovr <- as_tibble(pdps$agr_profiles) %>%
    rename(x = `_x_`, yhat = `_yhat_`, var = `_vname_`) %>%
    dplyr::select(var, x, yhat)
  
  # plot
  ggplot(pdp_ovr, aes(x = x, y = yhat)) +
    geom_line() +
    facet_wrap(~var, scales = "free_x") +
    theme_bw()
  
  # export PDP values
  saveRDS(pdp_ovr, 
          paste0("output/habmod/supp/brt_intensity_pdp_", study, ".rds"))
  
  
  #-------------------------------------------------------------------------------
  # Combine Model Predictions
  #-------------------------------------------------------------------------------
  
  # read in predictions
  pred_zero <- rast(paste0("output/habmod/predictions/brt_zero_pred_", study, ".tif"))
  pred_intensity <- rast(paste0("output/habmod/predictions/brt_intensity_pred_", study, ".tif"))
  plot(pred_zero)
  plot(pred_intensity)
  
  # multiply together
  final_pred <- pred_zero * pred_intensity
  plot(final_pred)
  
  # export
  writeRaster(final_pred, 
              paste0("output/habmod/predictions/brt_final_pred_", study, ".tif"),
              overwrite = T)
}
