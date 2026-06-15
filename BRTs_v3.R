#-------------------------------------------------------------------------------
# Boosted Regression Trees for each colony - incorporating chl and sst
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
library(spatialsample)
library(bonsai)
library(lightgbm)
library(themis)
library(DALEXtra)

# list all files to read in 
files <- list.files(path = "output/habmod/extraction/", pattern = "*.rds", full.names = T, recursive = T)
files <- files[c(2:3, 5:7)]
files

# loop over files for each colony
for(i in 5){
  
  if(i %in% c(1,3)){
    data1 <- readRDS(files[i])
    data2 <- readRDS(files[i+1])
    data <- rbind(data1, data2)
    rm(data1, data2)
  } else {
     data <- readRDS(files[i])
  }
  
  # get study name from file
  study <- tools::file_path_sans_ext(basename(files[i]))
  
  # remove numbers from study name
  study <- gsub("[0-9]", "", study)
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
    select(pa, depth, dist2col, dist2coast, slope, sst, chl, x, y)
  
  
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
  learn.rate <- c(0.005, 0.01, 0.5)
  tree.depth <- c(1, 3, 5)
  trees <- seq(1000, 5000, 1000)
  grid <- expand_grid(learn_rate = learn.rate, tree_depth = tree.depth, trees = trees)
  
  # create spatial cross-validation folds
  df_sf <- st_as_sf(df, coords = c("x", "y"), crs = 4326)
  folds <- spatial_block_cv(data = df_sf, v = 10)
  
  # define formula for modelling
  rec <- recipe(pa ~ depth + dist2col + dist2coast + slope + sst + chl, data = df) %>%
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
    filter(n == max(n))
  
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
          paste0("output/habmod/models/brt_zero_model_", study, "2.rds"))
  saveRDS(best_metrics, 
          paste0("output/habmod/metrics/brt_zero_metrics_", study, "2.rds"))
  
  
  #-------------------------------------------------------------------------------
  # Plot Model Predictions and Supplementary Info
  #-------------------------------------------------------------------------------
  
  # read in environmental rasters
  depth <- rast("output/habmod/rasters/depth.tif")
  slope <- rast("output/habmod/rasters/slope.tif")
  dist2coast <- rast("output/habmod/rasters/dist2coast.tif")
  if(study == "bullroad"){
    dist2col <- rast("output/habmod/rasters/dist2BR.tif")
    sst13 <- rast("output/habmod/rasters/br_sst_13.tif")
    sst14 <- rast("output/habmod/rasters/br_sst_14.tif")
    chl13 <- rast("output/habmod/rasters/br_chl_13.tif")
    chl14 <- rast("output/habmod/rasters/br_chl_14.tif")
  } else if(study == "cowbay"){
    dist2col <- rast("output/habmod/rasters/dist2CB.tif")
    sst13 <- rast("output/habmod/rasters/cb_sst_13.tif")
    sst14 <- rast("output/habmod/rasters/cb_sst_14.tif")
    chl13 <- rast("output/habmod/rasters/cb_chl_13.tif")
    chl14 <- rast("output/habmod/rasters/cb_chl_14.tif")
  } else if(study == "steeplejason"){
    dist2col <- rast("output/habmod/rasters/dist2SJ.tif")
    sst13 <- rast("output/habmod/rasters/sj_sst_13.tif")
    chl13 <- rast("output/habmod/rasters/sj_chl_13.tif")
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
  preds13 <- c(depth, slope, dist2coast, dist2col, sst13, chl13)
  #preds14 <- c(depth, slope, dist2coast, dist2col, sst14, chl14)
  
  # rename predictors
  names(preds13) <- c("depth", "slope", "dist2coast", "dist2col", "sst", "chl")
  #names(preds14) <- c("depth", "slope", "dist2coast", "dist2col", "sst", "chl")
  
  # crop predictors
  preds13 <- crop(preds13, e)
  #preds14 <- crop(preds14, e)
  
  # predict final model
  pred13 <- predict_raster(best_fit, preds13, type = "prob")
  pred13 <- pred13[[names(pred13) == ".pred_1"]]
  plot(pred13)
  
  # pred14 <- predict_raster(best_fit, preds14, type = "prob")
  # pred14 <- pred14[[names(pred14) == ".pred_1"]]
  # plot(pred14)
  
  # export prediction
  writeRaster(pred13, 
              paste0("output/habmod/predictions/brt_zero_pred_", study, "_13.tif"),
              overwrite = T)
  # writeRaster(pred14, 
  #             paste0("output/habmod/predictions/brt_zero_pred_", study, "_14.tif"),
  #             overwrite = T)
  
  # get variable importance
  vip <- vip::vi(best_fit)
  vip
  
  # export variable importance scores
  saveRDS(vip,
          paste0("output/habmod/supp/brt_zero_varimp_", study, "2.rds"))
  
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
          paste0("output/habmod/supp/brt_zero_pdp_", study, "2.rds"))
  
  
  #-------------------------------------------------------------------------------
  # Fit Boosted Regression Trees - Dive Intensity Model
  #-------------------------------------------------------------------------------
  
  # remove absences from dataframe
  df2 <- data %>%
    filter(pa == 1)
  
  # select only columns of interest
  df2 <- df2 %>%
    select(dive_intensity, depth, slope, dist2coast, chl, sst, dist2col, x, y)
  
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
  learn.rate <- c(0.005, 0.01, 0.5)
  tree.depth <- c(1, 3, 5)
  trees <- seq(1000, 5000, 1000)
  grid <- expand_grid(learn_rate = learn.rate, tree_depth = tree.depth, trees = trees)
  
  # create spatial cross-validation folds
  df_sf <- st_as_sf(df2, coords = c("x", "y"), crs = 4326)
  folds <- spatial_block_cv(data = df_sf, v = 10)
  
  # define formula for modelling
  rec <- recipe(dive_intensity ~ dist2col + dist2coast + depth + slope +
                  chl + sst, data = df2) 
  
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
    filter(n == max(n))
  
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
          paste0("output/habmod/models/brt_intensity_model_", study, "2.rds"))
  saveRDS(best_metrics, 
          paste0("output/habmod/metrics/brt_intensity_metrics_", study, "2.rds"))
  
  
  #-------------------------------------------------------------------------------
  # Plot Model Predictions and Supplementary Info
  #-------------------------------------------------------------------------------
  
  # read in environmental rasters
  depth <- rast("output/habmod/rasters/depth.tif")
  slope <- rast("output/habmod/rasters/slope.tif")
  dist2coast <- rast("output/habmod/rasters/dist2coast.tif")
  if(study == "bullroad"){
    dist2col <- rast("output/habmod/rasters/dist2BR.tif")
    sst13 <- rast("output/habmod/rasters/br_sst_13.tif")
    sst14 <- rast("output/habmod/rasters/br_sst_14.tif")
    chl13 <- rast("output/habmod/rasters/br_chl_13.tif")
    chl14 <- rast("output/habmod/rasters/br_chl_14.tif")
  } else if(study == "cowbay"){
    dist2col <- rast("output/habmod/rasters/dist2CB.tif")
    sst13 <- rast("output/habmod/rasters/cb_sst_13.tif")
    sst14 <- rast("output/habmod/rasters/cb_sst_14.tif")
    chl13 <- rast("output/habmod/rasters/cb_chl_13.tif")
    chl14 <- rast("output/habmod/rasters/cb_chl_14.tif")
  } else if(study == "steeplejason"){
    dist2col <- rast("output/habmod/rasters/dist2SJ.tif")
    sst13 <- rast("output/habmod/rasters/sj_sst_13.tif")
    chl13 <- rast("output/habmod/rasters/sj_chl_13.tif")
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
  preds13 <- c(depth, slope, dist2coast, dist2col, sst13, chl13)
  #preds14 <- c(depth, slope, dist2coast, dist2col, sst14, chl14)
  
  # rename predictors
  names(preds13) <- c("depth", "slope", "dist2coast", "dist2col", "sst", "chl")
  #names(preds14) <- c("depth", "slope", "dist2coast", "dist2col", "sst", "chl")
  
  # crop predictors
  preds13 <- crop(preds13, e)
  #preds14 <- crop(preds14, e)
  
  # predict final model
  pred13 <- predict_raster(best_fit, preds13)
  #pred14 <- predict_raster(best_fit, preds14)
  plot(pred13)
  #plot(pred14)
  
  # export prediction
  writeRaster(pred13, 
              paste0("output/habmod/predictions/brt_intensity_pred_", study, "_13.tif"),
              overwrite = T)
  # writeRaster(pred14, 
  #             paste0("output/habmod/predictions/brt_intensity_pred_", study, "_14.tif"),
  #             overwrite = T)
  
  # get variable importance
  vip <- vip::vi(best_fit)
  vip
  
  # export variable importance scores
  saveRDS(vip,
          paste0("output/habmod/supp/brt_intensity_varimp_", study, "2.rds"))
  
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
          paste0("output/habmod/supp/brt_intensity_pdp_", study, "2.rds"))
  
  
  #-------------------------------------------------------------------------------
  # Combine Model Predictions
  #-------------------------------------------------------------------------------
  
  # read in predictions
  pred_zero13 <- rast(paste0("output/habmod/predictions/brt_zero_pred_", study, "_13.tif"))
  pred_intensity13 <- rast(paste0("output/habmod/predictions/brt_intensity_pred_", study, "_13.tif"))
  
  # pred_zero14 <- rast(paste0("output/habmod/predictions/brt_zero_pred_", study, "_14.tif"))
  # pred_intensity14 <- rast(paste0("output/habmod/predictions/brt_intensity_pred_", study, "_14.tif"))
  
  # multiply together
  final_pred13 <- pred_zero13 * pred_intensity13
  plot(final_pred13)
  
  # final_pred14 <- pred_zero14 * pred_intensity14
  # plot(final_pred14)
  
  # export
  writeRaster(final_pred13, 
              paste0("output/habmod/predictions/brt_final_pred_", study, "_13.tif"),
              overwrite = T)
  # writeRaster(final_pred14, 
  #             paste0("output/habmod/predictions/brt_final_pred_", study, "_14.tif"),
  #             overwrite = T)
}
