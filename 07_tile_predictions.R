### run tiled predictions to avoid using the cluster

library(terra)
library(glmmTMB)
library(tidyverse)

## get the models and the raster tiles 
tmbs <- list.files("/home/hbronnvik/Documents/chapter3/model_fits", pattern = ".rds", full.names = T)
# rasters prepared in 04_env_rasts.R
lyrs <- rast("/home/hbronnvik/Documents/chapter3/env_rasts/full_terr.tif")
# pare down the layers
lyrs <- lyrs[[c("DEM", "obs", "max_blh")]]

# fill in the crossing points
crossing_coasts <- rast("/home/hbronnvik/Documents/chapter3/barrier_mask.tif")
filled_3x <- lapply(lyrs, function(r){
  # fill the empty spaces with 0s 
  r <- classify(r, cbind(NA, 0))
  # then remove them except for sea crossing points
  r <- mask(r, crossing_coasts)
  return(r)
})
lyrs <- rast(filled_3x)
plot(lyrs)

dim(lyrs) # rows, cols, layers, 11867  7654     6
# set splitting parameters
rc <- ceiling(dim(lyrs)[1:2] / c(4,4))
# make 16 tiles
x <- makeTiles(lyrs, rc, filename = "/home/hbronnvik/Documents/chapter3/tiles3x/tile3x_.tif")

# the annotated data for its attributes
ua_locs <- readRDS("/home/hbronnvik/Documents/chapter3/ua_data_60min_20250804.rds")

# ensure the backgrounds are scaled and transformed the same way the model data were
att <- ua_locs %>% 
  # the attributes of the annotated data to scale the full data with
  summarize(m_dem = mean(log(DEM+424), na.rm = T),
            sd_dem = sd(log(DEM+424), na.rm = T),
            m_blh = mean((1/max_blh), na.rm = T),
            sd_blh = sd((1/max_blh), na.rm = T),
            m_obs = mean(1/(1+obs)**2, na.rm = T),
            sd_obs = sd(1/(1+obs)**2, na.rm = T))
# keep the IDs 
IDs <- unique(ua_locs$individual.id)
# then clear space
rm(ua_locs)

## prep the values to predict across
## predict the model across the values
epsilon <- 0.001
tiles <- list.files("/home/hbronnvik/Documents/chapter3/tiles3x", full.names = T)

start_time <- Sys.time()
outputs <- lapply(1:length(tiles), function(y){
  print(paste0("Processing tile ", y, " of ", length(tiles), "."), quote = F)
  tile <- tiles[y]
  # scale the background data
  lb <- tile %>% 
    rast() %>% 
    as.data.frame(na.rm = F, xy = T) %>% 
    mutate(sqrt_dem_z = (log(DEM+424)-att$m_dem)/att$sd_dem,
           max_blh_z = ((1/max_blh)-att$m_blh)/att$sd_blh,
           obs_z = ((1/(1+obs)**2)-att$m_obs)/att$sd_obs,
           # set sl_ and ta_ to means
           sqrt_sl_z = 0,
           ta_z = 0,
           # set the grouping variables to NA as per glmm.TMB recommendation
           stratum_ID = NA,
           individual.id = sample(unique(IDs), n(), replace = T))
  
  # go through each model to predict across that tile
  results <- lapply(tmbs, function(x){
    # the model
    TMB_mod <- readRDS(x)
    new_lb <- lb
    # predict using the models
    new_lb$preds <- predict(TMB_mod, newdata = lb, type = "link")
    new_lb$mod <- names(tmbs)[x]
    gc()
    # convert the log odds to relative selection strengths (RSS)
    new_lb$rss <- exp(new_lb$preds)
    # rasterize
    output <- rast(tile)[[1]]
    values(output) <- new_lb$rss
    names(output) <- gsub("/home/hbronnvik/Documents/chapter3/model_fits/|.rds", "", x)
    # save that tile for each model
    return(output)
  })
  res <- rast(results)
  return(res)
})
Sys.time()-start_time # Time difference of 1.105256 hours

# saveRDS(outputs, file = "/home/hbronnvik/Documents/chapter3/tiled_local_RSS.rds")

preds <- sprc(outputs)
preds <- merge(preds)
plot(preds)

# invert the RSS to be a resistance
resistance <- 1/preds
plot(resistance)
plot(log(resistance))

lapply(1:nlyr(resistance), function(r){
  writeRaster(resistance[[r]], file = paste0("/home/hbronnvik/Documents/chapter3/Raven/resist/res_", names(resistance[[r]]), ".tif"))
})

# check <- rast(list.files("/home/hbronnvik/Documents/chapter3/Raven/Current", pattern = ".asc", full.names = T))
# obs <- log(check[["curconn_res_obsW_curmap"]])
# library(rnaturalearth)
# library(sf)
# world <- ne_coastline(scale = "large", returnclass = "sf") %>% 
#   st_transform(crs(obs))
# plot(obs)
# lines(world, col = "gray50")
# 
# theme_set(theme_classic()+
#             theme(axis.text = element_text(color = "black", size = 12), 
#                   text = element_text(size = 15)))
# colfunc <- colorRampPalette(c("#0081A7", "#0098b0", "#00AFB9", "#7fc4b8", "#FED9B7", "#f7a58f", "#f27e71", "#F07167", "#ee5e53"))
# 
# ggplot() +
#   tidyterra::geom_spatraster(data = obs) +
#   scale_fill_gradientn("log Current", colors = colfunc(100), na.value = NA) +
#   geom_sf(data = world, color = "gray50") +
#   coord_sf(xlim = c(-1771942.75609039, 5005879.4100524), 
#            ylim = c(-3639482.55008231, 6869064.04589582)) +
#   labs(subtitle = "Cumulative current through conspecific resistance")

# adjust the RSS to be on a scale of 0-1 to invert 
# rss <- lapply(1:nlyr(preds), function(p){
#   pred_lyr <- preds[[p]]
#   pred_lyr <- (pred_lyr-min(values(pred_lyr, na.rm = T)))/(max(values(pred_lyr, na.rm = T))-min(values(pred_lyr, na.rm = T)))
#   return(pred_lyr)
# })
# rss <- rast(rss)
# plot(log(1-rss))

