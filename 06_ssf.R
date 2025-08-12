### model route selection on the basis of 1, 2, or 3 variables and create a resistance surface.
### Hester Bronnvik
### 2025-03-27

# We will prepare a resistance layer as 1-selection by modeling route selection based on
# an increasing number of predictors, predicting from that model across the flyway, and
# then save out the result for Circuitscape.

## 1. prep environment
library(terra)
library(glmmTMB)
library(performance)
library(tidyverse)
theme_set(theme_classic()+
            theme(axis.text = element_text(color = "black", size = 12), 
                  text = element_text(size = 15),
                  strip.background = element_rect(color = "white"),
                  strip.text = element_text(size = 15)))
colfunc <- colorRampPalette(c("#0081A7", "#0098b0", "#00AFB9", "#7fc4b8", "#FED9B7", "#f7a58f", "#f27e71", "#F07167", "#ee5e53"))

## 2. access the migration data that have been burst and their alternative points (annotated)
a_data <- readRDS("/home/hbronnvik/Documents/chapter3/ua_data_60min_20250812.rds")

# the days classified as migratory, start and end migration dates, and outcome
metad <- readRDS("/home/hbronnvik/Documents/chapter3/migration_dates_2025-03-08.rds")
# the ages by track
ages <- metad %>% 
  group_by(trackID) %>% 
  slice(1) %>% 
  dplyr::select(trackID, journey_number, track_status) %>% 
  rename(track = trackID,
         age = journey_number)

# prep the data
a_data <- a_data %>% 
  separate(track, c("id", "season", "year"), sep = "_", remove = F) %>% 
  dplyr::select(-id, -year) %>% 
  # add on age
  left_join(ages) %>% 
  # only consider one season for now and use juvenile decisions
  filter(season == "fall" & age == 1 & track_status == "complete") %>% 
  mutate(#age_z = scale(age)[,1],
         stratum_ID = as.factor(stratum),
         stratum_ID = as.numeric(stratum_ID),
         individual.id = as.numeric(individual.id))

# IDs <- unique(a_data$individual.id)
# saveRDS(IDs, file = "/home/hbronnvik/Documents/chapter3/adata_ids.rds")

## 3. build the models
# make a list of model formulas
# f <- c(# elevation alone
#   used ~ -1 + sqrt_dem_z + sqrt_sl_z + ta_z + (1|stratum_ID) + (0 + sqrt_dem_z|individual.id),
#   # elevation + pbl
#   # used ~ -1 + sqrt_dem_z + log_blh_z + sqrt_sl_z + ta_z + (1|stratum_ID) + (0 + sqrt_dem_z|individual.id) + (0 + log_blh_z|individual.id),
#   # elevation + grass/crop/wetland
#   # used ~ -1 + sqrt_dem_z + forage_z + sqrt_sl_z + ta_z + (1|stratum_ID) + (0 + sqrt_dem_z|individual.id) + (0 + forage_z|individual.id),
#   # elevation + aspect 
#   used ~ -1 + sqrt_dem_z + aspect_z + sqrt_sl_z + ta_z + (1|stratum_ID) + (0 + sqrt_dem_z|individual.id) + (0 + aspect_z|individual.id),
#   # elevation + conspecifics
#   used ~ -1 + sqrt_dem_z + obs_z + sqrt_sl_z + ta_z + (1|stratum_ID) + (0 + sqrt_dem_z|individual.id) + (0 + obs_z|individual.id),
#   # elevation + grass/crop/wetland + aspect
#   # used ~ -1 + sqrt_dem_z + forage_z + aspect_z + sqrt_sl_z + ta_z + (1|stratum_ID) + (0 + sqrt_dem_z|individual.id) + (0 + forage_z|individual.id) + (0 + aspect_z|individual.id),
#   # elevation + conspecifics + aspect
#   used ~ -1 + sqrt_dem_z + obs_z + aspect_z + sqrt_sl_z + ta_z + (1|stratum_ID) + (0 + sqrt_dem_z|individual.id) + (0 + obs_z|individual.id) + (0 + aspect_z|individual.id),
#   # aspect alone
#   used ~ -1 + aspect_z + sqrt_sl_z + ta_z + (1|stratum_ID) + (0 + aspect_z|individual.id),
#   # conspecifics alone
#   used ~ -1 + obs_z + sqrt_sl_z + ta_z + (1|stratum_ID) + (0 + obs_z|individual.id),
#   # aspect + conspecifics
#   used ~ -1 + aspect_z + obs_z + sqrt_sl_z + ta_z + (1|stratum_ID) + (0 + aspect_z|individual.id) + (0 + obs_z|individual.id)
# ) 
# names(f) <- c("dem", "dem_asp", "dem_obs", "dem_obs_asp", "asp", "obs", "asp_obs")

# use a categorical aspect and remove ID from the slopes to compensate
# f <- c(
#   # elevation alone
#   used ~ -1 + sqrt_dem_z + sqrt_sl_z + ta_z + (1|stratum_ID),
#   # elevation + aspect 
#   used ~ -1 + sqrt_dem_z + compass + sqrt_sl_z + ta_z + (1|stratum_ID),
#   # elevation + conspecifics
#   used ~ -1 + sqrt_dem_z + obs_z + sqrt_sl_z + ta_z + (1|stratum_ID),
#   # elevation + conspecifics + aspect
#   used ~ -1 + sqrt_dem_z + obs_z + compass + sqrt_sl_z + ta_z + (1|stratum_ID),
#   # aspect alone
#   used ~ -1 + compass + sqrt_sl_z + ta_z + (1|stratum_ID),
#   # conspecifics alone
#   used ~ -1 + obs_z + sqrt_sl_z + ta_z + (1|stratum_ID),
#   # aspect + conspecifics
#   used ~ -1 + compass + obs_z + sqrt_sl_z + ta_z + (1|stratum_ID),
#   # PBL height alone
#   used ~ -1 + pnt_blh_z + sqrt_sl_z + ta_z + (1|stratum_ID),
#   # med PBL height alone
#   # used ~ -1 + med_blh_z + sqrt_sl_z + ta_z + (1|stratum_ID),
#   # elevation + PBL
#   used ~ -1 + sqrt_dem_z + pnt_blh_z + sqrt_sl_z + ta_z + (1|stratum_ID),
#   # aspect + PBL
#   used ~ -1 + pnt_blh_z + compass + sqrt_sl_z + ta_z + (1|stratum_ID),
#   # conspecifics + PBL
#   used ~ -1 + pnt_blh_z + obs_z + sqrt_sl_z + ta_z + (1|stratum_ID),
#   # aspect + conspecifics + PBL
#   used ~ -1 + pnt_blh_z + obs_z + compass + sqrt_sl_z + ta_z + (1|stratum_ID),
#   # elevation + conspecifics + PBL
#   used ~ -1 + pnt_blh_z + obs_z + sqrt_dem_z + sqrt_sl_z + ta_z + (1|stratum_ID),
#   # elevation + aspect + conspecifics + PBL
#   used ~ -1 + pnt_blh_z + obs_z + compass + sqrt_sl_z + ta_z + (1|stratum_ID)
# ) 
# names(f) <- c("dem", "dem_asp", "dem_obs", "dem_obs_asp", "asp", "obs", "asp_obs", "blh", 
#               "dem_blh", "asp_blh", "obs_blh", "asp_obs_blh", "blh_obs_dem", "dem_asp_obs_blh")

# a_data$obs_z <- NULL
# colnames(a_data)[grepl("obs_num_z", colnames(a_data))] <- "obs_z"

f <- c(
  # elevation alone
  used ~ -1 + sqrt_dem_z + sqrt_sl_z + ta_z + (1|stratum_ID),
  # conspecifics alone
  used ~ -1 + obs_z + sqrt_sl_z + ta_z + (1|stratum_ID),
  # PBL height alone
  used ~ -1 + max_blh_z + sqrt_sl_z + ta_z + (1|stratum_ID),
  # elevation + conspecifics
  used ~ -1 + sqrt_dem_z + obs_z + sqrt_sl_z + ta_z + (1|stratum_ID),
  # elevation + PBL
  used ~ -1 + sqrt_dem_z + max_blh_z + sqrt_sl_z + ta_z + (1|stratum_ID),
  # conspecifics + PBL
  used ~ -1 + max_blh_z + obs_z + sqrt_sl_z + ta_z + (1|stratum_ID),
  # elevation + conspecifics + PBL
  used ~ -1 + max_blh_z + obs_z + sqrt_dem_z + sqrt_sl_z + ta_z + (1|stratum_ID)
) 
names(f) <- c("dem", "obs", "blh", "dem_obs", "dem_blh", "obs_blh", "blh_obs_dem")

# the 4x files, 15 models using DEM, aspect in 16 bins, land cover, and maximum BLH
items <- c("max_blh_z", "sqrt_dem_z", "foo_z", "aspect")
combs <- unlist(lapply(1:length(items), function(k) {
  combn(items, k, simplify = FALSE)
}), recursive = FALSE)

f <- c(
  # elevation alone
  used ~ -1 + sqrt_dem_z + sqrt_sl_z + ta_z + (1|stratum_ID),
  # conspecifics alone
  used ~ -1 + obs_z + sqrt_sl_z + ta_z + (1|stratum_ID),
  # PBL height alone
  used ~ -1 + max_blh_z + sqrt_sl_z + ta_z + (1|stratum_ID),
  # elevation + conspecifics
  used ~ -1 + sqrt_dem_z + obs_z + sqrt_sl_z + ta_z + (1|stratum_ID),
  # elevation + PBL
  used ~ -1 + sqrt_dem_z + max_blh_z + sqrt_sl_z + ta_z + (1|stratum_ID),
  # conspecifics + PBL
  used ~ -1 + max_blh_z + obs_z + sqrt_sl_z + ta_z + (1|stratum_ID),
  # elevation + conspecifics + PBL
  used ~ -1 + max_blh_z + obs_z + sqrt_dem_z + sqrt_sl_z + ta_z + (1|stratum_ID)
) 
names(f) <- c("dem", "obs", "blh", "dem_obs", "dem_blh", "obs_blh", "blh_obs_dem")


# loop through the model formulas
tmbs <- lapply(1:length(f), function(l){
  print(names(f)[l])
  # begin by setting up rules for how SD will be handled on random variables
  # if(nchar(names(f)[l]) == 3){
  #   # for the single-variable models
  #   m <- c(NA, 1)
  #   s <- c(log(1e3), 0)
  # }else{
  #   if(nchar(names(f)[l]) %in% c(7, 8)){
  #     # for the 2-variable model
  #     m <- c(NA, 1, 1)
  #     s <- c(log(1e3), 0, 0)
  #   }else{
  #     # for the 3-variable model
  #     m <- c(NA, 1, 1, 1)
  #     s <- c(log(1e3), 0, 0, 0)
  #   }
  # }
  # without ID on the slopes:
  m <- c(NA)
  s <- c(log(1e3))
  
  # set up the model
  tmb_struc <- glmmTMB(f[[l]],
                       family = poisson,
                       data = a_data, doFit = FALSE,
                       #Tell glmmTMB not to change the first standard deviation, all other values are freely estimated (and are different from each other)
                       map = list(theta = factor(m)), # 3 is the n of random slopes
                       #Set the value of the standard deviation of the first random effect (here (1|startum_ID)):
                       start = list(theta = s)) #add a 0 for each random slope. in this case, 2
  # run the model
  TMB <- glmmTMB:::fitTMB(tmb_struc)
  rm(tmb_struc)
  saveRDS(TMB, paste0("/home/hbronnvik/Documents/chapter3/model_fits/", names(f)[l], ".rds"))
  gc()
  return(TMB)
})
# Warning messages:
#   1: In finalizeTMB(TMBStruc, obj, fit, h, data.tmb.old) :
#   Model convergence problem; false convergence (8). See vignette('troubleshooting'), help('diagnose')
# 2: In finalizeTMB(TMBStruc, obj, fit, h, data.tmb.old) :
#   Model convergence problem; false convergence (8). See vignette('troubleshooting'), help('diagnose')
tmbs <- lapply(list.files("/home/hbronnvik/Documents/chapter3/model_fits/with_bins", full.names = T), readRDS)
names(tmbs) <- sub(".rds", "", list.files("/home/hbronnvik/Documents/chapter3/model_fits/with_bins/"))
# save out to push to the cluster
saveRDS(tmbs, file = "/home/hbronnvik/Documents/chapter3/tmb_fits_bins7x.rds")

# a subset of the model fits
# fls <- data.frame(f = list.files("/home/hbronnvik/Documents/chapter3/model_fits",
#                                  full.names = T)) %>%
#   mutate(fn = sub(".rds", "", list.files("/home/hbronnvik/Documents/chapter3/model_fits"))) %>%
#   filter(!grepl("asp", fn))
# mods_static <- lapply(fls$f, readRDS)
# names(mods_static) <- fls$fn
# saveRDS(mods_static, file = "/home/hbronnvik/Documents/chapter3/tmb_fits_static.rds")

performances <- as.data.frame(performance(tmbs[[1]])) %>%
  rbind(performance(tmbs[[2]])) %>%
  rbind(performance(tmbs[[3]])) %>%
  rbind(performance(tmbs[[4]])) %>%
  rbind(performance(tmbs[[5]])) %>%
  rbind(performance(tmbs[[6]])) %>% 
  rbind(performance(tmbs[[7]])) %>% 
  mutate(model = names(tmbs)) %>% 
  arrange(RMSE)

ests <- lapply(1:length(tmbs), function(x){
  e <- confint(tmbs[[x]]) %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "pred")  %>% 
    filter(!grepl("id|sl|ta", pred)) %>% 
    rename(lower = "2.5 %",
           upper = "97.5 %") %>% 
    mutate(model = names(tmbs)[x])
  return(e)
}) %>% reduce(rbind)

ests %>% 
  mutate(pred = sub("sqrt_dem_z", "DEM",
                    sub("log_blh_z", "BLH", 
                        sub("grass_crop_z", "Grass/crop",
                            pred)))) %>% 
  ggplot(aes(Estimate, pred)) +
  geom_vline(xintercept = 0, lty = 2, linewidth = 1) +
  geom_pointrange(aes(xmin = lower, xmax = upper),
                  linewidth = 1, size = 1, color = "#0081A7") +
  facet_wrap(~model)

# png(paste0("/home/hbronnvik/Documents/chapter3/figures/fall_ssf_coefs.png"),
#     width = 8.2, height = 11.7, units = "in", res = 300)
ggpubr::ggarrange(coef_plots[[1]], coef_plots[[2]], coef_plots[[3]], nrow = 3)
# dev.off()

## 4. predict from the models across the flyway
# get a file of all lat/longs and DEMs for a raster 
pred_z_ls <- list.files("/home/hbronnvik/Documents/chapter3/pred_tiles5x_z", full.names = T)
IDs <- readRDS("/home/hbronnvik/Documents/chapter3/adata_ids.rds")

start_time <- Sys.time() # grep("asp", names(tmbs))
mod_preds <- lapply(c(1, 2, 4, 6, 7), function(x){
  # predict from one model at a time and clean up
  TMB_mod <- tmbs[[x]]
  
  print(paste0("Predicting selection probabilities for model: ", names(tmbs)[x], "."), quote = F)
  # go through each 500 km tile and predict using that model
  pred_tiles <- lapply(pred_z_ls, function(p){
    # print(p)
    template <- rast(p)
    names(template)[grepl("forage", names(template))] <- "forage_z"
    names(template)[grepl("locs", names(template))] <- "obs_z"
    # generate a new dataset for predictions
    lb <- template %>% 
      as.data.frame(xy = T, na.rm = F, cell = T) %>% 
      # rename(grass_crop_z = cropland) %>% 
      mutate(sqrt_sl_z = 0,
             ta_z = 0,
             #set the grouping variables to NA as per glmm.TMB recommendation
             stratum_ID = NA,
             individual.id = sample(IDs, n(), replace = T)) 
    # predict using the models
    lb$preds <- predict(TMB_mod, newdata = lb, type = "link")
    gc()
    # convert the log odds to probabilities
    lb$probs <- gtools::inv.logit(lb$preds)
    gc()
    # rasterize these probabilities
    lb_probs <- vect(lb, geom = c("x", "y"), crs = crs(template))
    lb_probs <- rast(x = lb_probs, extent = ext(template), resolution = res(template), 
                     crs = crs(template), names = names(tmbs)[x], 
                     vals = c(lb_probs$probs))
    # plot(lb_probs)
    return(lb_probs)
  })
  print("Merging model predictions.", quote = F)
  # merge the 500 km tiles into the flyway and return
  pred_tiles <- sprc(pred_tiles)
  pred_tiles <- merge(pred_tiles)
  saveRDS(pred_tiles, file = paste0("/home/hbronnvik/Documents/chapter3/ssf_mod_preds/west_", names(tmbs)[x], ".rds"))
  rm(TMB_mod);gc()
  return(pred_tiles)
})
Sys.time()-start_time # Time difference of 5.335511 hours (3 mods)
pred_files <- list.files("/home/hbronnvik/Documents/chapter3/ssf_mod_preds", pattern = "west", full.names = T)
mod_preds <- lapply(pred_files, rast)
pred_tiles <- c(mod_preds[[1]], mod_preds[[2]], mod_preds[[3]], mod_preds[[4]], mod_preds[[5]])
# pred_tiles <- merge(pred_tiles)
plot(pred_tiles)

## 6. convert to resistance layers and export
# defining resistance to movement as the opposite of probability of use
# the files off the computing cluster
preds <- list.files("/home/hbronnvik/Documents/chapter3/Raven", pattern = "prediction_rasts",
                    full.names = T)
resistances <- lapply(preds, function(x){
  res <- 1-rast(x)
  writeRaster(res, file = paste0("/home/hbronnvik/Documents/chapter3/Circuit/resist_", names(res),"_fall_strait.tif"), overwrite = T)
})
# writeRaster(res[[1]], file = paste0("/home/hbronnvik/Documents/chapter3/Circuit/resist_", names(res[[1]]),"_fall_strait.tif"), overwrite = T)
# writeRaster(res[[2]], file = paste0("/home/hbronnvik/Documents/chapter3/Circuit/resist_", names(res[[2]]),"_fall_strait.tif"), overwrite = T)
# writeRaster(res[[3]], file = paste0("/home/hbronnvik/Documents/chapter3/Circuit/resist_", names(res[[3]]),"_fall_strait.tif"), overwrite = T)
# writeRaster(res[[4]], file = paste0("/home/hbronnvik/Documents/chapter3/Circuit/resist_", names(res[[4]]),"_fall_strait.tif"), overwrite = T)
# writeRaster(res[[5]], file = paste0("/home/hbronnvik/Documents/chapter3/Circuit/resist_", names(res[[5]]),"_fall_strait.tif"), overwrite = T)

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>% 
  sf::st_transform("ESRI:53009")

colfunc2 <- colorRampPalette(c("#dedde4ff", "#e2c5d9ff", "#e498c2ff", "#e387b7ff", "#d9608fff", 
                               "#d45078ff", "#df7070ff", "#ea9565ff", "#efa55eff", "#f3b457ff", 
                               "#F5C050", "#F7C74C", "#F8CB49"))
# names(res) <- c("Elevation (m)", "Elevation and boundary layer height", "Elevation, B.L.H., & Greenness")
# res <- project(res, "EPSG:4326")
# png("/home/hbronnvik/Documents/chapter3/figures/resistance_layers.png",
#     height = 8.2, width = 11.7, units = "in", res = 400)
ggplot() +
  tidyterra::geom_spatraster(data = res) +
  scale_fill_gradientn("Relative\nresistance", colors = colfunc2(200), n.breaks = 7, na.value = "white") +
  # tidyterra::scale_fill_princess_c(palette = "maori", n.breaks = 13, na.value = "white", n = 150,
  #                                  guide = guide_legend(reverse = TRUE, title = "Relative\nresistance")) +
  geom_sf(data = rnaturalearth::ne_coastline(scale = "medium", returnclass = "sf") %>%
            sf::st_transform("ESRI:53009"), color = "grey40") +
  coord_sf(xlim = c(-1997508.40406021, 1697882.14345118), ylim = c(617232.175201235, 6423128.52866402)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_wrap(~lyr) +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())
# dev.off()


# off the cluster
preds <- readRDS("/home/hbronnvik/Documents/chapter3/ssf_preds_dfs_5x.rds")
names(preds)
head(preds[[3]])
summary(preds[[1]]$y)
preds[[1]]$probs <- gtools::inv.logit(preds[[1]]$preds)
ggplot2::ggplot(preds[[1]][preds[[1]]$y > 4971525,], ggplot2::aes(x, y, color = probs)) +
  ggplot2::geom_tile()


