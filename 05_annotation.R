### annotate burst data with terrain variables

library(terra)
library(sf)
library(geodata)
library(tidyverse)
theme_set(theme_classic()+
            theme(axis.text = element_text(color = "black", size = 12), 
                  text = element_text(size = 15)))
colfunc <- colorRampPalette(c("#0081A7", "#0098b0", "#00AFB9", "#7fc4b8", "#FED9B7", "#f7a58f", "#f27e71", "#F07167", "#ee5e53"))

# 1. the location data from 03_step_generation
ua_locs <- readRDS("/home/hbronnvik/Documents/chapter3/used_av_df_60min.rds") %>% 
  as.data.frame()
table(ua_locs$used)

# 2. rasters prepared in 04_env_rasts.R
lyrs <- rast("/home/hbronnvik/Documents/chapter3/env_rasts/full_terr.tif")
# the ODs from stork_ods.R
ods <-rast("/home/hbronnvik/Documents/chapter3/ctmm_ods_sum_fall1.tif")
uds <- rast("/home/hbronnvik/Documents/chapter3/ctmm_uds_sum_fall1.tif")
lyrs <- c(lyrs, ods, uds)

# pare down and reclassify the layers
lyrs <- lyrs[[c("DEM", "obs", "max_blh")]]

# obs <- lyrs[["obs"]]
# obs <- classify(obs, cbind(c(-1,1,5,10), c(1,5,10,140), c(1,2,3,4)))
# levels(obs) <- data.frame(ID = 1:4, Category = c("none", "few", "some", "many"))
# plot(obs)
# names(obs) <- "obs"

# lyrs <- c(lyrs[[c("DEM", "max_blh")]], obs)

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

# 3. annotate
ua_locs_ls <- ua_locs %>% 
  group_by(individual.id) %>% 
  group_split()
# too many locations, R will crash, so we do it in parts
ua_locs_anno <- lapply(ua_locs_ls, function(ua){
  # make a spatial object
  ua <- ua %>% 
    st_as_sf(coords = c("location.long", "location.lat"), crs = 4326) %>% 
    st_transform(crs = "ESRI:53009") 
  # extracts
  ua <- ua %>% 
    mutate(terra::extract(lyrs, ua))
  # additions
  ua <- ua %>% 
    mutate(obs = ifelse(obs >= 1, obs-1, obs),# remove the individual itself from each cell
           # add on the long/lat again (in m) and drop geometry for file size concerns
           location.long = st_coordinates(.)[,1],
           location.lat = st_coordinates(.)[,2]) %>% 
    st_drop_geometry()
  gc()
  return(ua)
})

ua_locs <- data.table::rbindlist(ua_locs_anno)

# discretize the observations
# create bins with no (0-1], few (1-5], intermediate number (5-10], or many (> 10)
ua_locs <- ua_locs %>% 
  mutate(obs_bin = ifelse(between(obs, -1, 0), "none",
                          ifelse(between(obs, 1, 4), "few",
                                 ifelse(between(obs, 5, 9), "some",
                                        ifelse(between(obs, 10, 200), "many", obs)))),
         obs_bin = factor(obs_bin, levels = c("none", "few", "some", "many")))

# ua_locs_ls <- ua_locs %>% 
#   mutate(stamp = round_date(timestamp, unit = "hour")) %>%
#   group_by(year(timestamp)) %>% 
#   group_split()

# the boundary layer heights
# pbls <- list.files("/home/hbronnvik/Documents/chapter3/hourly_blhs", full.names = T)

# start_time <- Sys.time()
# anno <- lapply(1:length(ua_locs_ls), function(s){
#   
#   steps <- ua_locs_ls[[s]]
#   # identify the year, month, and hours
#   d_year <- unique(year(steps$timestamp))
#   stamp <- round_date(unique(steps$timestamp), unit = "hour")
#   print(d_year)
#   # print(paste0("Extracting environmental data for ", month.name[d_month], " ", d_year, ". Item ",
#   #              s, " of ", length(data_ls), "."), 
#   #       quote = F)
#   
#   # the height of the boundary layer
#   file <- grep(d_year, pbls, value = T)
#   blh <- try(terra::rast(file))
#   # the time stamps from Copernicus don't match terra's format. Extract and add them.
#   ts <- as.POSIXct(depth(blh), tz = "UTC")
#   time(blh) <- ts
#   # then subset the raster stack and keep only the hours we need
#   blh <- blh[[which(terra::time(blh) %in% stamp)]]
#   # create a list to go through each hour individually
#   steps_ls <- steps %>% 
#     st_as_sf(coords = c("location.long", "location.lat"), crs = "ESRI:53009") %>% 
#     st_transform(crs = crs(blh)) %>% 
#     group_by(stamp) %>% 
#     group_split()
#   
#   # the boundary layer height at the points of interest (bi-linearly interpolated)
#   steps_ls <- lapply(steps_ls, function(df){
#     ts <- unique(df$stamp)
#     h <- blh[[which(terra::time(blh) == ts)]]
#     df$blh_exact <- terra::extract(h, df, method = "bilinear")[,2]
#     df %>% 
#       mutate(location.long = st_coordinates(.)[,1],
#              location.lat = st_coordinates(.)[,2]) %>% 
#       st_drop_geometry()
#   })
#   steps <- data.table::rbindlist(steps_ls)
#   return(steps)
# })
# Sys.time()-start_time # Time difference of 2.851933 hours
# ua_locs <- data.table::rbindlist(anno)
# 
# rm(ua_locs_anno, ua_locs_ls, anno);gc()

# quickly look at the relationship between use and these variables
# png(paste0("/home/hbronnvik/Documents/chapter3/figures/land_var_dists.png"),
#     height = 8.2, width = 11.7, units = "in", res = 300, bg = "transparent")
ua_locs %>% 
  dplyr::select(used, obs, max_blh, DEM) %>% 
  pivot_longer(cols = !used) %>% 
  ggplot(aes(as.factor(used), value, fill = as.factor(used))) +
  geom_boxplot() +
  labs(x = "Used", y = "") +
  facet_wrap(~name, scales= "free") + 
  theme(legend.position = "none")
# dev.off()

ms <- ua_locs %>% 
  filter(used == T) %>% 
  group_by(compass) %>% 
  summarize(m = mean(sqrt(DEM), na.rm = T))
# png(paste0("/home/hbronnvik/Documents/chapter3/figures/demXcompass_used.png"),
#     height = 8.2, width = 11.7, units = "in", res = 300, bg = "transparent")
ua_locs %>% 
  filter(used == 1) %>% 
  left_join(ms) %>% 
  drop_na(aspect) %>% 
  ggplot(aes(compass, sqrt(DEM), fill = m)) +
  geom_boxplot() +
  labs(x = "Compass direction of the slope", y = "sqrt DEM", fill = "Mean sqrt DEM")
# dev.off()
ua_locs %>% 
  filter(used == 1) %>% 
  left_join(ms) %>% 
  drop_na(DEM) %>% 
  ggplot(aes(compass, sqrt(slope))) +
  geom_boxplot() +
  labs(x = "Compass direction of the slope", y = "sqrt slope", fill = "Mean sqrt DEM")

# check correlations (max. is blhXdem = 0.333)
ua_locs %>% 
  dplyr::select(used, DEM, obs, aspect, slope, roughness, TRI, TPI) %>% 
  drop_na(DEM) %>% 
  drop_na(obs) %>% 
  drop_na(slope) %>% 
  drop_na(aspect) %>% 
  drop_na(roughness) %>% 
  drop_na(TRI) %>% 
  drop_na(TPI) %>% 
  cor()

                  # used          DEM          obs       aspect        slope    roughness
# used       1.000000000 -0.035248343  0.016591187 -0.012637737 -0.036655975 -0.038765634
# DEM       -0.035248343  1.000000000 -0.208402816  0.003243313  0.468283494  0.505757961
# obs        0.016591187 -0.208402816  1.000000000  0.003104891  0.005721162  0.004210816
# aspect    -0.012637737  0.003243313  0.003104891  1.000000000 -0.006436424 -0.009272129
# slope     -0.036655975  0.468283494  0.005721162 -0.006436424  1.000000000  0.969962471
# roughness -0.038765634  0.505757961  0.004210816 -0.009272129  0.969962471  1.000000000
# TRI       -0.039110922  0.508620503  0.001895856 -0.012196772  0.948034578  0.975302622
# TPI       -0.006570716  0.155844372 -0.002646889 -0.004541640  0.047589022  0.047932573
                   # TRI          TPI
# used      -0.039110922 -0.006570716
# DEM        0.508620503  0.155844372
# obs        0.001895856 -0.002646889
# aspect    -0.012196772 -0.004541640
# slope      0.948034578  0.047589022
# roughness  0.975302622  0.047932573
# TRI        1.000000000  0.095185494
# TPI        0.095185494  1.000000000

# We will not include slope, roughness, or TRI because of their correlation to DEM, TPI requires binning

# normalize and scale the predictors
hist(ua_locs$DEM) 
# skew the data to account for the Dead Sea
hist(log(ua_locs$DEM+424), main = "", xlab = "log DEM")

hist(ua_locs$step_length)
hist(sqrt(ua_locs$step_length), main = "", xlab = "sqrt sl_")

# ggplot(ua_locs, aes(turning_angle, fill = as.factor(used))) +
#   geom_density(color = "black", alpha = 0.5) +
#   labs(x = "Turning angle (deg)", y = "Density", fill = "Use")

# ua_locs <- ua_locs %>% 
#   # fix the turning angles
#   mutate(turning_angle = ifelse(used == 0, (turning_angle + 180)%%360-180, turning_angle))
# 
# ggplot(ua_locs, aes(turning_angle, fill = as.factor(used))) +
#   geom_density(color = "black", alpha = 0.5) +
#   labs(x = "Turning angle (deg)", y = "Density", fill = "Use")

# ggplot(ua_locs, aes(as.factor(hour(timestamp)), aspect, fill = as.factor(used))) +
#   geom_hline(yintercept = c(0, 90, 180, 270, 360)) +
#   geom_boxplot()

hist(ua_locs$turning_angle)

summary(ua_locs$obs)
hist(ua_locs$obs)
hist(1/(1+ua_locs$obs)**2, main = "", xlab = "1/obs^2")

summary(ua_locs$UDs)
hist(ua_locs$UDs)
hist(sqrt(ua_locs$UDs))

# transform and scale the predictors for the models
ua_locs <- ua_locs %>% 
  mutate(sqrt_dem_z = scale(log(DEM+424))[,1], 
         # aspect_z = scale(aspect)[,1],
         obs_num_z = scale(1/(1+obs)**2)[,1], 
         # ods_z = scale(ODs)[,1], 
         # sqrt_uds_z = scale(sqrt(UDs))[,1],
         # med_blh_z = scale(med_blh)[,1], 
         max_blh_z = scale(1/max_blh)[,1],
         # pnt_blh_z = scale(blh_exact)[,1],
         sqrt_sl_z = scale(sqrt(step_length))[,1],
         ta_z = scale(turning_angle)[,1]) %>% 
  rename(obs_z = obs_bin)

# saveRDS(ua_locs, file = "/home/hbronnvik/Documents/chapter3/ua_data_60min_20250804.rds")


# get the layers as a data frame to predict on
filled_lyrs <- as.data.frame(lyrs, na.rm = F, xy = T)
# ensure they are scaled and transformed the same way the model data were
att <- ua_locs %>% 
  # the attributes of the annotated data to scale the full data with
  summarize(m_dem = mean(log(DEM+424), na.rm = T),
            sd_dem = sd(log(DEM+424), na.rm = T),
            m_blh = mean((1/max_blh), na.rm = T),
            sd_blh = sd((1/max_blh), na.rm = T))

# scale the data
filled_lyrs <- filled_lyrs %>% 
  mutate(sqrt_dem_z = (log(DEM+424)-att$m_dem)/att$sd_dem,
         max_blh_z = ((1/max_blh)-att$m_blh)/att$sd_blh)

# saveRDS(filled_lyrs, file = "/home/hbronnvik/Documents/chapter3/filled_lyrs_df.rds")










