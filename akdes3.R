### Estimate population distributions and model
### Hester Bronnvik
### 2023-07-13
### hbronnvik@ab.mpg.de

library(move)
library(ctmm)
library(lubridate)
library(raster)
library(terra)
library(corrr)
library(performance)
library(glmmTMB)
library(data.table)
library(cowplot)
library(tidyverse)

# functions: 
wind_support <- function(u,v,heading) {
  angle <- atan2(u,v) - heading/180*pi
  return(cos(angle) * sqrt(u*u+v*v))
}
cross_wind <- function(u,v,heading) {
  angle <- atan2(u,v) - heading/180*pi
  return(sin(angle) * sqrt(u*u+v*v))
}
wind_speed <- function(u,v) {
  return(sqrt(u*u+v*v))
}
CubeRoot<-function(x){
  sign(x)*abs(x)^(1/3)
}
convective <- function(data) {
  
  df <- data # the data input
  z <- df$blh # the boundary layer height
  phi <- df$phi # the geopotential at the layer below the boundary layer
  s_flux <- df$sflux # the surface sensible heat flux
  m_flux <- df$mflux # the moisture flux
  T_k_2m <- df$T2m # the temperature at ground level
  p0 <- df$p0 # the pressure at ground level (Pascals)
  p0 <- p0/100  # the pressure at ground level (mbar)
  q <- df$q # the humidity below the boundary layer height
  p1 <- as.numeric(df$level) # the pressure below the boundary layer height (mbar)
  T_K_blh <- df$temp # the temperature below the boundary layer height
  
  k <- 0.2854 # Poisson constant
  g <- 9.80665 # acceleration due to gravity
  T_c_2m <- T_k_2m - 273.15 # surface temperature in degrees Celsius
  c_p <- 1012 # the isobaric mass heat capacity of air at common conditions
  p <- 1.225 # the density of air at sea level and 273 K
  
  Theta_k_z <- T_K_blh*((p0/p1)^k)
  Thetav_k_z <- Theta_k_z*(1+(0.61*q))
  wT <- (s_flux * -1)/c_p/p # reverse the sign. ECMWF upward fluxes are negative
  wq <- (m_flux * -1)/p # reverse the sign
  
  wthetav <- wT + 0.61 * T_c_2m * wq
  
  w_star <- CubeRoot(phi*(wthetav/Thetav_k_z))
  
  return(w_star)
  
}# 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

## 1) Get the data

# read in the saved data for the 509 birds for all locations during migration
locs <- readRDS("/home/hbronnvik/Documents/storkSSFs/soc_info_migration_locations_2023_06_12.rds")
locs_add <- readRDS("/home/hbronnvik/Documents/storkSSFs/soc_info_migration_locations_Sarralbe_2023_08_21.rds") %>% 
  mutate("manually.marked.valid" = NA)
locs <- locs %>% 
  rbind(locs_add)
# the total number of migrations attempted and completed by each animal in each season 
meta <- locs %>%
  mutate(season = ifelse(grepl("fall", trackID), "post", "pre")) %>% 
  group_by(individual.id, season) %>% 
  count(trackID) %>% 
  mutate(journey_number = row_number()) %>% 
  ungroup() %>% 
  group_by(individual.id) %>%
  mutate(total_journeys = n()) %>% 
  ungroup() 
# add the number of journeys
locs <- locs %>% 
  rowwise() %>% 
  mutate(journey_number = meta$journey_number[which(meta$trackID == trackID)]) %>% 
  ungroup()
# align these locations in time
locs$datestamp <- locs$timestamp
year(locs$datestamp) <- 2024
second(locs$datestamp) <- 00
locs$datestamp <- round_date(locs$datestamp, "hour")

# read in the saved data for the 158 birds that completed migrations
# these have already been burst (see step_generation.R) and annotated with uplift (see wstar.R)
a_data <- readRDS("/home/hbronnvik/Documents/storkSSFs/full_annotated_230623.rds")
# align these locations in time
a_data$alignment <- round_date(a_data$timestamp, unit = "hour")
year(a_data$alignment) <- 2024
second(a_data$alignment) <- 00

# get the outlines of land to remove the possibility of information over water
# tmax <- raster::getData('worldclim', var = "tmax", res = 10)
# mask <- crop(raster(tmax, 1), raster::extent(min(first_locs$location.long), max(first_locs$location.long), min(first_locs$location.lat), max(first_locs$location.lat)))
# raster::plot(mask)
# # the Straits of Gibraltar are at most 60km across
# mm <- buffer(mask, 30000)
# raster::plot(mm)
# writeRaster(mm, "C:/Users/hbronnvik/Documents/storkSSFs/buffered_water_ras.grd")

# call in a buffered raster of land outlines + 30km
mm <- raster::raster("/home/hbronnvik/Documents/storkSSFs/buffered_water_ras.grd")
outlines <- raster::rasterToPolygons(mm, dissolve=TRUE)

# now, take one of the 158 animals at a time
# take it out of the 367 birds that migrated, and use the rest to build one AKDE per hour
# run the one bird across those AKDEs, annotating its locations with the AKDE for that hour
# lather, rinse, repeat 
start_time <- Sys.time()
HR <- lapply(unique(a_data$individual.id), function(x){
  # isolate the locations of a single, focal bird
  ind <- a_data %>% 
    filter(individual.id == x)
  # take all the locations that are not from that bird
  soc <- locs %>% 
    filter(individual.id != x) %>% 
    mutate(individual.id = as.character(datestamp),
           individual.local.identifier = as.character(datestamp)) %>% 
    group_by(individual.id) %>% 
    mutate(count = n()) %>% 
    ungroup() %>% 
    # only use hours with more than 5 locations
    filter(count >= 5)%>% 
    filter(datestamp %in% unique(a_data$alignment))
  # compress into a leap year
  year(soc$timestamp) <- 2024
  
  # take all of the migratory locations within the hours of migrations by the focal bird
  soc <- soc %>% 
    filter(datestamp %in% unique(ind$alignment))
  # make a list of telemetry objects
  soc_tele <- suppressWarnings(suppressMessages(as.telemetry(soc, projection = "ESRI:54009")))
  # make a list of model fits
  fits <- lapply(soc_tele, function(x){
    tryCatch({
      guess <- ctmm.guess(x, interactive=FALSE)
      # fit the models
      fit <- ctmm.select(x, guess)
      return(fit)
    }, error = function(e){print(geterrmessage(), quote = F)})
  }) %>% suppressWarnings() # warning: Duplicate timestamps require an error model.
  # (we ignore this warning because we are deliberately using only one timestamp and a population not an individual)
  # estimate utilization distributions of all (the non-focal) storks in each hour 
  # cropped to the buffered map, and on a grid of set extent and 30km resolution (approximately matching the weather data)
  kdes <- akde(soc_tele, fits, SP = outlines, grid = list(dr = 30000, extent = extent(-2468787, 3006683, 0, 6876759)))
  # go through the KDEs to annotate the tracks
  akde_data <- lapply(kdes, function(y){
    # the hour covered by this KDE
    kde_hour <- y@info$identity
    # the locations in that hour
    tracks <- ind %>% 
      filter(as.character(alignment) == kde_hour)
    if(nrow(tracks) > 0){
      # convert to a spatial object
      coordinates(tracks) <- ~ long + lat
      proj4string(tracks) <- CRS("EPSG:4326")
      tracks <- spTransform(tracks, "ESRI:54009")
      # the KDE as a raster ("PDF" gives the average probability density per cell)
      ud <- raster(y, DF = "PDF")
      # extract the UD value
      vals <- raster::extract(ud, tracks)
      # append
      df <- ind %>% 
        filter(as.character(alignment) == kde_hour) %>% 
        mutate(UD_PDF = vals)
      return(df)
    }
  }) %>% 
    discard(is.null)
  akde_data <- rbindlist(akde_data)
  print(paste0("Annotated data for individual ", which(unique(a_data$individual.id) == x), " of ", 
               length(unique(a_data$individual.id)), "."))
  return(akde_data)
}) %>% rbindlist()
Sys.time() - start_time # Time difference of 18.4626 hours

a_data <- a_data %>% 
  full_join(HR)

# saveRDS(a_data, file = "/home/hbronnvik/Documents/storkSSFs/annotations/HR_230713.rds")

a_data <- readRDS("/home/hbronnvik/Documents/storkSSFs/annotations/HR_230713.rds")

summary(a_data$UD_PDF)

# estimate w* for the full burst data
# first, determine phi (geopotential below the top of the boundary layer)
a_data$phi <- lapply(1:nrow(a_data), function(x){
  l <- a_data$level[x]
  geoH <- a_data[x, which(colnames(a_data) == l)]
  phi <- geoH*9.80665
}) %>% unlist()
# then, run the function (assuming all column names match)
a_data$w_star <- convective(a_data)

summary(a_data$w_star);hist(a_data$w_star, breaks = 100, main = "", xlab = "Convective velocity scale")

fac_labs <- c("Fall", "Spring")
names(fac_labs) <- c("post", "pre")
ggplot(a_data, aes(x = as.factor(migrations), y = w_star, fill = as.factor(used))) +
  geom_boxplot() +
  scale_fill_manual(values = c("#0096CC", "#FE6D5D")) +
  labs(x = "Migrations", y = "Uplift (m/s)", fill = "Use") +
  theme_classic() +
  facet_wrap(~season, labeller = labeller(season = fac_labs))

# we still need one more variable -- wind
# determine which level of wind to use by looking at the overall flight heights of white storks
# all of the full, raw data downloaded from Movebank
files <- list.files("/home/hbronnvik/Documents/storkSSFs/full_data", pattern = ".rds", full.names = T)
# the geoid heights from EGM2008 https://www.agisoft.com/downloads/geoids/
egm <- terra::rast("/home/hbronnvik/Documents/storkSSFs/us_nga_egm2008_1.tif")

# take the geoid height (AKA geoid undulation) out of the height above ellipsoid to get height above geoid (AKA msl)
# start with a simple collection of all the heights above ellipsoid
# there can be a lot of error in these readings, so we remove the very large or very small values (although error exists in the plausible ones as well)
heights <- lapply(files, function(file){
  data <- readRDS(file)
  data <- data %>%
    # take out the bursts
    mutate(td = as.numeric(difftime(timestamp, lag(timestamp), units = "secs"))) %>%
    filter(td >= 300) %>%
    mutate(seq15 = round_date(timestamp, unit = "15 minutes")) %>%
    group_by(seq15) %>%
    slice(1) %>%
    ungroup() %>%
    dplyr::select(-seq15, -td) %>% 
    mutate(distance = distVincentyEllipsoid(cbind(location.long, location.lat), cbind(lag(location.long), lag(location.lat))),
           timediff = as.numeric(difftime(timestamp, lag(timestamp), units = "secs")),
           ground_speed_15 = distance/timediff) %>% 
    # use only locations that are plausible and in-flight
    filter(ground_speed_15 < 50 & ground_speed_15 > 2)
  hae <- data %>% 
    filter(height.above.ellipsoid > -100 & height.above.ellipsoid < 10000) %>% 
    dplyr::select(timestamp, location.long, location.lat, height.above.ellipsoid)
  hae
})
heights <- data.table::rbindlist(heights)

# next, extract the geoid undulations at each of the locations that has HAE
heights$geoid_height <- terra::extract(egm, vect(heights, geom = c("location.long", "location.lat")))[ ,2] 

# add on the height above mean sea level
heights$ha_msl <- heights$height.above.ellipsoid - heights$geoid_height

# find the average flight height above mean sea level, 
# but try to use only the trustworthy ones by not allowing the birds to fly above the boundary layer
heights <- heights %>% 
  filter(ha_msl < max(na.omit(a_data$blh)))

mean(heights$ha_msl)

# compare this mean to the mean heights of each pressure level in the annotations
log_heights <- lapply(15:28, function(x){
  level_h <- a_data[, x]
  level <- colnames(a_data)[x]
  info <- data.frame(level = level, mean = mean(na.omit(log(level_h))), sd = sd(na.omit(log(level_h))))
}) %>% reduce(rbind) %>% arrange(as.numeric(level))

log_heights
mean(na.omit(log(heights$ha_msl)))

# read out files to annotate with the wind data from Movebank Env-DATA service at the given level
# https://www.movebank.org/cms/movebank-content/env-data
wind_file <- a_data %>% 
  dplyr::select(timestamp, long, lat, individual.id) %>% 
  mutate(timestamp = paste0(timestamp, ".000"),
         group = c(rep(c("a"), times = n()/4),
                   rep(c("b"), times = n()/4),
                   rep(c("c"), times = n()/4),
                   rep(c("d"), times = n()/4))) %>% 
  rename("location-long" = long,
         "location-lat" = lat)

# lapply(split(wind_file, wind_file$group), function(s){
#   write.csv(s, file = paste0("/home/hbronnvik/Documents/storkSSFs/wind_to_fill_", s$group[1],".csv"))
# })

# after annotation, read in the files from Movebank Env-DATA
wind_files <- list.files("/home/hbronnvik/Documents/storkSSFs/ecmwf/winds", pattern = ".csv", full.names = T)

# make the column names convenient
winds <- lapply(wind_files, read.csv) %>% reduce(rbind)
winds$X <- NULL
colnames(winds)
colnames(winds)[2:3] <- c("long", "lat")
colnames(winds)[6:11] <- c("u_950", "v_950", "u_10m", "u_100m", "v_100m", "v_10m")
# turn the stamp back to an R format
winds <- winds %>% 
  mutate(timestamp = as.POSIXct(timestamp, tz = "UTC"))

# also make the long/lat format match the a_data so that the join will recognize them (this is why I don't work in .csv)
a_data <- a_data %>% 
  mutate(long = as.character(long),
         lat = as.character(lat)) %>% 
  mutate(long = as.numeric(long),
         lat = as.numeric(lat))

# add the wind data to the steps with their social density and w*
a_data <- full_join(a_data, winds) %>% 
  mutate(cross_wind = cross_wind(u_950, v_950, heading),
         wind_support = wind_support(u_950, v_950, heading),
         wind_speed = wind_speed(u_950, v_950),
         wind_speed_100m = wind_speed(u_100m, v_100m),
         wind_speed_ground = wind_speed(u_10m, v_10m))

# png(filename = "/home/hbronnvik/Documents/storkSSFs/figures/wind_correlations.png", res = 1000, 
#     width = 11.69, height = 8.27, units = "in")
# print(ggplot(a_data, aes(wind_speed_ground, wind_speed)) +
#         geom_point(alpha = 0.1)+
#         geom_point(data = a_data, aes(wind_speed_ground, wind_speed_100m), pch = 15, alpha = 0.1)+
#         geom_smooth(method = "lm", color = "red") +
#         geom_smooth(data = a_data, aes(wind_speed_ground, wind_speed_100m), method = "lm", color = "orange") +
#         labs(x = "Wind speed at 10m above ground", y = "Wind speed at pressure") +
#         theme_classic()
# )
# dev.off()

meta <- meta %>% 
  rename(track = trackID) %>% 
  dplyr::select(individual.id, track, journey_number, total_journeys)

# add the number of journeys
a_data <- a_data %>% 
  arrange(timestamp) %>% 
  left_join(meta)

# scale the predictors
a_data <- a_data %>% 
  rename(ud_pdf = UD_PDF,
         migrations = journey_number) %>% 
  mutate_at(c("migrations", "wind_support", "cross_wind", "w_star", "wind_speed", "ud_pdf", "step_length", "turning_angle", "blh"),
            list(z = ~(scale(.)))) %>% 
  # remove the geopotential heights
  dplyr::select(-colnames(.)[15:28])

# square root transform the social influence proxy for normality
a_data <- a_data %>% 
  mutate(sqrt_ud = sqrt(ud_pdf),
         sqrt_ud_z = scale(sqrt_ud))

# add on whether locations are over sea or over land using the buffered raster we used for the AKDEs
spatMM <- terra::rast(mm)
poi_vals <- terra::extract(spatMM, vect(a_data, geom = c("long", "lat")))[,2]
unique(poi_vals)
a_data$poi_vals <- poi_vals
a_data$land <- ifelse(is.na(poi_vals), "water", "land")

# ex <- a_data %>%
#   filter(used == 1)%>%
#   select(long, lat, stratum, land)
# ex <- st_as_sf(ex, coords=1:2, crs=4326)
# mapview::mapview(ex, zcol = "land")

# there are no actual, observed points classified as being over water
table(is.na(a_data$poi_vals), a_data$used)

a_data %>% 
  filter(is.na(poi_vals)) %>% 
  group_by(stratum) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  select(n, stratum, used, step_length) %>% 
  arrange(desc(n))

# remove the "available" locations that we generated over the sea and pare down the data
# so that all strata are the same size (not strictly needed for Bayesian approaches)
used_ad <- a_data %>% 
  filter(used == 1)
avail_ad <- a_data %>% 
  filter(used == 0 & land == "land") %>% 
  group_by(stratum) %>% 
  slice_sample(n = 49) %>% 
  ungroup()

a_data <- used_ad %>% 
  rbind(avail_ad) %>% 
  arrange(individual.id, stratum)

a_data %>% 
  # filter(land == F) %>% 
  group_by(stratum) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  arrange(desc(n)) %>% 
  select(n) %>% 
  unique()

# save out the data 
saveRDS(a_data, file = paste0("/home/hbronnvik/Documents/storkSSFs/a_data_", Sys.Date(), ".rds"))
a_data <- a_data %>% 
  mutate(season = ifelse(grepl("fall", track), "post", "pre")) %>% 
  select(timestamp, long, lat, used, individual.id, stratum, migrations_z, 
         wind_speed_z, w_star_z, sqrt_ud_z, step_length_z, turning_angle_z, season) 
saveRDS(a_data, file = paste0("/home/hbronnvik/Documents/storkSSFs/a_data_clust_", Sys.Date(), ".rds"))

a_data <- readRDS("/home/hbronnvik/Desktop/cluster_INLA/a_data_clust_2023-07-26.rds")

a_data <- readRDS("/home/hbronnvik/Documents/storkSSFs/a_data_2023-07-26.rds")%>% 
  mutate(season = ifelse(grepl("fall", track), "post", "pre"))

colfunc <- colorRampPalette(c("#0081A7", "#0098b0", "#00AFB9", "#7fc4b8", "#FED9B7", "#f7a58f", "#f27e71", "#F07167", "#ee5e53"))

a_data$datestamp <- a_data$timestamp
year(a_data$datestamp) <- 2024

png(filename = "/home/hbronnvik/Documents/storkSSFs/figures/migration_timing.png",
    height = 18, width = 20, units = "in", res = 300)
a_data %>% 
  mutate(yd = date(datestamp)) %>% 
  group_by(yd, migrations) %>% 
  mutate(count = n()) %>% 
  slice(1) %>%
  ungroup() %>% 
  arrange(migrations) %>% 
  ggplot(aes(x = yd, y = count, group = as.factor(migrations), color = as.factor(migrations))) +
  geom_segment(aes(alpha = 0.7, x=yd, xend=yd, y=0, yend=count), linewidth = 1.5) +
  geom_point(size = 3) +
  labs(x = "Day", y = "Observations", color = "Migration") +
  scale_color_manual(values = c("#0081A7", "#0098B0", "#00AFB9", "#7FC4B8", "#F7A58F", "#F27E71", "#F07167", "#ED5145", "#EE5E53")) +
  guides(alpha = "none") +
  scale_x_date(date_breaks = "1 month", date_labels =  "%m-%d") +
  theme_classic() +
  theme(text = element_text(size = 25),
        axis.line = element_line(color = "black", linewidth = 1.5),
        axis.text = element_text(color = "black"),
        legend.key.size = unit(2, 'cm'))
dev.off()

png(filename = "/home/hbronnvik/Documents/storkSSFs/figures/seasonal_uplifts.png",
    height = 18, width = 20, units = "in", res = 300)
ggplot(a_data %>% mutate(season = ifelse(season == "pre", "Spring", "Fall")), aes(as.factor(migrations), w_star, fill = season)) +
  geom_boxplot(lwd = 1.5, aes(alpha = forcats::fct_rev(as.factor(migrations)))) +
  scale_fill_manual(values = c("#FE6D5D", "#0096CC"))  +
  guides(alpha = "none") +
  labs(x = "Migrations", y = "Uplift (m/s)", fill = "Season") +
  theme_classic() +
  theme(text = element_text(size = 25),
        axis.line = element_line(color = "black", linewidth = 1.5),
        axis.text = element_text(color = "black"))
dev.off()

fac_labs <- c("Fall", "Spring")
names(fac_labs) <- c("post", "pre")
ggplot(a_data, aes(x = as.factor(migrations), y = ud_pdf, fill = as.factor(used))) +
  geom_boxplot(lwd = 1.5, aes(alpha = forcats::fct_rev(as.factor(migrations)))) +
  scale_fill_manual(values = c("#0096CC", "#FE6D5D")) + 
  guides(alpha = "none") +
  labs(x = "Migrations", y = "Conspecific density", fill = "Use") +
  theme_classic() +
  theme(text = element_text(size = 25),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black")) +
  facet_wrap(~season, labeller = labeller(season = fac_labs))

# the number of journeys that passed all the filtering
a_data %>% 
  filter(season == "pre") %>% 
  group_by(track) %>% 
  slice(1) %>% 
  ungroup() %>% 
  dplyr::select(migrations) %>% 
  table()
a_data %>% 
  filter(season == "post") %>% 
  group_by(track) %>% 
  slice(1) %>% 
  ungroup() %>% 
  dplyr::select(migrations) %>% 
  table()

# look at seasons
res_aov <- aov(w_star ~ factor(season), data = a_data)
hist(res_aov$residuals)
# not at all normal
car::qqPlot(res_aov$residuals)
wilcox.test(w_star ~ factor(season), data = a_data)
png(filename = "/home/hbronnvik/Documents/storkSSFs/figures/seasons_uplift.png",
    height = 18, width = 20, units = "in", res = 300)
print(ggplot(a_data %>% mutate(season = ifelse(season == "post", "Fall", "Spring")), aes(season, w_star, fill = season)) +
        geom_boxplot(linewidth = 2) +
        scale_color_manual(values = c("#007BA7", "#FE6F5E")) +
        labs(x = "Season", y = "Uplift potential (m/s)", fill = "Season") +
        theme_classic() +
        theme(axis.text = element_text(color = "black"),
              axis.line = element_line(linewidth = 1.5),
              text = element_text(size = 25)))
dev.off()
# look at correlation
a_data %>% 
  dplyr::select(c(ud_pdf, w_star, wind_support, step_length, turning_angle, migrations)) %>% 
  correlate() # wind_support and PDF = 0.106

# STEP 1: run the model ------------------------------------------------------------------ 
#this is based on Muff et al:
#https://conservancy.umn.edu/bitstream/handle/11299/204737/Otters_SSF.html?sequence=40&isAllowed=y#glmmtmb-1

a_d1 <- a_data1 %>% 
  mutate(stratum_ID = as.factor(stratum),
         startum_ID = as.numeric(stratum_ID),
         individual.id = as.numeric(individual.id))

TMB_struc <- glmmTMB(used ~ -1 + sqrt_ud_z*w_star_z*migrations_z + step_length_z + turning_angle_z + 
                       (1|stratum_ID) + 
                       (0 + sqrt_ud_z | individual.id) + 
                       (0 + w_star_z | individual.id) + 
                       (0 + migrations_z | individual.id), 
                     family = poisson,
                     data = a_d1, doFit = FALSE,
                     #Tell glmmTMB not to change the first standard deviation, all other values are freely estimated (and are different from each other)
                     map = list(theta = factor(c(NA,1:3))), # 3 is the n of random slopes
                     #Set the value of the standard deviation of the first random effect (here (1|startum_ID)):
                     start = list(theta = c(log(1e3),0,0,0))) #add a 0 for each random slope. in this case, 2

start_time <- Sys.time()
TMB_M <- glmmTMB:::fitTMB(TMB_struc)
Sys.time() - start_time # 3.499091 mins without interactions, 3.641424 mins with
summary(TMB_M) # fall AIC: 428226.7, spring AIC: 192566

TMB_struc <- glmmTMB(used ~ -1 + sqrt_ud_z*w_star_z*migrations_z + step_length_z + turning_angle_z + 
                       (1|stratum_ID) + 
                       (0 + sqrt_ud_z | individual.id) + 
                       (0 + w_star_z | individual.id) + 
                       (0 + migrations_z | individual.id), 
                     family = poisson, ziformula = ~1,
                     data = a_d1, doFit = FALSE,
                     #Tell glmmTMB not to change the first standard deviation, all other values are freely estimated (and are different from each other)
                     map = list(theta = factor(c(NA,1:3))), # 3 is the n of random slopes
                     #Set the value of the standard deviation of the first random effect (here (1|startum_ID)):
                     start = list(theta = c(log(1e3),0,0,0))) #add a 0 for each random slope. in this case, 2

start_time <- Sys.time()
TMB_M_zi <- glmmTMB:::fitTMB(TMB_struc)
Sys.time() - start_time 
summary(TMB_M_zi) # fall AIC: 428319.6, spring AIC: 192568.0

# the AIC is slightly better without considering zero-inflation, so we proceed with that model

lapply(split(a_data, a_data$season), function(data){
  season <- unique(data$season)
  # prep the data
  prep_d <- data %>% 
    mutate(stratum_ID = as.factor(stratum),
           stratum_ID = as.numeric(stratum_ID),
           individual.id = as.numeric(individual.id))
  # define the structure of the model
  TMB_struc <- glmmTMB(used ~ -1 + sqrt_ud_z*w_star_z*migrations_z + #step_length_z + turning_angle_z + 
                         (1|stratum_ID) + 
                         (0 + sqrt_ud_z | individual.id) + 
                         (0 + w_star_z | individual.id) + 
                         (0 + migrations_z | individual.id), 
                       family = poisson,
                       data = prep_d, doFit = FALSE,
                       #Tell glmmTMB not to change the first standard deviation, all other values are freely estimated (and are different from each other)
                       map = list(theta = factor(c(NA,1:3))), # 3 is the n of random slopes
                       #Set the value of the standard deviation of the first random effect (here (1|startum_ID)):
                       start = list(theta = c(log(1e3),0,0,0))) #add a 0 for each random slope. in this case, 2
  # fit the model
  TMB_M <- glmmTMB:::fitTMB(TMB_struc)
  # save the model
  saveRDS(TMB_M, file = paste0("/home/hbronnvik/Documents/storkSSFs/glmm_TMB_M_", season, "_", Sys.Date(), ".rds"))
  
  # the second part demands that we create values and have the model predict use of them
  
  ## make a grid of the length of the data to fill in with predicted values (one for each predictor)
  #to make sure the predictions cover the parameter space, create a dataset with all possible combinations. one per interaction term. merge later on
  grd_up <- expand.grid(x = (1:max(prep_d$migrations)),
                        y = seq(from = quantile(prep_d$w_star, 0.025, na.rm = T), to = quantile(prep_d$w_star, 0.975, na.rm = T), length.out = 15)) %>% # quantile(prep_d$w_star, .9, na.rm = T)
    rename(migrations = x,
           w_star = y) %>% 
    mutate(sqrt_ud = mean(prep_d$sqrt_ud, na.rm = T), #set other variables to their mean
           # migrations = mean(prep_d$journey_number),
           step_length = mean(prep_d$step_length, na.rm = T),
           turning_angle = mean(prep_d$turning_angle, na.rm = T),
           interaction = "uplift_migration")
  
  grd_soc <- expand.grid(x = (1:max(prep_d$migrations)),
                         y = seq(from = min(prep_d$sqrt_ud, na.rm = T), to = max(prep_d$sqrt_ud, na.rm = T), length.out = 15)) %>% # n = 135
    rename(migrations = x,
           sqrt_ud = y) %>% 
    mutate(w_star = mean(prep_d$w_star), #set other variables to their mean
           # migrations = mean(prep_d$journey_number),
           step_length = mean(prep_d$step_length, na.rm = T),
           turning_angle = mean(prep_d$turning_angle, na.rm = T),
           interaction = "ud_migration")
  
  grd_soc_up <- expand.grid(x = seq(from = min(prep_d$w_star, na.rm = T), to = max(prep_d$w_star, na.rm = T), length.out = max(prep_d$migrations)),
                            y = seq(from = min(prep_d$sqrt_ud, na.rm = T), to = max(prep_d$sqrt_ud, na.rm = T), length.out = 15)) %>% # n = 135
    rename(w_star = x,
           sqrt_ud = y) %>% 
    mutate(#set other variables to their mean
      migrations = mean(prep_d$migrations),
      step_length = mean(prep_d$step_length, na.rm = T),
      turning_angle = mean(prep_d$turning_angle, na.rm = T),
      interaction = "ud_up")
  
  grd_3 <- expand.grid(x = (1:max(prep_d$migrations)),
                       y = seq(from = quantile(prep_d$w_star, 0.025, na.rm = T), to = quantile(prep_d$w_star, 0.975, na.rm = T), length.out = 25),
                       z = seq(from = quantile(prep_d$sqrt_ud, 0.025, na.rm = T), to = quantile(prep_d$sqrt_ud, 0.975, na.rm = T), length.out = 25)) %>% # quantile(prep_d$w_star, .9, na.rm = T)
    rename(migrations = x,
           w_star = y,
           sqrt_ud = z) %>% 
    mutate(step_length = mean(prep_d$step_length, na.rm = T),
           turning_angle = mean(prep_d$turning_angle, na.rm = T),
           interaction = "uplift_migration_ud")
  
  grd_all <- bind_rows(grd_up, grd_soc, grd_soc_up, grd_3) 
  
  set.seed(770)
  n <- nrow(grd_all)
  
  new_data_only <- prep_d %>%
    group_by(stratum_ID) %>% 
    slice_sample(n = 1) %>% #randomly selects one row (from each stratum)
    ungroup() %>% 
    slice_sample(n = n, replace = F) %>% #randomly select n of these strata. Make sure n is less than n_distinct(data$stratum)
    # only keep the columns that I need
    dplyr::select(c("stratum_ID", "individual.id")) %>% 
    bind_cols(grd_all) %>% 
    #calculate z-scores. get the mean(center) and sd(scale) from the previous z transformation. for consistency
    mutate(w_star_z = (w_star - mean(a_data$w_star))/(sd(a_data$w_star)),
           sqrt_ud_z = (sqrt_ud - mean(a_data$sqrt_ud))/(sd(a_data$sqrt_ud)),
           migrations_z = (migrations - mean(a_data$migrations))/(sd(a_data$migrations)),
           step_length_z = (step_length - mean(a_data$step_length, na.rm = T))/(sd(a_data$step_length, na.rm = T)), 
           turning_angle_z = (step_length - mean(a_data$turning_angle, na.rm = T))/(sd(a_data$turning_angle, na.rm = T)))
  
  new_data <- prep_d %>% 
    mutate(interaction = "OG_data") %>% 
    dplyr::select(names(new_data_only)) %>%  #only keep the columns that are necessary for the model
    bind_rows(new_data_only) %>% 
    #accoring to the predict.glmmTMB help file: "To compute population-level predictions for a given grouping variable 
    #(i.e., setting all random effects for that grouping variable to zero), set the grouping variable values to NA."
    mutate(stratum_ID = NA)#,
           # individual.id = NA)
  
  # now that we have the values to predict, run the model on them
  preds <- predict(TMB_M, newdata = new_data, type = "link")
  
  preds_pr <- new_data %>% 
    mutate(preds = preds) %>% 
    rowwise() %>% 
    mutate(probs = gtools::inv.logit(preds)) #https://rpubs.com/crossxwill/logistic-poisson-prob
  
  inter_preds <- preds_pr %>% 
    filter(interaction != "OG_data") 
  
  saveRDS(inter_preds, file = paste0("/home/hbronnvik/Documents/storkSSFs/glmm_preds_", season, "_", Sys.Date(), ".rds"))
})

# 2023-08-23 considers individual variation, 2023-08-17 put NA for individual.id
# 2023-08-26 includes 3-way predictions at a higher resolution
# 2023-08-28 does the same as above but limits the predictions space to 2.5% to 97.5% of the values
pre_mod <- readRDS("/home/hbronnvik/Documents/storkSSFs/glmm_TMB_M_pre_2023-08-28.rds")
post_mod <- readRDS("/home/hbronnvik/Documents/storkSSFs/glmm_TMB_M_post_2023-08-28.rds")

pre_preds <- readRDS("/home/hbronnvik/Documents/storkSSFs/glmm_preds_pre_2023-08-28.rds")
post_preds <- readRDS("/home/hbronnvik/Documents/storkSSFs/glmm_preds_post_2023-08-28.rds")

# STEP 2: model validation ------------------------------------------------------------------ 
# check the outcomes

diagnose(pre_mod)
# "large negative and/or positive components in binomial or Poisson
# conditional parameters suggest (quasi-)complete separation"
# "Complete or quasi-complete separation occurs when there is a combination of 
# regressors in the model whose value can perfectly predict one or both outcomes." 
# (https://www.zeileis.org/news/biasreduction/)

# extract coefficient estimates and confidence intervals
confint(pre_mod)

# extract individual-specific random effects:
# ranef(pre_mod)[[1]]$individual.id

#calculate the RMSE
performance_rmse(pre_mod) # 0.139612

confint(post_mod)

# extract individual-specific random effects:
# ranef(pre_mod)[[1]]$individual.id

#calculate the RMSE
performance_rmse(post_mod) # 0.1381841

## Plot the interaction terms for each season

# pansy, wisteria, apricot, tangerine
colfunc <- colorRampPalette(c("#9672D5", "#C5B0E8", "#FDC4AF", "#FB8F67"))
# flame, carrot, xanthous, yellow green, apple, avocado
colfunc <- colorRampPalette(c("#6A8532", "#87A330", "#A1C349", "#F3C053", "#F9A03F", "#EB5E28"))
# cerulean, Munsell, verdigris, Tiffany, light orange, melon, salmon, bittersweet
colfunc <- colorRampPalette(c("#0081A7", "#0098b0", "#00AFB9", "#7fc4b8", "#FED9B7", "#f7a58f", "#f27e71", "#F07167", "#ee5e53"))

png(filename = "/home/hbronnvik/Documents/storkSSFs/figures/preds_glmm_spring_3_08-26.png",
    height = 18, width = 20, units = "in", res = 300)
ggplot(pre_preds %>%  
         filter(interaction == "uplift_migration_ud") %>% 
         mutate(label = paste0("Spring ", migrations))) +
  geom_raster(aes(x = w_star, y = sqrt_ud, fill = probs, group = probs), interpolate = F) +
  # geom_contour(aes(x = z, y = y, z = value), color = "black") +
  scale_x_continuous(expand = c(0, 0), n.breaks = 11) +
  scale_y_continuous(expand = c(0, 0), n.breaks = 6) +
  labs(x = "Uplift (m/s)", y = "Conspecific density") +
  scale_fill_gradientn("Selection", colours = colfunc(135)) +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 20),
        axis.line = element_line(linewidth = 1.2),
        text = element_text(size = 25),
        legend.text = element_text(size = 25),
        strip.text.x = element_text(size = 30)) +
  facet_wrap(~as.factor(label))
dev.off()

png(filename = "/home/hbronnvik/Documents/storkSSFs/figures/preds_glmm_fall_3_08-26.png",
    height = 18, width = 20, units = "in", res = 300)
ggplot(post_preds %>%  
         filter(interaction == "uplift_migration_ud") %>% 
         mutate(label = paste0("Fall ", migrations))) +
  geom_raster(aes(x = w_star, y = sqrt_ud, fill = probs, group = probs), interpolate = F) +
  # geom_contour(aes(x = z, y = y, z = value), color = "black") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Uplift (m/s)", y = "Conspecific density") +
  scale_fill_gradientn("Selection", colours = colfunc(10)) +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 20),
        axis.line = element_line(linewidth = 1.2),
        text = element_text(size = 25),
        legend.text = element_text(size = 25),
        strip.text.x = element_text(size = 30)) +
  facet_wrap(~as.factor(label))
dev.off()

inter_plots <- function(data, x_id, y_id, xlab, ylab, xmin, xmax){
  ggplot(data, aes_string(x_id, y_id, fill = "probs")) +
    geom_tile(color = "white", lwd = 1.5, linetype = 1) +
    scale_fill_gradientn(colors = colfunc(135)) +
    scale_y_continuous(expand=c(0, 0))+
    scale_x_continuous(expand=c(0, 0),
                       breaks=round(seq(from = xmin, to = xmax, length.out = 9), digits = 1))+
    labs(x = xlab, y = ylab, fill = "Selection probability") +
    theme_classic() +
    theme(axis.text = element_text(color = "black", size = 15),
          axis.line = element_line(linewidth = 1.2),
          text = element_text(size = 20),
          legend.text = element_text(size = 15),
          legend.position = "none")
}

# spring interactions
w_age_s <- inter_plots(pre_preds %>% filter(interaction == "uplift_migration"),
                       "migrations", "w_star",
                       "Number of spring migrations", "Uplift (m/s)", 
                       1, 9) + labs(title = " ")
ud_age_s <- inter_plots(pre_preds %>% filter(interaction == "ud_migration"),
                        "migrations", "sqrt_ud",
                       "Number of spring migrations", "Social density", 
                       1, 9) + labs(title = "Spring")
ud_up_s <- inter_plots(pre_preds %>% filter(interaction == "ud_up"),
                       "w_star", "sqrt_ud",
                       "Uplift (m/s)", "Social density", 
                       -1.9, 3.6) + labs(title = " ")
# fall interactions
w_age_f <- inter_plots(post_preds %>% filter(interaction == "uplift_migration"),
                       "migrations", "w_star",
                       "Number of fall migrations", "Uplift (m/s)", 
                       1, 9) + labs(title = " ")
ud_age_f <- inter_plots(post_preds %>% filter(interaction == "ud_migration"),
                        "migrations", "sqrt_ud",
                        "Number of fall migrations", "Social density", 
                        1, 9) + labs(title = "Fall")
ud_up_f <- inter_plots(post_preds %>% filter(interaction == "ud_up"),
                       "w_star", "sqrt_ud",
                       "Uplift (m/s)", "Social density", 
                       -2.3, 3.9) + labs(title = " ")
# combine the plots
blank <- ggplot(post_preds %>% filter(interaction == "ud_up"), aes(w_star, sqrt_ud, fill = probs)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) +
  scale_fill_gradientn(colors = colfunc(135)) +
  labs(x = "", y = "", fill = "Selection probability") +
  theme(text = element_text(size = 20),
        legend.text = element_text(size = 15),
        plot.margin = margin(c(0,1,0,0)))

legend <- get_legend(
  # create some space to the left of the legend
  blank + theme(legend.box.margin = margin(0, 0, 0, 10))
)

top_row <- plot_grid(ud_age_f, w_age_f, ud_up_f, labels = c("1", "2", "3"), 
                     label_x = 0.1, label_y = 0.9,
                     ncol = 4,align = 'v', axis = 'l')
bottom_row <- plot_grid(ud_age_s, w_age_s, ud_up_s, labels = c("1", "2", "3"), 
                        label_x = 0.1, label_y = 0.9, 
                        ncol = 4, align = 'v', axis = 'l')
full_plot <-  plot_grid(top_row, bottom_row, labels = c("A", "B"), ncol = 1,
                        align = 'v', axis = 'l')
full_plot + draw_grob(legend, scale = 0.7, x = 0.4)

png(filename = "/home/hbronnvik/Documents/storkSSFs/figures/preds_glmm_sqrt_ud_08-23.png",
    height = 18, width = 20, units = "in", res = 300)
print(full_plot + draw_grob(legend, scale = 0.7, x = 0.4))
dev.off()

library(plotly)
plot_ly(z = as.matrix(a_data[, c("migrations", "w_star", "sqrt_ud_z")])) %>% add_surface()

# interpolate data onto grid
test_plotly <- with(a_data[a_data$season == "post",], akima::interp(migrations, w_star, sqrt_ud,
                                        duplicate = "mean"))
# plot surface over grid
axx <- list(title = "Migrations")
axy <- list(title = "Uplift (m/s)")
axz <- list(title = "Conspecific density")
plot_ly(x = test_plotly$x, y = test_plotly$y, z = test_plotly$z,
        type = "surface") %>% 
  layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz))

p <- inter_preds %>%  
  filter(interaction == "uplift_migration_ud") %>% 
  ggplot(aes(x = "migrations", y = "w_star", z = "sqrt_ud")) +  
  geom_contour() +  
  scale_fill_distiller(palette = "Spectral", direction = -1)
ggplotly(p)

grd_3 <- expand.grid(x = (1:max(prep_d$migrations)),
                      y = seq(from = min(prep_d$w_star, na.rm = T), to = max(prep_d$w_star, na.rm = T), length.out = 15),
                     z = seq(from = min(prep_d$sqrt_ud, na.rm = T), to = max(prep_d$sqrt_ud, na.rm = T), length.out = 15)) %>% # quantile(prep_d$w_star, .9, na.rm = T)
  rename(migrations = x,
         w_star = y,
         sqrt_ud = z) %>% 
  mutate(step_length = mean(prep_d$step_length, na.rm = T),
         turning_angle = mean(prep_d$turning_angle, na.rm = T),
         interaction = "uplift_migration_ud")

plot3D::scatter3D(x = i_data$migrations, y = i_data$w_star, z = i_data$sqrt_ud, colvar = i_data$probs)
rgl::surface3d(x = i_data$migrations, y = i_data$w_star, z = i_data$sqrt_ud)

i_data <- inter_preds %>%  
  filter(interaction == "uplift_migration_ud") %>% 
  dplyr::select("migrations", "sqrt_ud", "w_star", "probs") %>% 
  rename(x = "migrations", 
         y = "sqrt_ud",
         z = "w_star",
         value = "probs")

p <- ggplot(i_data, aes(z, y, fill = value)) +
  geom_tile(color = "white", lwd = 0, linetype = 1) +
  scale_fill_gradientn(colors = colfunc(135)) +
  # scale_y_continuous(expand=c(0, 0))+
  # scale_x_continuous(expand=c(0, 0),
                     # breaks=round(seq(from = min(y), to = max(y), length.out = 9), digits = 1))+
  labs(x = "Uplift (m/s)", y = "Conspecific density", fill = "Selection probability") +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 15),
        axis.line = element_line(linewidth = 1.2),
        text = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.position = "none") +
  facet_wrap(~x, strip.position="left") #+
  # gganimate::transition_states(x,
  #                   transition_length = 4,
  #                   state_length = 0.2,
  #                   wrap = F) +
  # gganimate::exit_fade() +
  # labs(title = "Migration {closest_state}") 

png(filename = "/home/hbronnvik/Documents/storkSSFs/figures/facet_spring_3d.png",
    height = 18, width = 20, units = "in", res = 500)
plot_gg(p,multicore = TRUE,width=10,height=8,scale=250,windowsize=c(1400,866),
        zoom = 0.8, phi = 60, theta = -30)
render_snapshot(clear = TRUE)
dev.off()

library(raster)
xyz <- rasterFromXYZ(i_data[1:8, 1:3], res=c(NA,NA), crs="", digits=5)
xyz <- raster(xmn = min(i_data$x), xmx = max(i_data$x),
       ymn = min(i_data$y), ymx = max(i_data$y),
       ncols = 8, nrows = 15)
raster::setValues(xyz, i_data$z)
rasterVis::plot3D(xyz)
# install.packages("rayshader")
library(rayshader)
library(ggplot2)
heatmap_fill <- ggplot(i_data %>% group_by(x, y) %>% slice_sample(n = 1) %>% ungroup()) +
  geom_tile(aes(x, y, fill = value, color = z)) +
  scale_fill_gradientn(colors = colfunc(135)) +
  theme_classic() +
  theme(legend.position = "none")
heatmap_height <- ggplot(i_data %>% group_by(x, y) %>% slice_sample(n = 1) %>% ungroup()) +
  geom_tile(aes(x, y, fill = z, color = z)) +
  scale_fill_gradientn(colors = colfunc(135)) +
  theme_classic() +
  theme(legend.position = "none")
# rayshader::plot_gg(plot, width = 7, height = 4, raytrace = FALSE, preview = TRUE)
rayshader::plot_gg(heatmap_fill, heatmap_height, multicore = TRUE, raytrace = TRUE, width = 7, height = 4, 
        scale = 300, windowsize = c(1400, 866), zoom = 0.6, phi = 30, theta = 30);
render_snapshot(clear = TRUE)



ex <- expand.grid(x = 1:8,
                  y = seq(from = -1.85, to = 3.58, length.out = 15),
                  z = seq(from = 0, to = 2.41, length.out = 15)) %>% 
  mutate(value = sample(seq(from = 0, to = 1, length.out = 3000), 1800))

vol <- simplify2array(by(i_data, i_data$x, as.matrix))

library(misc3d)
con <- misc3d::computeContour3d(vol, max(vol), 1)
misc3d::drawScene(makeTriangles(con))
contour3d(vol, 1)



pre_graph <- confint(pre_mod) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Factor") %>% 
  filter(!grepl("id", Factor)) %>% 
  mutate(Variable = c("Conspecific density", "Uplift", "Spring migrations", #"Step length", "Turning angle", 
                      "Conspecifics X Uplift", "Conspecifics X Migrations", "Uplift X Migrations", "Conspecifics X Uplift X Migrations"))
colnames(pre_graph)[2:3] <- c("Lower", "Upper") 

pre_graph$p <- as.data.frame(summary(pre_mod)[[6]][1])[,4]
pre_graph$significance <- ifelse(pre_graph$p < 0.001, "***", 
                                  ifelse(pre_graph$p > 0.001 & pre_graph$p < 0.01, "**",
                                         ifelse(pre_graph$p > 0.01 & pre_graph$p < 0.05, "*",
                                                "")))
pre_graph <- pre_graph %>% 
  mutate(Variable = factor(Variable, levels = c("Conspecifics X Uplift X Migrations", "Uplift X Migrations",
                                                "Conspecifics X Uplift", "Conspecifics X Migrations", "Spring migrations", 
                                                "Uplift", "Conspecific density"))) %>% 
  arrange(Variable)

post_graph <- confint(post_mod) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Factor") %>% 
  filter(!grepl("id", Factor)) %>% 
  mutate(Variable = c("Conspecific density", "Uplift", "Fall migrations", #"Step length", "Turning angle", 
                      "Conspecifics X Uplift", "Conspecifics X Migrations", "Uplift X Migrations", "Conspecifics X Uplift X Migrations"))
colnames(post_graph)[2:3] <- c("Lower", "Upper") 

post_graph$p <- as.data.frame(summary(post_mod)[[6]][1])[,4]
post_graph$significance <- ifelse(post_graph$p < 0.001, "***", 
                                  ifelse(post_graph$p > 0.001 & post_graph$p < 0.01, "**",
                                         ifelse(post_graph$p > 0.01 & post_graph$p < 0.05, "*",
                                                "")))
post_graph <- post_graph %>% 
  mutate(Variable = factor(Variable, levels = c("Conspecifics X Uplift X Migrations", "Uplift X Migrations",
                                          "Conspecifics X Uplift", "Conspecifics X Migrations", "Fall migrations", 
                                           "Uplift", "Conspecific density"))) %>% 
  arrange(Variable)

post_coefs <- ggplot(post_graph %>% filter(Variable != "Fall migrations"), aes(Estimate, Variable)) +
  geom_vline(xintercept = 0, lty = 2, linewidth = 1) +
  geom_pointrange(aes(xmin = Lower, xmax = Upper), #fill = significance_05, color = significance_05), 
                  linewidth = 1, size = 1) +
  scale_color_manual(values = c("#bd68ee", "#f48c3d")) +
  labs(y = "", fill = "Significant", color = "Significant", title = "Fall") +
  theme_classic() +
  annotate(geom="text", x=16, y=c(1:6), label=post_graph$significance[post_graph$Variable != "Fall migrations"], color="black", size = 6) +
  theme(axis.text = element_text(color = "black", size = 15),
        axis.line = element_line(linewidth = 1),
        text = element_text(size = 15),
        legend.text = element_text(size = 10))
pre_coefs <- ggplot(pre_graph %>% filter(Variable != "Spring migrations"), aes(Estimate, Variable)) +
  geom_vline(xintercept = 0, lty = 2, linewidth = 1) +
  geom_pointrange(aes(xmin = Lower, xmax = Upper), #fill = significance_05, color = significance_05), 
                  linewidth = 1, size = 1) +
  scale_color_manual(values = c("#bd68ee", "#f48c3d")) +
  labs(y = "", fill = "Significant", color = "Significant", title = "Spring") +
  theme_classic() +
  annotate(geom="text", x=5, y=c(1:6), label=pre_graph$significance[pre_graph$Variable != "Spring migrations"], color="black", size = 6) +
  theme(axis.text = element_text(color = "black", size = 15),
        axis.line = element_line(linewidth = 1),
        text = element_text(size = 15),
        legend.text = element_text(size = 10))

coefs_plot <-  plot_grid(post_coefs, pre_coefs, labels = c("A", "B"), nrow = 2,
                         align = 'v', axis = 'l')
coefs_plot

png(filename = "/home/hbronnvik/Documents/storkSSFs/figures/coefs_glmm.png",
    height = 150, width = 300, units = "mm", res = 500)
coefs_plot
dev.off()


# look at homogeneity
fac_labs <- c("Fall", "Spring")
names(fac_labs) <- c("post", "pre")
w_star_var <- ggplot(a_data, aes(as.factor(migrations), w_star, group = as.factor(migrations), 
                                 fill = season, alpha = forcats::fct_rev(as.factor(migrations))))+
  geom_boxplot() +
  scale_fill_manual(values = c("#FE6F5E", "#007BA7")) +
  scale_alpha_discrete(
    range = c(0.3, 0.9),
    guide = guide_legend(override.aes = list(fill = "black"))) +
  labs(x = "", y = "Uplift potential (m/s)", fill = "Season") +
  theme_classic() +
  theme(text = element_text(size = 20, color = "black"),
        axis.text = element_text(size = 15, color = "black"),
        legend.position = "none") +
  facet_wrap(~season, labeller = labeller(season = fac_labs))
dens_var <- ggplot(a_data, aes(as.factor(migrations), ud_pdf, group = as.factor(migrations), 
                     fill = season, alpha = forcats::fct_rev(as.factor(migrations)))) +
  geom_boxplot() +
  scale_fill_manual(values = c("#FE6F5E", "#007BA7")) +
  scale_alpha_discrete(
    range = c(0.3, 0.9),
    guide = guide_legend(override.aes = list(fill = "black"))) +
  labs(x = "Migration", y = "Conspecific density", fill = "Season") +
  theme_classic() +
  theme(text = element_text(size = 20, color = "black"),
        axis.text = element_text(size = 15, color = "black"),
        legend.position = "none") +
  facet_wrap(~season, labeller = labeller(season = fac_labs))

vars_plot <-  plot_grid(w_star_var, dens_var, labels = c("A", "B"), ncol = 1,
                         align = 'v', axis = 'l')

png(filename = "/home/hbronnvik/Documents/storkSSFs/figures/distributions.png",
    height = 200, width = 300, units = "mm", res = 500)
vars_plot
dev.off()

# check the stratum variance
v <- a_data %>% 
  mutate(date = yday(alignment)) %>% 
  group_by(stratum) %>% 
  mutate(varU = var(na.omit(w_star)),
         varS = var(na.omit(ud_pdf)),
         logVarU = log(varU),
         logVarS = log(varS)) %>% 
  ungroup() %>% 
  group_by(date) %>% 
  mutate(n = length(unique(stratum))) %>% 
  ungroup() %>% 
  mutate(semU = varU/sqrt(n),
         semS = varS/sqrt(n))
# png("C:/Users/hbronnvik/Documents/storkSSFs/stratum_sdsqrtn_blhz_per_date.png", 
# width = 33.87, height = 19.05, units = "cm", res = 500)
print(ggplot(v, aes(migrations, varU, group = migrations)) +
        geom_violin() +
        # geom_bar(stat = "identity") +
        labs(x = "Migration", y = "Stratum variance") +
        theme_classic() +
        theme(text = element_text(size = 20, color = "black"),
              axis.text = element_text(size = 15, color = "black")))
# dev.off()


ggplot(v, aes(as.factor(migrations), logVarU, group = as.factor(migrations), 
                   fill = season, alpha = forcats::fct_rev(as.factor(migrations))))+
  geom_boxplot() +
  scale_fill_manual(values = c("#FE6F5E", "#007BA7")) +
  scale_alpha_discrete(
    range = c(0.3, 0.9),
    guide = guide_legend(override.aes = list(fill = "black"))) +
  labs(x = "", y = "Log variance in available uplift", fill = "Season") +
  theme_classic() +
  theme(text = element_text(size = 20, color = "black"),
        axis.text = element_text(size = 15, color = "black"),
        legend.position = "none") +
  facet_wrap(~season, labeller = labeller(season = fac_labs))

library(ggridges)

ggplot(v, aes(y = forcats::fct_rev(as.factor(migrations)), x = logVarU, group = as.factor(migrations),
              fill = forcats::fct_rev(as.factor(migrations)), alpha = forcats::fct_rev(as.factor(migrations)))) +
  stat_density_ridges(quantile_lines = TRUE, alpha = 0.75) +
  scale_fill_manual(values = c("#0081A7", "#0092AD", "#00A3B4", "#1FB4B8", "#7FC4B8", "#DED3B7", "#FABEA3", "#F59884", "#F07167")) +
  labs(x = "Log variance in available uplift", y = "Migration", fill = "Migration") +
  theme_classic()

lapply(split(v[which(v$season == "post"),], v$migrations[which(v$season == "post")]), function(x){
  print(sd(x$varU))
})

ggplot(v %>% group_by(migrations, season) %>% mutate(m = mean(varU)) %>% slice(1) %>% ungroup(), 
       aes(migrations, m)) +
  geom_point() + 
  geom_smooth(method = "lm") +
  labs(x = "Migrations", y = "Mean per-stratum variance") +
  theme_classic() + 
  facet_wrap(~season, labeller = labeller(season = fac_labs))
vv <- v %>% 
  group_by(migrations, season) %>% 
  mutate(mU = mean(varU),
         sdU = sd(varU),
         mS = mean(varS),
         sdS = sd(varS)) %>% 
  slice(1) %>% 
  ungroup()

r <- c(format(summary(lm(mU ~ migrations, vv[vv$season == "post",]))$r.squared, digits = 3),
       format(summary(lm(mU ~ migrations, vv[vv$season == "pre",]))$r.squared, digits = 3),
       format(summary(lm(mS ~ migrations, vv[vv$season == "post",]))$r.squared, digits = 3),
       format(summary(lm(mS ~ migrations, vv[vv$season == "pre",]))$r.squared, digits = 3))
graph_labels <- data.frame(season = c("post", "pre", "post", "pre"),
                          r_val = paste("R =", c(r[1], r[2], r[3], r[4])))


cor.test(x=v$migrations[v$season == "post"], y=v$varU[v$season == "post"], method = 'spearman')
cor.test(x=v$migrations[v$season == "pre"], y=v$varU[v$season == "pre"], method = 'spearman')
cor.test(x=v$migrations[v$season == "post"], y=v$varS[v$season == "post"], method = 'spearman')
cor.test(x=v$migrations[v$season == "pre"], y=v$varS[v$season == "pre"], method = 'spearman')
graph_labels <- data.frame(season = c("post", "pre", "post", "pre"),
                           r_val = paste("rho =", c(0.036, -0.025, -0.035, -0.305)),
                           p_val = paste("p < ", c("2.2e-16", "2.2e-16", "2.2e-16", "2.2e-16")))
w <- ggplot(vv, aes(migrations, mU)) +
  geom_errorbar(aes(ymin = mU-sdU, ymax = mU+sdU)) +
  geom_point() + 
  geom_smooth(method = "lm", aes(color = season)) +
  scale_color_manual(values = c("#FE6F5E", "#007BA7")) +
  labs(x = "Migrations", y = "Per-stratum variance in w*") +
  theme_classic() + 
  theme(legend.position = "none") +
  facet_wrap(~season, labeller = labeller(season = fac_labs)) + 
  geom_text(data = graph_labels[1:2,], aes(x = 8.5, y = 0.3, label = p_val)) + 
  geom_text(data = graph_labels[1:2,], aes(x = 8.5, y = 0.35, label = r_val))
s <- ggplot(vv, aes(migrations, mS)) +
  geom_errorbar(aes(ymin = mS-sdS, ymax = mS+sdS)) +
  geom_point() + 
  geom_smooth(method = "lm", aes(color = season)) +
  scale_color_manual(values = c("#FE6F5E", "#007BA7")) +
  labs(x = "Migrations", y = "Per-stratum variance in social density") +
  theme_classic() + 
  theme(legend.position = "none") +
  facet_wrap(~season, labeller = labeller(season = fac_labs)) + 
  geom_text(data = graph_labels[3:4,], aes(x = 8.5, y = max(vv$sdS), label = p_val)) + 
  geom_text(data = graph_labels[3:4,], aes(x = 8.5, y = max(vv$sdS)+.4e-25, label = r_val))
png("/home/hbronnvik/Documents/storkSSFs/figures/stratum_vars_migration.png",
    width = 22.58, height = 12.7, units = "cm", res = 500)
plot_grid(w, s, labels = c("A", "B"), ncol = 1,
          align = 'v', axis = 'l', label_x = 0.1, label_y = 1)
dev.off()
# rank test
v <- v %>% 
  mutate(varU_scaled = scale(varU),
         varS_scaled = scale(varS))
cor.test(x=v$migrations, y=v$varS, method = 'spearman')
ggplot(v, aes(x=migrations, y=varU)) + 
  geom_point(color='#2980B9', size = 4) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, color='#2C3E50') +
  theme_classic()
ggplot(v, aes(migrations, varU, group = migrations)) +
  geom_boxplot() +
  theme_classic()

lme4::lmer(w_star~migrations+(1|stratum), data = v)

nlme::gls(varU ~ migrations, v)
summary(nlme::gls(varU ~ as.factor(migrations), v[v$season == "post",]))

#Welch’s Test for Unequal Variances
t.test(v$varU[v$journey_number == 1], v$varU[v$journey_number == 5])
#Bartlett's test 
bartlett.test(varU ~ journey_number, data = v)
#perform Welch's ANOVA
model <- oneway.test(blh ~ factor(journey_number), data = mod_data, var.equal = FALSE)
#perform a normal ANOVA
model <- aov(varU ~ factor(migrations), data = v[v$season == "post",])
#load DescTools package
library(DescTools)
#perform Scheffe's test
ScheffeTest(model)
# compute Kruskal-Wallis test
kruskal.test(varU ~ factor(migrations), data = v)
pairwise.wilcox.test(v$varU, factor(v$migrations),
                     p.adjust.method = "BH")


























## make a grid of the length of the data to fill in with predicted values (one for each predictor)
#to make sure the predictions cover the parameter space, create a dataset with all possible combinations. one per interaction term. merge later on
grd_up <- expand.grid(x = (1:9),
                       y = seq(from = min(a_data1$w_star, na.rm = T), to = quantile(a_data1$w_star, .9, na.rm = T), length.out = 15)) %>% # n = 135
  rename(migrations = x,
         w_star = y) %>% 
  mutate(sqrt_ud = mean(a_data1$sqrt_ud, na.rm = T), #set other variables to their mean
         # migrations = mean(a_data1$journey_number),
         step_length = mean(a_data1$step_length, na.rm = T),
         turning_angle = mean(a_data1$turning_angle, na.rm = T),
         interaction = "uplift_migration")

grd_soc <- expand.grid(x = (1:9),
                      y = seq(from = min(a_data1$sqrt_ud, na.rm = T), to = quantile(a_data1$sqrt_ud, .9, na.rm = T), length.out = 15)) %>% # n = 135
  rename(migrations = x,
         sqrt_ud = y) %>% 
  mutate(w_star = mean(a_data1$w_star), #set other variables to their mean
         # migrations = mean(a_data1$journey_number),
         step_length = mean(a_data1$step_length, na.rm = T),
         turning_angle = mean(a_data1$turning_angle, na.rm = T),
         interaction = "ud_migration")

grd_soc_up <- expand.grid(x = seq(from = min(a_data1$w_star, na.rm = T), to = quantile(a_data1$w_star, .9, na.rm = T), length.out = 9),
                       y = seq(from = min(a_data1$sqrt_ud, na.rm = T), to = quantile(a_data1$sqrt_ud, .9, na.rm = T), length.out = 15)) %>% # n = 135
  rename(w_star = x,
         sqrt_ud = y) %>% 
  mutate(#set other variables to their mean
         migrations = mean(a_data1$migrations),
         step_length = mean(a_data1$step_length, na.rm = T),
         turning_angle = mean(a_data1$turning_angle, na.rm = T),
         interaction = "ud_up")

grd_all <- bind_rows(grd_up, grd_soc, grd_soc_up) 

set.seed(770)
n <- nrow(grd_all)

new_data_only <- a_data1 %>%
  group_by(stratum) %>% 
  slice_sample(n = 1) %>% #randomly selects one row (from each stratum)
  ungroup() %>% 
  slice_sample(n = n, replace = F) %>% #randomly select n of these strata. Make sure n is less than n_distinct(data$stratum)
  #only keep the columns that I need
  dplyr::select(c("stratum", "individual.id")) %>% 
  bind_cols(grd_all) %>% 
  #calculate z-scores. get the mean(center) and sd(scale) from the previous z transformation. for consistency
  mutate( w_star_z = (w_star - mean(a_data$w_star))/(sd(a_data$w_star)),
          sqrt_ud_z = (sqrt_ud - mean(a_data$sqrt_ud))/(sd(a_data$sqrt_ud)),
          migrations_z = (migrations - mean(a_data$migrations))/(sd(a_data$migrations)),
          step_length_z = (step_length - mean(a_data$step_length, na.rm = T))/(sd(a_data$step_length, na.rm = T)), 
          turning_angle_z = (step_length - mean(a_data$turning_angle, na.rm = T))/(sd(a_data$turning_angle, na.rm = T)))

new_data <- a_data1 %>% 
  mutate(interaction = "OG_data") %>% 
  dplyr::select(names(new_data_only)) %>%  #only keep the columns that are necessary for the model
  bind_rows(new_data_only) %>% 
  #accoring to the predict.glmmTMB help file: "To compute population-level predictions for a given grouping variable 
  #(i.e., setting all random effects for that grouping variable to zero), set the grouping variable values to NA."
  mutate(stratum = NA,
         individual.id = NA)

# now that we have the values to predict, run the model on them
preds <- predict(TMB_M, newdata = new_data, type = "link")

preds_pr <- new_data %>% 
  mutate(preds = preds) %>% 
  rowwise() %>% 
  mutate(probs = gtools::inv.logit(preds)) #https://rpubs.com/crossxwill/logistic-poisson-prob

inter_preds <- preds_pr %>% 
  filter(interaction != "OG_data") %>% 
  mutate(migrations = as.factor(migrations),
         # w_star = as.factor(round(w_star, digits = 2)),
         ud_pdf = as.factor(round(ud_pdf*1000000000000, digits = 2)))

colfunc <- colorRampPalette(c("#9672D5", "#FB8F67"))
colfunc(13)


build <- inter_preds %>% 
  select("w_star", "migrations", "ud_pdf", "probs") %>% 
  mutate(ud_pdf = as.numeric(ud_pdf),
         w_star = as.numeric(w_star),
         migrations = as.numeric(migrations)) %>% 
  pivot_longer(!probs,
               names_to = "variable", 
               values_to = "probs")

png(filename = "/home/hbronnvik/Documents/storkSSFs/figures/interaction_up_falls_glmm.png",
    width = 210, height = 297, units = "mm", res = 200)
w_age <- ggplot(inter_preds %>% filter(interaction == "uplift_migration"), aes(migrations, w_star, fill = probs)) +
        geom_tile(color = "white",
                  lwd = 1.5,
                  linetype = 1) +
        # scale_fill_gradient(low = "cornflowerblue",
        #                     # mid = "white",
        #                     high = "firebrick") +
        
        scale_fill_gradientn(colors = colfunc(135)) +
        # coord_fixed() +
        scale_y_continuous(expand=c(0, 0))+
        scale_x_continuous(expand=c(0, 0),
                           breaks=seq(from = 1, to = 9, length.out = 9))+
        labs(x = "Number of fall migrations", y = "Uplift (m/s)", fill = "Selection probability") +
        theme_classic() +
        theme(axis.text = element_text(color = "black", size = 15),
              axis.line = element_line(linewidth = 1.2),
              text = element_text(size = 20),
              legend.text = element_text(size = 15),
              legend.position = "none")
dev.off()

png(filename = "/home/hbronnvik/Documents/storkSSFs/figures/interaction_ud_falls_glmm.png",
    width = 210, height = 297, units = "mm", res = 200)
ud_age <- ggplot(inter_preds %>% filter(interaction == "ud_migration"), aes(migrations, sqrt_ud, fill = probs)) +
        geom_tile(color = "white",
                  lwd = 1.5,
                  linetype = 1) +
        # scale_fill_gradient(low = "cornflowerblue",
        #                     # mid = "white",
        #                     high = "firebrick") +
        
        scale_fill_gradientn(colors = colfunc(135)) +
        # coord_fixed(ratio = 1) +
        scale_y_continuous(expand=c(0, 0))+
        scale_x_continuous(expand=c(0, 0),
                           breaks=seq(from = 1, to = 9, length.out = 9))+
        labs(x = "Number of fall migrations", y = "Social density", fill = "Selection probability") +
        theme_classic() +
        theme(axis.text = element_text(color = "black", size = 15),
              axis.line = element_line(linewidth = 1.2),
              text = element_text(size = 20),
              legend.text = element_text(size = 15),
              legend.position = "none")
dev.off()

png(filename = "/home/hbronnvik/Documents/storkSSFs/figures/interaction_up_ud_glmm.png",
    width = 210, height = 297, units = "mm", res = 200)
ud_up <- ggplot(inter_preds %>% filter(interaction == "ud_up"), aes(w_star, sqrt_ud, fill = probs)) +
        geom_tile(color = "white",
                  lwd = 1.5,
                  linetype = 1) +
        # scale_fill_gradient(low = "cornflowerblue",
        #                     # mid = "white",
        #                     high = "firebrick") +
        scale_fill_gradientn(colors = colfunc(135)) +
        # coord_fixed() +
        scale_y_continuous(expand=c(0, 0))+
        scale_x_continuous(expand=c(0, 0),
                         breaks=seq(from= -2, to = 4, length.out = 9))+
        labs(x = "Uplift (m/s)", y = "Social density", fill = "Selection probability") +
        theme_classic() +
        theme(axis.text = element_text(color = "black", size = 15),
              axis.line = element_line(linewidth = 1.2),
              text = element_text(size = 20),
              legend.text = element_text(size = 15),
              legend.position = "none")
dev.off()

png(filename = "/home/hbronnvik/Documents/storkSSFs/figures/fall_glmm.png",
    height = 150, width = 300, units = "mm", res = 200)
coef_p <- ggplot(graph, aes(Variable, Estimate)) +
  geom_hline(yintercept = 0, lty = 2, linewidth = 1.5) +
  geom_pointrange(mapping = aes(ymin = Lower, ymax = Upper), 
                  linewidth = 2, size = 1.5) +
  # scale_color_manual(values = c("#9672D5", "#FB8F67")) +
  labs(y = "Estimate", x = "", fill = "Significant", color = "Significant") +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 15),
        axis.line = element_line(linewidth = 1.2),
        text = element_text(size = 20),
        legend.text = element_text(size = 15)) + 
  annotate(geom="text", x=c(1,2,4,6), y=25, label="***", color="black", size = 10) + 
  # annotate(geom="text", x=c(2), y=25, label="**", color="black", size = 10) + 
  # annotate(geom="text", x=c(3), y=25, label="*", color="black", size = 10) + + 
  annotate(geom="text", x=c(3), y=25, label=".", color="black", size = 10)+ 
  annotate(geom="text", x=c(1,2,3,4,5,6,7), y=22, label=c("11.2", "0.417", "-0.046",
                                                          "-0.099", "0.891", "0.615",
                                                          "0.026"), 
           color="black", size = 8)
dev.off()

## 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

blank <- ggplot(inter_preds %>% filter(interaction == "ud_up"), aes(w_star, sqrt_ud, fill = probs)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) +
  scale_fill_gradientn(colors = colfunc(135)) +
  coord_cartesian(xlim = c(0.4,0.5)) +
  geom_rect(aes(xmin = min(w_star)-1, xmax = max(w_star)+1, ymin = min(sqrt_ud)-1, ymax = max(sqrt_ud)+1), fill = "white") +
  labs(x = "", y = "", fill = "Selection probability") +
  theme_void() +
  theme(text = element_text(size = 20),
        legend.text = element_text(size = 15),
        plot.margin = margin(c(0,1,0,0)))

library(cowplot)
legend <- get_legend(
  # create some space to the left of the legend
  blank + theme(legend.box.margin = margin(0, 0, 0, 10))
)
bottom_row <- plot_grid(ud_age, w_age, ud_up, legend, labels = c("1", "2", "3"), ncol = 4, 
          align = 'v', axis = 'l')
full_plot <-  plot_grid(coef_p, bottom_row, labels = "AUTO", ncol = 1,
          align = 'v', axis = 'l')

png(filename = "/home/hbronnvik/Documents/storkSSFs/figures/combi_glmm_sqrt_ud.png",
    height = 18, width = 20, units = "in", res = 300)
print(full_plot)
dev.off()

# individual coefficients
ind_coefs <- coef(TMB_M_zi)$cond$individual.id
ind_coefs <- ind_coefs %>% 
  rownames_to_column()
ggplot(ind_coefs, aes(w_star_z, rowname)) +
  geom_point() +
  labs(x = "w*", y = "Individual") +
  theme_classic()
ggplot(ind_coefs, aes(sqrt_ud_z, rowname)) +
  geom_point() +
  labs(x = "PDF", y = "Individual") +
  theme_classic()
