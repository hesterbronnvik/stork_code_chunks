### Estimate population distributions and model
### Hester Bronnvik
### 2023-07-13
### hbronnvik@ab.mpg.de

library(move)
library(ctmm)
library(lubridate)
library(raster)
library(terra)
library(tidyverse)
library(data.table)

# functions: 
wind_support <- function(u,v,heading) {
  angle <- atan2(u,v) - heading/180*pi
  return(cos(angle) * sqrt(u*u+v*v))
  14
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
  
  w_star <- CubeRoot(g*z*(wthetav/Thetav_k_z))
  
  return(w_star)
  
}

## 1) Get the data

# read in the saved data for the 509 birds for all locations during migration
locs <- readRDS("/home/hbronnvik/Documents/storkSSFs/soc_info_migration_locations_2023_06_12.rds")

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
a_data$w_star <- convective(a_data)

summary(a_data$w_star)

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
colMeans(a_data[,c(15:28)])
level <- colnames(a_data[,c(15:28)][which(abs(colMeans(a_data[,c(15:28)]) - mean(heights$ha_msl)) == min(abs(colMeans(a_data[,c(15:28)]) - mean(heights$ha_msl))))])
# this is the level to use:
level
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
wind_files <- list.files("/home/hbronnvik/Documents/storkSSFs/ecmwf", pattern = "wind_filled", full.names = T)



