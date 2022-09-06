### Finding multi-year birds for SSF analysis
### Hester Bronnvik
### 10.08.2022

library(move)
library(amt)
library(lubridate)
library(sf)
library(raster)
library(rgeos)
library(rasterVis)
library(countrycode)
library(tidyverse)

# set wd and load login data
setwd("C:/Users/heste/Desktop/HB_storks")
load("loginStored.rdata")

# collect the names of the stork studies to use
studies <- c(24442409, 212096177, 76367850, 21231406, 1176017658, 173641633)

#import required functions
NCEP.loxodrome.na <- function (lat1, lat2, lon1, lon2) {
  deg2rad <- pi/180
  acot <- function(x) {
    return(atan(1/x))
  }
  lat1 <- deg2rad * lat1
  lat2 <- deg2rad * lat2
  lon1 <- deg2rad * lon1
  lon2 <- deg2rad * lon2
  deltaLon <- lon2 - lon1
  pi4 <- pi/4
  Sig1 <- log(tan(pi4 + lat1/2))
  Sig2 <- log(tan(pi4 + lat2/2))
  deltaSig <- Sig2 - Sig1
  if (deltaLon == 0 && deltaSig > 0) {
    head <- 0
  }
  else if (deltaLon == 0 && deltaSig < 0) {
    head <- 180
  }
  else if (deltaSig == 0 && deltaLon > 0) {
    head <- 90
  }
  else if (deltaSig == 0 && deltaLon < 0) {
    head <- 270
  }
  else if (deltaSig < 0 && deltaLon < 0) {
    head <- acot(deltaSig/deltaLon) * 180/pi + 180
  }
  else if (deltaSig < 0 && deltaLon > 0) {
    head <- acot(deltaSig/deltaLon) * 180/pi + 180
  }
  else if (deltaSig > 0 && deltaLon > 0) {
    head <- acot(deltaSig/deltaLon) * 180/pi
  }
  else if (deltaSig > 0 && deltaLon < 0) {
    head <- acot(deltaSig/deltaLon) * 180/pi + 360
  }
  else {
    head <-NA}
  return(head)
}
  
# lapply(studies, function(x){
#   cols <- getMovebankSensorsAttributes(x, login = loginStored)
#   print(x)
#   print("height_above_msl"  %in% cols$short_name)
# })

## download DEM data for the countries along the flyway 
countries <- c("Niger", "Denmark", "Western Sahara", "Germany", "Spain", "Morocco", "Portugal", "Belgium", "Netherlands",
               "Algeria", "Liechtenstein", "Luxembourg", "Mauritania", "Mali", "Senegal", 
               "Gambia", "Ivory Coast", "Ghana", "Burkina Faso", "Switzerland", "France", "Benin") %>% 
  countrycode(origin = 'country.name', destination = 'iso3c')

shp <- shapefile("C:/Users/heste/Desktop/HB_storks/srtm/tiles.shp")
plot(shp)

#Get country geometry first
outlines <- lapply(countries, function(x){
  print(x)
  temp <- getData("GADM", country = x, level=0)
  crs(temp) <- crs(shp)
  
  return(temp)
}) %>% reduce(rbind)


#Intersect country geometry and tile grid
# intersects <- gIntersects(outlines, shp, byid=T, prepared=F)
# tiles <- shp[intersects[,1],]
tiles <- raster::intersect(outlines, shp)

srtm_list  <- list()
for(i in 1:length(tiles)) {
  lon <- extent(tiles[i,])[1]  + (extent(tiles[i,])[2] - extent(tiles[i,])[1]) / 2
  lat <- extent(tiles[i,])[3]  + (extent(tiles[i,])[4] - extent(tiles[i,])[3]) / 2
  
  tile <- getData('SRTM', lon=lon, lat=lat)
  
  srtm_list[[i]] <- tile
}

#Mosaic tiles
srtm_list$fun <- mean 
#srtm_mosaic   <- do.call(mosaic, srtm_list)
srtm_mosaic <- raster::raster("srtm_mosaic")


#Crop tiles to country borders
srtm_crop     <- mask(srtm_mosaic, outlines)

#Plot
p <- levelplot(srtm_mosaic)
p + latticeExtra::layer(sp.lines(outlines, lwd=0.8, col='darkgray'))


## determine the identities of the nestling birds
nestlings <- lapply(studies, function(x){
  
  print(x)
  
  birds <- getMovebankReferenceTable(study = x, login = loginStored) %>% 
    filter(sensor_type_id == 653) %>% 
    filter(!is.na(animal_id))
  
  if(length(unique(birds$animal_id[grep("juv|chick|young", birds$animal_comments, fixed = F)])) > 0){
    nesties1 <- unique(birds$animal_id[grep("juv|chick|young", birds$animal_comments, fixed = F)])
  } else {nesties1 <- NA}
  
  if("animal_life_stage" %in% colnames(birds) & length(unique(birds$animal_id[which(birds$animal_life_stage == "nestling")])) > 0){
    nesties2 <- unique(birds$animal_id[which(birds$animal_life_stage == "nestling")])
  } else {nesties2 <-  NA}
  
  nesties <- na.omit(unique(c(unique(nesties1), unique(nesties2))))
  
  return(nesties)
  
}) %>% unlist()

## determine approximate nest locations for each individual (adapted from Flack dBBMMs)
nest_sites <- lapply(studies, function(x){
  
  # get the names of the animals to download
  days <- getMovebankAnimals(study =  x, login = loginStored) %>% 
    # remove errors and duplicated information
    filter(timestamp_start != "" & sensor_type_id == 653 & individual_id %in% nestlings) %>% 
    # reduce the variables
    dplyr::select(individual_id, timestamp_start) %>% 
    # calculate the time each bird transmitted for
    mutate(three_days = as.POSIXct(timestamp_start, tz = "UTC") + (60*60*24*3))   
  
  nests <- lapply(unique(days$individual_id), function(y){
    print(which(unique(days$individual_id) == y))
    three_days <- getMovebankData(study = x, animalName =  y, sensorID = "GPS", login = loginStored,
                                     timestamp_start = as.POSIXct(days$timestamp_start[which(days$individual_id == y,)], tz = "UTC"),
                                     timestamp_end = days$three_days[which(days$individual_id == y)],
                                  removeDuplicatedTimestamps =T)
    
    nest <- geomean(as.matrix(coordinates(three_days)))
    
    nest_loc <- data.frame(individual.id = y, location.long = nest[,1], location.lat = nest[,2])
    
    return(nest_loc)
    
  }) %>% reduce(rbind)
  
  return(nests)
    
}) %>% reduce(rbind)



## find the initiation of migration in each year for each individual
#starts <- lapply(list, function(X){})




tracks <- lapply(studies, function(x){
  
  # get the names of the animals to download
  birds <- getMovebankAnimals(study =  x, login = loginStored)# %>% 
    # # reduce the variables
    # dplyr::select(individual_id, timestamp_start, timestamp_end) %>% 
    # filter(timestamp_start != "") %>% 
    # # calculate the time each bird transmitted for
    # mutate(duration = difftime(timestamp_end, timestamp_start, units = "days", tz = "UTC")) %>% 
    # # select only birds that transmitted more than 300 days of data
    # filter(duration > 300)
  birds <- unique(birds$individual_id)
  
  # download the data
  locations <- lapply(birds, function(y){
    
    print(y)
    
    ind_locs <- getMovebankData(study = x, animalName =  y, sensorID = "GPS", 
                                login = loginStored, removeDuplicatedTimestamps = T) %>% 
      # make the move object manipulable
      as.data.frame(row.names = NULL) %>% 
      # remove missing locations
      drop_na(location.long) %>% 
      # thin the data
      mutate(hourly = round_date(timestamp, "1 hour")) %>% 
      group_by(hourly) %>% 
      slice(1) %>% 
      ungroup() %>% 
      group_by(date(timestamp)) %>% 
      # the Haversine distance between first and last locations of the day
      mutate(daily_dist = distVincentyEllipsoid(c(head(location.long,1), head(location.lat, 1)), c(tail(location.long,1), tail(location.lat, 1))),
             # the rhumbline bearing between first and last locations of the day
             daily_direction = bearingRhumb(c(head(location.long,1), head(location.lat, 1)), c(tail(location.long,1), tail(location.lat, 1))),
             # finally a binary category denoting migratory or not, it is unnecessary but simplifies code
             migratory = ifelse(daily_dist > 1e5, 1, 0),
             compass_direction = ifelse(daily_direction > 90 & daily_direction < 270, "southward", "northward"), 
             phase = ifelse(migratory == 1 & compass_direction == "southward", "autumn_migration", 
                            # migrating north 
                            ifelse(migratory == 1 & compass_direction == "northward", "spring_migration", NA))) %>% 
      filter(migratory == 1) %>% 
      ungroup()
    
    return(ind_locs)
    
  }) %>% reduce(rbind)
  
  #save(locations, file = paste0("thinned_migration_multi-year/", x, ".Rdata"))
  
  print(paste0("Downloaded data from study ", which(studies == x), ", ID ", x, "."))
  
  return(locations)
  
}) %>% reduce(rbind)

# read in the saved data and combine into a single df
track_files <- list.files("thinned_migration_multi-year/", full.names = T)

tracks <- lapply(track_files, function(x){
  load(x)
  print(as.character(x))
  if("acceleration.raw.x" %in% names(locations)){
    locations <- locations %>% dplyr::select(-"acceleration.raw.x", -"acceleration.raw.y", -"acceleration.raw.z", -
                                               "barometric.height", -"battery.charge.percent", -"battery.charging.current", -
                                               "external.temperature", -"gps.hdop", -"gps.time.to.fix",  -"height.above.msl", -
                                               "light.level", -"ornitela.transmission.protocol", -
                                               "tag.voltage")
  }
  if("manually.marked.valid" %in% names(locations)){
    locations <- locations %>% dplyr::select(-"manually.marked.valid")
  }

  return(locations)
}) %>% reduce(rbind)

# define birds that took eastern routes as ones that are ever east of 12 longitude
eastern_birds <- unique(tracks$individual.id[tracks$location.long > 12])

# remove the eastern birds and add a phase column
tracks <- tracks %>% 
  filter(!individual.id %in% eastern_birds)%>% 
  mutate(yr_phase = paste0(year(timestamp), "_", phase), 
         track_id = paste0(individual.id, "_", yr_phase))

# number the migrations according to phase
sorting <- tracks  %>% 
  group_by(individual.id, phase) %>% 
  count(yr_phase) %>% 
  mutate(no. = row_number(),
         number = ifelse(no. == 1, "first", 
                         ifelse(no. == 2, "second", 
                                ifelse(no. == 3, "third", 
                                       ifelse(no. == 4, "fourth", 
                                              ifelse(no. == 5, "fifth", 
                                                     ifelse(no. == 6, "sixth", 
                                                            ifelse(no. == 7, "seventh", 
                                                                   ifelse(no. == 8, "eighth", 
                                                                          ifelse(no. == 9, "ninth", "tenth"))))))))),
         migration = paste0(number, "_", sub("_migration", "", phase)))

# add those numbers onto the full data set
tracks <- tracks %>% 
  group_by(individual.id, yr_phase) %>% 
  mutate(migration = sorting$migration[sorting$yr_phase == unique(yr_phase) & sorting$individual.id == unique(individual.id)]) %>% 
  ungroup()

# plot the tracks
countries <- c("Germany", "Switzerland", "France", "Spain", "Morocco", "Italy", "Portugal", "Belgium", "Netherlands",
               "Tunisia", "Algeria", "Liechtenstein", "Luxembourg", "Mauritania", "Mali", "Libya", "Austria", "Slovenia",
               "Czech Republic", "Western Sahara", "Senegal", "Gambia", "Guinea Bissau", "Guinea", "Mali", "Sierra Leone",
               "Liberia", "Ivory Coast", "Togo", "UK", "Ghana", "Burkina Faso", "Benin", "Niger", "Nigeria", "Chad", "Cameroon", "Denmark")


ggplot(tracks, aes(x = location.long, y = location.lat, fill = phase, color = phase)) + 
  borders(region = countries, fill = "gray80") +
  geom_point() +
  #geom_polygon(alpha = 0.4) +
  labs(x = " Longitude", y = "Latitude") +
  #transition_time(time = date)+
  theme_classic() +
  theme(legend.position="none", axis.text = element_text(color = "black"))


# check death dates
last_times <- lapply(studies, function(x){
  times <- getMovebankAnimals(study =  x, login = loginStored) %>% 
    select(individual_id, timestamp_end) %>% 
    distinct()
}) %>% reduce(rbind)

# the last dates of the birds in the data 
last_times <- last_times[last_times$individual_id %in% unique(tracks$individual.id),]

# the last dates during migration (none)
last_times[last_times$timestamp_end %in% tracks$timestamp,]

# build tracks for each bird (amt)
burst_tracks <- lapply(unique(tracks$track_id), function(x){
  
  print(x)
  
  ind <- tracks %>% 
    # select only one track
    filter(track_id == x) %>%
    # make a track object
    make_track(location.long, location.lat, timestamp, track_id) %>% 
    # add the times between reads
    mutate(height.above.ellipsoid = tracks$height.above.ellipsoid[tracks$track_id == unique(track_id)],
           gap = difftime(lead(t_), t_), 
           change = ifelse(round(gap) < 45, 0, 1)) %>%
    group_by(date(t_), change) %>% 
    # add a unique ID to each group
    mutate(burst_ = cur_group_id()) %>% 
    filter_min_n_burst(3) %>% 
    ungroup() %>% 
    select(-gap, -change, -`date(t_)`)
  
  # if(nrow(ind) > 3){
  #   
  #   ssr <- summarize_sampling_rate(ind, time_unit = "min")
  #   
  # } else {ind <- data.frame()}
  # 
  # if(ssr$median < 75 & ssr$median > 45){
  #   ind <- ind   %>%
  #     steps_by_burst(keep_cols = "start")
  # } else {ind <- data.frame()}
  
  if(nrow(ind) > 3){
    ind <- ind %>%
      steps_by_burst(keep_cols = "start") %>% 
      random_steps(n_control = 100) %>%
      rowwise() %>%
      mutate(direction = NCEP.loxodrome.na(y1_, y2_, x1_, x2_), 
             timestamp = paste0(t1_, ".000"))  %>%
      ungroup() %>% 
      rename("location-long" = x2_, 
             "location-lat" = y2_)
  } else {ind <- data.frame()}
  
  return(ind)
  
}) %>% reduce(rbind)

length(unique(gsub("\\_.*", "", burst_tracks$track_id)))

ggplot(burst_tracks, aes(sl_, fill = factor(track_id)))+
  geom_density(alpha = 0.4) +
  labs(x = "Step length (degrees)", y = "Density") +
  theme_classic() +
  theme(legend.position="none", axis.text = element_text(color = "black")) +
  facet_wrap(~case_)

ggplot(burst_tracks, aes(ta_*180/pi, fill = factor(track_id)))+
  geom_density(alpha = 0.4) +
  labs(x = "Turning angle (degrees)", y = "Density") +
  theme_classic() +
  theme(legend.position="none", axis.text = element_text(color = "black")) +
  facet_wrap(~case_)


# centroid_spdf <- SpatialPointsDataFrame(
#   tracks[,c(4,5)], proj4string=srtm@crs, tracks)
# 
# cs <- centroid_spdf[centroid_spdf$location.long > -10 & centroid_spdf$location.long < -5 & centroid_spdf$location.lat > 20 & centroid_spdf$location.lat < 25,]
# cent_max <- raster::extract(srtm,             # raster layer
#                             cs,               # SPDF with centroids for buffer
#                             buffer = 0.25,    # buffer size, units depend on CRS
#                             fun=mean,         # what value to extract
#                             df=TRUE)          # return a data frame?
# unique(cent_max[,2])
# 
# library(geodata)
# srtm <- elevation_3s(lon=6, lat=23,  path=tempdir())
# 
# cent_max <- terra::extract(srtm, tracks[,c(5,4)], fun=NULL, method="simple")

## extract the orthometric height for the given locations

srtm_mosaic <- raster::raster("srtm_mosaic.grd")

centroid_spdf <- SpatialPointsDataFrame(burst_tracks[,c(2,4)], proj4string= srtm_mosaic@crs, burst_tracks)

cent_max <- raster::extract(srtm_mosaic,       # raster layer
                            centroid_spdf,     # SPDF with centroids for buffer
                            buffer = 0.25,    # buffer size, units depend on CRS
                            method = "simple", # whether to return one value or interpolated values
                            fun=mean,         # what value to extract
                            df=TRUE)           # return a data frame?

burst_tracks$srtm <- cent_max$layer
burst_tracks$srtm[is.na(burst_tracks$srtm)] <- 0
burst_tracks$"height-above-msl" <- burst_tracks$height.above.ellipsoid - burst_tracks$srtm
burst_tracks$divider <- c(rep("A", times = 1000000), rep("B", times = 1000000), rep("C", times = 1000000), rep("D", times = 1000000), rep("E", times = 1000000), rep("F", times = nrow(burst_tracks)-5000000))
burst_tracks$divider <- c(rep("A", times = 584921), rep("B", times = 584921), rep("C", times = 584921), rep("D", times = 584921), rep("E", times = 584921), 
                          rep("F", times = 584921), rep("G", times = 584921), rep("H", times = 584921), rep("I", times = 584921), rep("J", times = 584924))
#list2env(split(burst_tracks, burst_tracks$divider), .GlobalEnv)

invisible(lapply(unique(burst_tracks$divider), function(x){
  df <- burst_tracks %>% 
    dplyr::filter(divider == x)# %>% 
    #dplyr::select(2, 4, 5, 6, 7, 8, 13, 14, 15, 16, 19)
  
  write.csv(df, file = paste0("C:/Users/heste/Desktop/HB_storks/tracks/burst_tracks_", x, "_", Sys.Date(), ".csv"))
}))


check2 <- burst_tracks  %>% 
     select(2,4, 7, timestamp)%>% 
     filter(track_id == "1176031140_2020_autumn_migration")
write.csv(check2, file = "check2.csv")



# Specify the data set
request <- list(
  "dataset_short_name" = "reanalysis-era5-single-levels",
  "product_type"   = "reanalysis",
  "variable"       = "boundary_layer_height",
  #"pressure_level" = "925",
  "year"           = "2015",
  "month"          = "09",
  "day"            = "19",
  "time"           = "12:00",
  "area"           = "60/-20/20/0",
  "format"         = "netcdf",
  "target"         = "era5-demo.nc"
)
# Start downloading the data, the path of the file
# will be returned as a variable (ncfile)
ncfile <- wf_request(
  user = "hbronnvik",
  request = request,   
  transfer = TRUE,  
  path = "~",
  verbose = FALSE
)
# Open NetCDF file and plot
library("ncdf4")
r <- raster::raster("adaptor.mars.internal-1662025875.7416215-9183-17-dd109f60-ba66-409c-9159-836ed5e85a5d.nc")
raster::plot(r, main = "ERA-5 Reanalysis (Boundary Layer Height 09.2015)")
maps::map("world", add = TRUE)

fit_distr <- function(x, dist_name, na.rm = TRUE) {
  
  checkmate::check_numeric(x)
  checkmate::check_character(dist_name, len = 1)

  switch(dist_name,
         gamma = {
           if (any(x == 0)) {
             sl_min <- min(x[x !=0])
             x[x == 0] <- sl_min
             base::message(paste0("Steps with length 0 are present. This will lead to an error when fitting a gamma distribution. 0 step lengths are replaced with the smallest non zero step length, which is: ", sl_min))
           }
           
           #  get closed form estimates as starting values for the optimization
           n <- length(x)
           t1 <- n * sum(x * log(x)) - sum(log(x)) * sum(x)
           shape_closed <- (n * sum(x)) / t1
           scale_closed <- 1/n^2 * t1
           
           fit <- fitdistrplus::fitdist(
             x, "gamma", keepdata = FALSE,
             start = list(scale = scale_closed, shape = shape_closed))
           
           make_gamma_distr(
             shape = unname(fit$estimate["shape"]),
             scale = unname(fit$estimate["scale"]),
             vcov = stats::vcov(fit)
           )
         },
         vonmises = {
           xx <- circular::as.circular(
             x, type = "angles", units = "radians", template = "none",
             modulo = "asis", zero = 0, rotation = "counter")
           fit <- circular::mle.vonmises(xx)
           make_vonmises_distr(kappa = fit$kappa, vcov = matrix(fit$se.kappa))
         }
  )
}

# build tracks for the birds (freestyle)
step_tracks <- lapply(unique(tracks$track_id), function(x){
  print(x)
  ind <- tracks %>% 
    filter(track_id == x) %>% 
    mutate(step_length = c(distHaversine(cbind(lag(location.long), lag(location.lat)), cbind(location.long, location.lat))),
           turn_angle = c(bearing(cbind(lag(location.long), lag(location.lat)), cbind(location.long, location.lat)))) 
  
  if(nrow(ind) > 3){
    # generate gamma and von mises distributions
    sl_distr = fit_distr(ind$step_length, "gamma")
    ta_distr = fit_distr(ind$turn_angle, "vonmises")
  } else {  # generate gamma and von mises distributions
    sl_distr = NA
    ta_distr = NA
  }
  
  rand_sl = random_numbers(sl_distr, n = 1e5)
  rand_ta = random_numbers(ta_distr, n = 1e5)

  
}) %>% reduce(rbind)

obs_sl <- burst_tracks %>%
  ggplot(aes(sl_, fill = factor(track_id))) + ggtitle("Observed") +
  geom_density(alpha = 0.4) +
  labs(x = "Step length (degrees)", y = "Density") +
  theme_minimal()
obs_ta <- burst_tracks %>%
  ggplot(aes(ta_, fill = factor(track_id))) + ggtitle("Observed") +
  geom_density(alpha = 0.4) +
  labs(x = "Turn angle (degrees)", y = "Density") +
  theme_minimal()
# rand_sl <- autumn_track[autumn_track$case_ == FALSE,] %>%
#   ggplot(aes(sl_, fill = factor(id))) + ggtitle("Alternative") +
#   geom_density(alpha = 0.4) +
#   labs(x = "Step length (degrees)", y = "Density") +
#   theme_minimal()
# rand_ta <- autumn_track[autumn_track$case_ == FALSE,] %>%
#   ggplot(aes(ta_, fill = factor(id))) + ggtitle("Alternative") +
#   geom_density(alpha = 0.4) +
#   labs(x = "Turn angle (radians)", y = "Density") +
#   theme_minimal()
ggarrange(obs_sl, rand_sl)#, obs_ta, rand_ta, ncol=2, nrow=2, legend = "none")

