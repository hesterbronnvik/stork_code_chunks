### CTMM AKDE building
### Hester Bronnvik 
### hbronnvik@ab.mpg.de
### 2023-01-23

library(move)
library(ctmm)
library(lubridate)
library(terra)
library(tidyverse)

## Get the data
# read in the saved data
locs <- readRDS("C:/Users/hbronnvik/Documents/storkSSFs/migration_locations_40kmCompass_2023-02-15.rds")

# get the track IDs from the migratory birds
rs_ids <- locs %>% 
  group_by(individual.id, track_id) %>% 
  slice(1) %>% 
  dplyr::select(individual.id, track_id) %>% 
  ungroup()

# the total number of migrations attempted and completed by each animal in each season 
meta <- locs %>%
  group_by(individual.id, season) %>% 
  count(track_id) %>% 
  mutate(journey_number = row_number()) %>% 
  ungroup() %>% 
  group_by(individual.id) %>%
  mutate(total_journeys = n()) %>% 
  ungroup()

# add the number of journeys
locs <- locs %>% 
  rowwise() %>% 
  mutate(journey_number = meta$journey_number[which(meta$track_id == track_id)]) %>% 
  ungroup()

# align the locations in time
locs$datestamp <- locs$timestamp
year(locs$datestamp) <- 2024
second(locs$datestamp) <- 00
locs$datestamp <- round_date(locs$datestamp, "hour")

# get the outlines of land to remove the possibility of information over water
# tmax <- raster::getData('worldclim', var = "tmax", res = 10)
# mask <- crop(raster(tmax, 1), raster::extent(min(first_locs$location.long), max(first_locs$location.long), min(first_locs$location.lat), max(first_locs$location.lat)))
# raster::plot(mask)
# # the Straits of Gibraltar are at most 60km across
# mm <- buffer(mask, 30000)
# raster::plot(mm)
# writeRaster(mm, "C:/Users/hbronnvik/Documents/storkSSFs/buffered_water_ras.grd")
mm <- raster::raster("C:/Users/hbronnvik/Documents/storkSSFs/buffered_water_ras.grd")
outlines <- rasterToPolygons(mm, dissolve=TRUE)

# make a list of the unique hours
dates <- split(locs, locs$datestamp)
# keep only list elements with at least 10 data points so that the model can fit
dates <- dates[lapply(dates, nrow) > 10]

# Error in if (ANY) { : missing value where TRUE/FALSE needed
# Error in eigen(hess) : infinite or missing values in 'x'
x <- locs[locs$datestamp == "2024-06-20 06:00:00",] 
# Error in if (ANY) { : missing value where TRUE/FALSE needed
# "2024-06-16 09:00:00 UTC"
# build AKDEs
start_time <- Sys.time()
full_akdes <- lapply(dates, function(x){
  tryCatch({
    print(unique(x$datestamp))
    # create a telemetry object to use the ctmm functions
    ind <- x %>%
      # erase identities because ctmm automatically detects them
      mutate(individual.id == 1,
             individual.local.identifier = 1) %>% 
      # align the locations in space
      as.telemetry(projection = "ESRI:54009")
    # fit a model to the data
    guess <- ctmm.guess(ind, interactive=FALSE)
    fit <- ctmm.select(ind, guess)
    # calculate the KDE using the land outlines and the model
    UD <- akde(ind, fit, SP = outlines, grid = list(dr = 30000))
    return(UD)
    # gc()
  }, error = function(e){print(geterrmessage())})
})
Sys.time() - start_time

ind <- locs %>% 
  arrange(timestamp) %>% 
  mutate(individual.id = as.character(datestamp),
         individual.local.identifier = as.character(datestamp)) %>% 
  group_by(individual.id) %>% 
  mutate(count = n()) %>% 
  ungroup() %>% 
  filter(count >= 10) %>% 
  filter(individual.id %in% unique(individual.id)[1:10])

ind_tele <- as.telemetry(ind, projection = "ESRI:54009")
start_time <- Sys.time()
it <- lapply(ind_tele, function(x){
  tryCatch({
    year(x$timestamp) <- 2024
    guess <- ctmm.guess(x, interactive=FALSE)
    # fit the models
    fit <- ctmm.select(x, guess)
    print(paste0("Fitted a model for ", x$timestamp[1], "."), quote = F)
    pair <- list(x, fit)
    saveRDS(pair, file = paste0("C:/Users/hbronnvik/Documents/storkSSFs/ctmm_fits/", gsub("-|:| ", "",x$timestamp[1]), ".rds"))
    return(fit)
  }, error = function(e){print(geterrmessage(), quote = F)})
})
Sys.time() - start_time
# "2024-06-21 14:30:06."
# Error in if (ANY) { : missing value where TRUE/FALSE needed
# [1] "2024-06-22 06:40:07."
# [1] "infinite or missing values in 'x'"
# [1] "Fitted a model for 2024-06-16 07:40:06."
# [1] "Fitted a model for 2024-07-05 10:30:07."
# [1] "Fitted a model for 2024-07-06 08:30:06."

fitz <- list.files("C:/Users/hbronnvik/Documents/storkSSFs/ctmm_fits/", full.names = T)
kdes <- lapply(1:length(it), function(x){
  fit <- readRDS(fitz[x])
  UD <- akde(fit[[1]], fit[[2]], SP = outlines, grid = list(dr = 30000))
})
# plot(kdes[[1]])


# create a list of hours so that an AKDE can be built for each hour
# create a list of sampling periods so that an AKDE can be built for each one
# first_ls <- first_locs %>% 
#   # temporarily reduce the data to test the function:
#   slice(1:1000) %>% 
#   mutate(day_time = paste(month(timestamp), day(timestamp), sampling_period, sep = "_"))
first_ls <- split(locs[which(locs$journey_number == 1),], locs$datestamp[which(locs$journey_number == 1)])
first_ls <- first_ls[lapply(first_ls, nrow) > 10]
# test 
late_ls <- late_ls[which(names(late_ls) %in% names(first_ls))]

start_time <- Sys.time()
first_akdes <- lapply(first_ls, function(x){
  # create a telemetry object to use the ctmm functions
  ind <- x %>%
    # erase identities because ctmm automatically detects them
    mutate(individual.id == 1,
           individual.local.identifier = 1) %>% 
    # align the locations in space
    as.telemetry(projection = "ESRI:54009")
  # fit a model to the data
  fit <- ctmm.fit(ind)
  # calculate the KDE using the land outlines and the model
  UD <- akde(ind, fit, SP = outlines)
  return(UD)
})

# split the data according to stage
# late_locs <- locs %>% 
#   filter(stage == "adult") %>% 
#   # align the locations in hour
#   mutate(datestamp = round_date(datestamp, unit = "hour"),
#          hour = hour(timestamp),
#          sampling_period = ifelse(hour >= 0 & hour < 6, 1, 
#                                   ifelse(hour >= 6 & hour < 12, 2, 
#                                          ifelse(hour >= 12 & hour < 18, 3, 4))))


# create a list of hours so that an AKDE can be built for each hour
late_ls <- split(late_locs, late_locs$datestamp)
# create a list of sampling periods so that an AKDE can be built for each one
# late_ls <- late_locs %>% 
#   # temporarily reduce the data to test the function:
#   slice(1:1000) %>% 
#   mutate(day_time = paste(month(timestamp), day(timestamp), sampling_period, sep = "_"))
late_ls <- split(locs[which(locs$journey_number != 1),], locs$datestamp[which(locs$journey_number != 1)])
late_ls <- late_ls[lapply(late_ls, nrow) > 10]
# test
late_ls <- late_ls[which(names(late_ls) %in% names(first_ls))]

late_akdes <- lapply(late_ls, function(x){
  # create a telemetry object to use the ctmm functions
  ind <- x %>%
    # erase identities because ctmm automatically detects them
    mutate(individual.id == 1,
           individual.local.identifier = 1) %>% 
    # align the locations in space
    as.telemetry(projection = "ESRI:54009")
  # fit a model to the data
  fit <- ctmm.fit(ind)
  # calculate the KDE using the land outlines and the model
  UD <- akde(ind, fit, SP = outlines)
  return(UD)
})
Sys.time() - start_time

first_akdes <- readRDS("C:/Users/hbronnvik/Documents/storkSSFs/first_akdes_02.02.rds")
late_akdes <- readRDS("C:/Users/hbronnvik/Documents/storkSSFs/late_akdes_02.03.rds")

# combine the adult and juvenile akdes in time
akdes <- lapply(1:length(names(late_akdes)), function(x){
  # the juvenile akde that matches an adult akde
  id <- which(names(first_akdes) == names(late_akdes)[x])
  if(length(id) == 1){# conversion to a raster of the PDF
    f <- raster(first_akdes[[id]], DF = "PDF")
    # using the terra package for reprojection
    f <- terra::rast(f)
    f <- terra::project(f, "EPSG:4326")
    l <- raster(late_akdes[[x]], DF = "PDF")
    l <- terra::rast(l)
    l <- terra::project(l, "EPSG:4326")
    l <- terra::project(l, f)
    # combining the two akdes
    fl <- terra::merge(f, l)
    # adding the date as the name
    names(fl) <- names(first_akdes)[[id]]
    print(paste0("Combined ", names(late_akdes)[x], "."), quotes = F)
    return(fl)}
})
# the hours that have information from both adults and juveniles:
akdes <- akdes[lapply(akdes, is.null) == F]

only_first <- which(!names(first_akdes) %in% names(late_akdes))
only_late <- which(!names(late_akdes) %in% names(first_akdes))

akdes_first <- lapply(1:length(first_akdes[only_first]), function(x){
  # do the same transformations to the juvenile only akdes as done to the combination above
  f <- raster::raster(first_akdes[only_first][[x]], DF = "PDF")
  # using the terra package for reprojection
  f <- terra::rast(f)
  f <- terra::project(f, "EPSG:4326")
  # adding the date as the name
  names(f) <- names(first_akdes[only_first])[x]
  return(f)
})

akdes_late <- lapply(1:length(late_akdes[only_late]), function(x){
  # do the same transformations to the adult only akdes as done to the combination above
  f <- raster::raster(late_akdes[only_late][[x]], DF = "PDF")
  # using the terra package for reprojection
  f <- terra::rast(f)
  f <- terra::project(f, "EPSG:4326")
  # adding the date as the name
  names(f) <- names(late_akdes[only_late])[x]
  return(f)
})

akdes <- c(akdes, akdes_first)#, akdes_late)
# saveRDS(akdes, file = "C:/Users/hbronnvik/Documents/storkSSFs/akdes_combined_2023-02-03.rds")

# plot all of the akdes to see how they move through time  
# p <- vect(outlines)
# png(file="%02d.png", width=480, height=480)
# 
# for (i in 1:length(akdes)){
#   par(mar = c(0, 0, 0, 0))
#   ud <- akdes[[i]]
#   ud <- subst(ud, 0, NA)
#   terra::plot(p, main = names(akdes[[i]]))
#   terra::plot(ud, add = T, buffer = T)
#   terra::plot(p, add = T)
# }
# dev.off()
# 
# library(av)
# library(gtools)
# imgs <- list.files("C:/Users/hbronnvik/Documents/storkSSFs/jan_ras_plot", full.names = TRUE)
# av_encode_video(mixedsort(sort(imgs)), framerate = 115, output = "test2.mp4")

library(terra)
terra::plot(late_akdes[[1]])
projection(late_akdes[[1]])
check <- raster(late_akdes[[1]], DF = "PDF")
check <- rast(check)
check <- terra::project(check, "EPSG:4326")
p <- vect(outlines)
check <- akdes[[6000]]
check <- akdes[[1904]]
check <- subst(check, 0, NA)
terra::plot(p, main = names(akdes[[1904]]))
terra::plot(check, add = T, buffer = T)
terra::plot(p, add = T)

f <- raster(first_akdes[[1]], DF = "PDF")
f <- terra::rast(f)
f <- terra::project(f, "EPSG:4326")
l <- raster(late_akdes[[1]], DF = "PDF")
l <- terra::rast(l)
l <- terra::project(l, "EPSG:4326")
l <- terra::project(l, f)
fl <- terra::merge(f, l)


locs <- locs %>% 
  arrange(track_id)
start_time <- Sys.time()
interpolated_tracks <- lapply(split(locs, locs$track_id), function(x){
  
  # print(paste0("Interpolating the track: ", x, ", ", which(unique(fall_tracks$track_id) == x),
  #              " of ", length(unique(fall_tracks$track_id)), "."))
  print(paste0("Interpolating the track: ", unique(x$track_id), ", ", which(unique(locs$track_id) == unique(x$track_id)),
               " of ", length(unique(locs$track_id)), "."))
  
  # ind <- fall_tracks %>%
  #   filter(track_id == x)
  
  ind <- as.telemetry(x)#, keep = T)
  
  GUESS1 <- ctmm.guess(ind, interactive = FALSE)
  
  print("fitting model")
  
  FIT1_pHREML <- ctmm.select(ind, GUESS1, method = "pHREML", verbose = TRUE)
  
  print("filling the gaps")
  
  filled <- predict(ind, CTMM = FIT1_pHREML[[1]], res = .25)#quarter hour
  
  cn <- colnames(filled)
  
  filled <- as.data.frame(filled@.Data)
  
  colnames(filled) <- cn
  
  filled$track_id <- unique(x$track_id)
  
  return(filled)
}) #%>% reduce(rbind)
Sys.time()-start_time






# pull out autumn data
fall_tracks <- tracks %>% 
  filter(phase == "autumn_migration")  %>% 
  # remove erroneous first fix of the day (only a half hour earlier than the next location)
  group_by(individual.id, year(timestamp), date(timestamp)) %>% 
  slice(-1) %>% 
  ungroup() %>% 
  # get the tracks in order of their ids
  arrange(track_id) %>% 
  group_by(track_id) %>%
  mutate(locs = n()) %>%
  ungroup() %>%
  filter(locs > 99)

# reduce the data to check results
fall_tracks <- fall_tracks %>%
  filter(!track_id %in% unique(fall_tracks$track_id)[101:length(unique(fall_tracks$track_id))])

## interpolate the missing hours in the tracks of each individual

# only 20 birds on their autumn migrations in August and September
locs <- locations_thin %>% 
  filter(month(datestamp) %in% c(8,9))


# fit the akde for each individual
# predict the missing locations using the total time from start to end in hourly intervals
# fit akde for each hour using one location for each individual
start_time <- Sys.time()
interpolated_tracks <- lapply(split(fall_tracks, f = factor(fall_tracks$track_id)), function(x){
  
  # print(paste0("Interpolating the track: ", x, ", ", which(unique(fall_tracks$track_id) == x),
  #              " of ", length(unique(fall_tracks$track_id)), "."))
  print(paste0("Interpolating the track: ", unique(x$track_id), ", ", which(unique(fall_tracks$track_id) == unique(x$track_id)),
               " of ", length(unique(fall_tracks$track_id)), "."))
  
  # ind <- fall_tracks %>%
  #   filter(track_id == x)
  
  ind <- as.telemetry(x)#, keep = T)
  
  GUESS1 <- ctmm.guess(ind, interactive = FALSE)
  
  print("fitting model")
  
  FIT1_pHREML <- ctmm.select(ind, GUESS1, method = "pHREML", verbose = TRUE)
  
  print("filling the gaps")
  
  filled <- predict(ind, CTMM = FIT1_pHREML[[1]], dt = 3600)
  
  cn <- colnames(filled)
  
  filled <- as.data.frame(filled@.Data)
  
  colnames(filled) <- cn
  
  filled$track_id <- unique(x$track_id)
  
  return(filled)
}) #%>% reduce(rbind)
Sys.time()-start_time
#sk <- interpolated_tracks
#save(interpolated_tracks, file = "interpolated_tracks_100inds_21.09.22.RData")

interpolated_tracks$t <- as.POSIXct(interpolated_tracks$t, tz = "UTC", origin = '1970-01-01')

wgs <- "+proj=tpeqd +lat_1=38.4747235276509 +lon_1=-4.10793186111652 +lat_2=45.0136355260082 +lon_2=4.8721292939578 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

i_tracks_mv <- lapply(split(interpolated_tracks, f = interpolated_tracks$track_id), function(df){
  
  mv <- move(x = df$x, y = df$y, time = df$t, proj = wgs, data = df)
  
})

## make stack from filtered data
i_ms <- moveStack(i_tracks_mv, force_tz = TRUE)
#i_ms <- spTransform(fall_ms, center = T)

# # begin with a single id to check run
# ind <- fall_tracks %>% 
#   filter(track_id == unique(fall_tracks$track_id)[1])
# 
# inds <- as.telemetry(ind)
# 
# plot(inds,col=rainbow(length(inds)))
# 
# GUESS1 <- ctmm.guess(inds, interactive = FALSE)
# 
# start_time <- Sys.time()
# FIT1_pHREML <- ctmm.select(inds, GUESS1, method = "pHREML", verbose = TRUE)
# Sys.time() - start_time # Time difference of 1.630057 mins
# summary(FIT1_pHREML)

# SVF <- variogram(inds)
# plot(SVF, CTMM = FIT1_pHREML[[1]],
#      units = TRUE, fraction = 0.5, level = c(0.95, 0.50), 
#      col = "black", col.CTMM = "red")

# AKDE1_pHREML <- akde(inds, FIT1_pHREML, debias = TRUE)
# summary(AKDE1_pHREML, level.UD = 0.95)$CI

# newEXT <- extent(AKDE1_pHREML)
# plot(inds, UD = AKDE1_pHREML, ext = newEXT)
# title(expression("AKDEc"))
# 
# check <- predict(inds, CTMM = FIT1_pHREML[[1]])

# OD <- occurrence(DATA, FITS[[1]])
# SIM <- simulate(DATA, FITS[[1]], dt = 5 %#% "min")
# # should allow an occurrence distribution which sort of smooths over the tracks, then can be 
# # related to underlying environmental variables to see which covariates were sampled by the animal
# # (by multiplying PDF OD (probability of the animal in that cell) with a get value raster layer)
# # over the sampling time regardless of gappy data and not overweighting oversampled locations


## build kdes for each hour
check <- interpolated_tracks 
check$datestamp <- check$t
year(check$datestamp) <- 1995
check$datestamp <- round_date(datestamp)

length(unique(check$datestamp))

ind <- as.telemetry(check)
GUESS <- ctmm.guess(ind, interactive = FALSE)
FIT_pHREML <- ctmm.select(ind, GUESS, method = "pHREML", verbose = TRUE)
AKDE <- akde(ind, FIT_pHREML)








# fresh data
ind_locs <- getMovebankData(study = 24442409, animalName =  24450590, sensorID = "GPS", 
                            login = loginStored, removeDuplicatedTimestamps = T)
unique(year(ind_locs$timestamp))

# the distance covered in each day for the full data (bursts and all)
locs_df <- ind_locs %>% 
  as.data.frame() %>% 
  filter(year(timestamp) == 2015) %>% 
  group_by(date(timestamp)) %>% 
  mutate(distance = distVincentyEllipsoid(c(head(location_long,1), head(location_lat, 1)), 
                                          c(tail(location_long,1), tail(location_lat, 1)))) %>% 
  ungroup()

# the first time after August that it moved at least 100 km . day
on <- locs_df %>% 
  filter(timestamp > as.POSIXct("2015-08-01 00:00:00", tz = "UTC", origin = "1970-01-01") & distance > 100*1000) %>% 
  slice(1) %>% 
  select(timestamp)

# the last time after August that it moved at least 100 km / day
off <- locs_df %>% 
  filter(timestamp > as.POSIXct("2015-08-01 00:00:00", tz = "UTC", origin = "1970-01-01") & distance > 100*1000) %>% 
  slice(n()) %>% 
  select(timestamp)

# filter the data down to between on and off of migration
locs_df <- locs_df %>% 
  filter(timestamp > on & timestamp < off)



load("C:/Users/hbronnvik/Documents/loginStored.rdata")
studies <- c(24442409, 212096177, 76367850, 21231406, 1176017658, 173641633)

## determine the identities of the nestling birds (remove any care center adults)
nestlings <- lapply(studies, function(x){
  
  print(x)
  
  birds <- getMovebankReferenceTable(study = x, login = loginStored) %>% 
    drop_na(animal_id) %>% 
    filter(sensor_type_id == 653)
  
  if(length(unique(birds$animal_id[grep("juv|chick", birds$animal_comments, fixed = F)])) > 0){
    nesties1 <- unique(birds$animal_id[grep("juv|chick", birds$animal_comments, fixed = F)])
  } else {nesties1 <- NA}
  
  if("animal_life_stage" %in% colnames(birds) & length(unique(birds$animal_id[which(birds$animal_life_stage == "nestling")])) > 0){
    nesties2 <- unique(birds$animal_id[which(birds$animal_life_stage == "nestling")])
  } else {nesties2 <-  NA}
  
  nesties <- na.omit(unique(c(unique(nesties1), unique(nesties2))))
  
  return(nesties)
  
}) %>% unlist()


birds <- lapply(studies, function(x){
  birds <- getMovebankAnimals(study =  x, login = loginStored) %>% 
    filter(sensor_type_id == 653) %>% # even birds that die provide information before that (and possibly during)
    dplyr::select(individual_id) %>% 
    mutate(study_id = x)
  return(birds)
}) %>% reduce(rbind)

birds <- birds[which(birds$individual_id %in% nestlings),]

subset <- birds[1:100,]

start_time <- Sys.time()
fall_tracks <- lapply(split(subset, f = subset$individual_id), function(x){
  
  # fresh data
  ind_locs <- getMovebankData(study = x$study_id, animalName =  x$individual_id, sensorID = "GPS", 
                              login = loginStored, removeDuplicatedTimestamps = T)
  # unique(year(ind_locs$timestamp))
  
  # the distance covered in each day for the full data (bursts and all)
  locs_df <- ind_locs %>% 
    as.data.frame() %>% 
    filter(year(timestamp) == unique(year(timestamp))[1]) %>% 
    group_by(date(timestamp)) %>% 
    mutate(distance = distVincentyEllipsoid(c(head(location_long,1), head(location_lat, 1)), 
                                            c(tail(location_long,1), tail(location_lat, 1)))) %>% 
    ungroup()
  
  # the first time after August that it moved at least 100 km/day
  on <- locs_df %>% 
    filter(timestamp > as.POSIXct("2015-07-01 00:00:00", tz = "UTC", origin = "1970-01-01") & distance > 100*1000) %>% 
    slice(1) %>% 
    dplyr::select(timestamp)
  
  # the last time after August that it moved at least 100 km/day
  off <- locs_df %>% 
    filter(timestamp > as.POSIXct("2015-07-01 00:00:00", tz = "UTC", origin = "1970-01-01") & distance > 100*1000) %>% 
    slice(n()) %>% 
    dplyr::select(timestamp)
  
  # filter the data down to between on and off of migration
  locs_df <- locs_df %>% 
    filter(timestamp > on[,1] & timestamp < off[,1]) %>% 
    mutate(seq15 = round_date(timestamp, "15 minutes")) %>% 
    group_by(seq15) %>% 
    slice(1)
  
  return(locs_df)
}) #%>% reduce(rbind)
Sys.time() - start_time
#save(fall_tracks, file = "fall_tracks_byOnOff_15min100ind.RData")
#load("fall_tracks_byOnOff_15min100ind.RData")

fall_tracks <- lapply(fall_tracks, function(x){
  x <- x[, c("location_long","location_lat", "ground_speed", "height_above_ellipsoid",
             "timestamp", "individual_id", "distance", "seq15")]
  
  return(x)
}) %>% reduce(rbind)

# the first location in a given 15 minute interval (rm the bursts)
thin_tracks <- fall_tracks %>% 
  #   mutate(seq15 = round_date(timestamp, "15 minutes")) %>% 
  #   group_by(individual_id, seq15) %>% 
  #   slice(1) %>% 
  #   ungroup() %>% 
  # remove erroneous first fix of the day (only 5 minutes earlier than the next location)
  group_by(individual_id, date(timestamp)) %>% 
  slice(-1)

# define birds that took eastern routes as ones that are ever east of 12 longitude
eastern_birds <- unique(thin_tracks$individual_id[thin_tracks$location_long > 12])

# remove the eastern birds and add a phase column
thin_tracks <- thin_tracks %>% 
  filter(!individual_id %in% eastern_birds)

countries <- c("Niger", "Denmark", "Western Sahara", "Germany", "Spain", "Morocco", "Portugal", "Belgium", "Netherlands",
               "Algeria", "Liechtenstein", "Luxembourg", "Mauritania", "Mali", "Senegal", 
               "Gambia", "Ivory Coast", "Ghana", "Burkina Faso", "Switzerland", "France", "Benin")
ggplot(thin_tracks, aes(x=location_long, y=location_lat)) + 
  borders(region = countries, fill = "gray80") + geom_point()+ facet_wrap(~individual_id)+theme_minimal()

thin_tracks <- thin_tracks %>% 
  mutate(alignment = paste(month(seq15), day(seq15), hour(seq15), minute(seq15), sep = "_"))

thin_tracks$alignment <- thin_tracks$seq15
year(thin_tracks$alignment) <- 1995

check <- thin_tracks %>% 
  group_by(date(alignment)) %>% 
  summarise(rown = n()) %>% 
  group_by(rown) %>% 
  summarise(stamps = n())
hist(check$rown, breaks = nrow(check))


# the number of locations per day (one ID at each time stamp) is much higher than per 15 mins
# this may improve with the addition of the remaining 269 birds?

# fit the akde for each individual
# predict the missing locations using the total time from start to end in hourly intervals
# fit akde for each hour using one location for each individual
thin_tracks <- thin_tracks %>% 
  mutate(track_id = paste(individual_id, year(timestamp), sep = "_")) %>% 
  arrange(track_id)
start_time <- Sys.time()
interpolated_tracks <- lapply(split(thin_tracks, f = factor(thin_tracks$track_id)), function(x){
  
  print(paste0("Interpolating the track: ", unique(x$track_id), ", ", which(unique(thin_tracks$track_id) == unique(x$track_id)),
               " of ", length(unique(thin_tracks$track_id)), "."))
  
  # sdf <- split(thin_tracks, f = factor(thin_tracks$track_id))
  # x <- sdf[[3]]
  
  ind <- as.telemetry(x)#, keep = T)
  
  GUESS1 <- ctmm.guess(ind, interactive = FALSE)
  
  print("fitting model")
  
  FIT1_pHREML <- ctmm.select(ind, GUESS1, method = "pHREML", verbose = TRUE)
  
  print("filling the gaps")
  
  filled <- predict(ind, CTMM = FIT1_pHREML[[1]], dt = 3600)
  
  cn <- colnames(filled)
  
  filled <- as.data.frame(filled@.Data)
  
  colnames(filled) <- cn
  
  filled$track_id <- unique(x$track_id)
  
  return(filled)
}) #%>% reduce(rbind)
Sys.time()-start_time

# install.packages("maptools")
library(maptools)
data(wrld_simpl)
ws <- crop(wrld_simpl, extent(-20,20,0,60))
# plot(ws)
tmax <- raster::getData('worldclim', var = "tmax", res = 10)
mask <- crop(raster(tmax, 1), extent(ws))
raster::plot(mask)
mm <- buffer(mask, 30000)
raster::plot(mm)
# writeRaster(mm, "buffered_water_ras.grd", overwrite = T)
mm <- raster("~/storkSSFs/buffered_water_ras.grd")
outlines <- rasterToPolygons(mm, dissolve=TRUE)

# each time step gets its own kde
thin_tracks <- tracks %>% 
  mutate(hourly = round_date(timestamp, "hour")) %>% 
  group_by(individual.id, hourly) %>% 
  slice(1) %>% 
  ungroup() %>% 
  mutate(alignment = round_date(timestamp, "minute"))
year(thin_tracks$alignment) <- 1995
kde_tracks <- split(thin_tracks, thin_tracks$alignment)
kde_tsub <- kde_tracks[sapply(kde_tracks, function(x) nrow(x) > 30)]
start_time <- Sys.time()
kdes <- lapply(kde_tsub, function(ind){
  print(unique(ind$alignment))
  #ind <- kde_tsub[[1]] # a testing line
  
  ind$individual.id <- "a_name"
  ind$individual.local.identifier <- "a_name"
  
  ind <- as.telemetry(ind)#, keep = T)
  
  projection(ind) <-  "ESRI:54009" # mollwiede projektion to keep all the data map-able and comparable
  
  fit <- ctmm.fit(ind) # IID
  KDE1 <- akde(ind, fit, SP = outlines)
  #plot(ind, UD = KDE1, ext = extent(KDE1))
  return(KDE1)
})
Sys.time() - start_time
#save(kdes, file = "kdes_ls_22.9.RData")
raster::plot(kdes[[200]])

library(rgdal)
library(raster)
library(terra)

srtm_mosaic <- raster::raster("srtm_mosaic.grd")
projection(srtm_mosaic) <- "ESRI:54009"
binary <- srtm_mosaic
values(binary[!is.na(binary)]) <- 1
binary[is.na(srtm_mosaic)] <- 0 
plot(binary) #"+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs"
sp_land <- rasterToPoints(srtm_mosaic, spatial = TRUE)

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
ws <- crop(wrld_simpl, extent(outlines))
regions <- raster(extend(extent(ws), c(1000, 1000)), res = 1000, crs = CRS("+proj=longlat +datum=WGS84 +no_defs"), vals = 0)
#ws@data[,1] <- runif(nrow(ws))
plot(ws, col = ws$REGION)
check <- raster::rasterize(ws, regions, field = "FIPS", background = 0, fun = "first")
plot(check)


tmax <- raster::getData('worldclim', var = "tmax", res = 10)
mask <- crop(raster(tmax, 1), extent(ws))
raster::plot(mask)
mm <- buffer(mask, 60000)
raster::plot(mm)
#writeRaster(mm, "C:/Users/heste/Desktop/HB_storks/buffered_water_ras.grd")

outlines <- rasterToPolygons(mm, dissolve=TRUE)



uds <- lapply(kdes, function(x){
  hr <- ctmm::SpatialPolygonsDataFrame.UD(x, proj4string = datproj) %>%
    sp::spTransform(sp::CRS("+proj=longlat")) %>%
    fortify() %>%
    bind_rows()
  return(hr)
})

# plug in the data and the column names
ggplot(uds, aes(x = location.long, y = location.lat)) + 
  # if you want borders, change the fill shade
  borders(regions = countries, fill = "gray50") +
  # this is an alternative to listing all the countries, but I find it behaves oddly:
  #borders(database = "world", xlim = c(-10, 20), ylim = c(10, 60), fill = "gray50") + 
  # I don't know whether you would want to color by ID, year, season, etc.
  geom_polygon(aes(color = year(timestamp))) +  
  # optional to change colors, but requires the viridis library:
  #scale_fill_viridis(discrete=TRUE, option="A") + 
  # fix the labels
  labs(x = " Longitude", y = "Latitude")+
  # clean up the gray background
  theme_classic() +
  # remove the legend
  theme(legend.position="none", axis.text = element_text(color = "black")) # +
# and if you want to animate that:
#transition_time(time = timestamp)
#transition_states(states = year(timestamp))

datproj <- sp::CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs")
# Build images -> save them at .png format
png(file="%02d.png", width=480, height=480)

for (i in 1:length(kdes)){
  hr <- ctmm::SpatialPolygonsDataFrame.UD(kdes[[i]], proj4string = datproj) %>%
    sp::spTransform(sp::CRS("+proj=longlat"))
  #hr <- uds[[i]]
  par(mar = c(0, 0, 0, 0))  
  sp::plot(flyway, col="grey80", border= F, ylim = c(-11,65)) 
  sp::plot(hr, add=T)
  #mtext(paste(day(dates[i]), "-", month(dates[i])), side=3)
}
dev.off()

library(magick)
imgs <- list.files("C:/Users/heste/Desktop/HB_storks/kdegif/", full.names = TRUE)
img_list <- lapply(gtools::mixedsort(imgs), image_read)

## join the images together
img_joined <- image_join(img_list)

## animate at 2 frames per second
img_animated <- image_animate(img_joined, fps = 100)

image_write(image = img_animated,
            path = "C:/Users/heste/Desktop/HB_storks/kernels_ctmm_15min.gif")




Lines_ws <- lapply(split(fall_tracks,fall_tracks$track_id), function(x){
  points <- sp::SpatialPoints(cbind(x$long, x$lat), proj4string = wgs)
  sp_line <- as(points,"SpatialLines")
})


# bltr
par(mar = c(0, 0, 0, 0))  
sp::plot(flyway, col="grey70", border= F, ylim = c(-11,65)) 

#add latitude lines
lines(x = c(-20,45.9), y = c(0,0),lty = 2,lwd = 1, col = "grey50")
lines(x = c(-20,45.9), y = c(30,30),lty = 2,lwd = 1, col = "grey50")
lines(x = c(-20,45.9), y = c(60,60),lty = 2,lwd = 1, col = "grey50")
text(x = -17, y = 30.7, "30° N", col = "grey50", cex = 1.1, font = 3)
text(x = -17, y = 60.7, "60° N", col = "grey50", cex = 1.1, font = 3)
text(x = -17, y = 0.70, "0° N", col = "grey50", cex = 1.1, font = 3)

invisible(lapply(1:length(Lines_ws), function(x){
  l <- Lines_hb[[x]]
  lines(l, lty= 1, lwd =  2.5, col= rainbow)
}))

### 11.05.2022
library(maptools)
data("wrld_simpl")
ws <- crop(wrld_simpl, extent(-20,20,0,60))
raster::plot(ws)
raster::plot(raster(kdes[[1]]))
ws <- spTransform(ws, crs(raster(kdes[[1]])))
raster::plot(ws, add = T)

build <- raster(kdes[[1]])
build[build > 0.95] <- NA
raster::plot(build, col = heat.colors(8, alpha = 1), legend = F)
raster::plot(ws, add = T)

tracks <- tracks %>% 
  mutate(datestamp = timestamp) 

year(tracks$datestamp) <- 1995

kde_tsub[[1]]$alignment[1]

tsub_stamps <- lapply(1:length(kde_tsub), function(x){
  stamp <- kde_tsub[[x]]$alignment[1]
  return(stamp)
}) %>% 
  unlist() %>% 
  as.POSIXct(tz = "UTC", origin = "1970-01-01")


bt1 <- read.csv("bt1.csv")
bt2 <- read.csv("bt2.csv")
bt3 <- read.csv("bt3.csv")

burst_tracks <- rbind(bt1, bt2, bt3)

burst_tracks <- burst_tracks %>% 
  mutate(timestamp = as.POSIXct(gsub("\\.000", "", timestamp), tz = "UTC"), 
         alignment = round_date(timestamp, "minute"),
         location.long = as.numeric(location.long),
         location.lat = as.numeric(location.lat))
year(burst_tracks$alignment) <- 1995

t <- lapply(tsub_stamps, function(x){
  # print(x)
  df <- burst_tracks %>% 
    filter(alignment == x)
  
  if(nrow(df > 3)){
    # for each unique hour in the data, call the list item
    id <- which(tsub_stamps == x)
    kd <- raster(kdes[[id]])
    
    # make the data an SP object and transform to the same CRS as the AKDE
    spdf <- df
    coordinates(spdf) <- ~location.long + location.lat
    crs(spdf) <- crs("+proj=longlat +datum=WGS84 +no_defs +type=crs")
    spdf <- spTransform(spdf, crs(kd))
    
    # extract the value for each location
    vals <- raster::extract(kd, spdf)#[, c("location.long", "location.lat")])
    
    # append
    df <- df %>% 
      mutate(UD = vals)
    return(df)
  }
  
  
}) %>% reduce(rbind)

t <- t %>% 
  mutate(season = str_extract(track_id, "fall|spring"),
         individual_id = str_extract(track_id, "[^_]+"))

# number the migrations according to phase
sorting <- t  %>% 
  group_by(individual_id, season) %>% 
  count(track_id) %>% 
  mutate(no. = row_number(),
         number = ifelse(no. == 1, "first", 
                         ifelse(no. == 2, "second", 
                                ifelse(no. == 3, "third", 
                                       ifelse(no. == 4, "fourth", 
                                              ifelse(no. == 5, "fifth", 
                                                     ifelse(no. == 6, "sixth", 
                                                            ifelse(no. == 7, "seventh", 
                                                                   ifelse(no. == 8, "eighth", 
                                                                          ifelse(no. == 9, "ninth", "tenth")))))))))) %>% 
  ungroup()

# add those numbers onto the full data set
t <- full_join(t, sorting) %>% 
  mutate(migration = paste0(number, "_", season)) 
t$number <- factor(t$number, levels = c("first", "second", 'third', "fourth", "fifth", "sixth", "seventh", "eighth", "ninth"))

ggplot(t, aes(number, UD)) + 
  geom_boxplot() +
  theme_classic() +
  facet_wrap(~case_)

library(survival)
modelUD <- function(df) {
  clogit(case_ ~ scaled_UD + scaled_sl + scaled_ta + strata(step_id_), data = df)
}

SSF_results <- t %>%
  mutate(scaled_UD = scale(UD),
         scaled_sl = scale(sl_),
         scaled_ta = scale(ta_)) %>% 
  group_by(track_id) %>%
  filter(length(unique(UD)) > 3) %>% 
  nest() %>%
  mutate(ssf_modelUD = purrr::map(data, modelUD),
         ssf_coefsUD = purrr::map(ssf_modelUD, coef),
         AIC_TV = map_dbl(ssf_modelUD, ~AIC(.)))

# flatten the coefficient column
ssf_coefs <- unnest(SSF_results, ssf_coefsUD) %>% 
  group_by(track_id) %>% 
  # take the UD coefs
  slice(1) %>% 
  ungroup() %>% 
  mutate(individual_id = str_extract(track_id, "[^_]+")) %>% 
  full_join(sorting) %>% 
  group_by(number) %>% 
  mutate(samp = length(unique(track_id))) %>% 
  ungroup() %>% 
  mutate(migration = paste(number, season, sep = "_"))

library(EnvStats)
ggplot(ssf_coefs[which(abs(ssf_coefs$ssf_coefsUD) < 40),], aes(as.factor(migration), ssf_coefsUD, group = migration, fill = season, label = number)) +
  geom_boxplot() +
  labs(x = "Migration", y = "Coefficient") +
  theme_classic() +
  theme(text = element_text(size = 20),
        axis.text = element_text(colour = "black"))  +
  stat_n_text()

