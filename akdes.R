### CTMM kernel distributions during migration of white storks
### Hester Bronnvik 
### hbronnvik@ab.mpg.de
### 20.09.2022


library(move)
library(ctmm)
library(lubridate)
library(tidyverse)

## Get the data
# read in the saved data and combine into a single df
track_files <- list.files("C:/Users/heste/Desktop/HB_storks/thinned_migration_multi-year/", full.names = T)

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


setwd("C:/Users/heste/Desktop/HB_storks")
load("loginStored.rdata")

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

# each time step gets its own kde
kde_tracks <- split(thin_tracks, thin_tracks$alignment)
kde_tsub <- kde_tracks[sapply(kde_tracks, function(x) nrow(x) > 30)]
kdes <- lapply(kde_tsub, function(ind){
  print(unique(ind$alignment))
  #ind <- kde_tsub[[1]] # a testing line
  
  ind$individual_id <- "a_name"

  ind <- as.telemetry(ind)#, keep = T)
  
  projection(ind) <-  "ESRI:54009" # mollwiede projektion to keep all the data map-able and comparable
  
  fit <- ctmm.fit(ind) # IID
  KDE1 <- akde(ind, fit, SP = outlines)
  #plot(ind, UD = KDE1, ext = extent(KDE1))
  return(KDE1)
})
#save(kdes, file = "kdes_ls_22.9.RData")
plot(kdes[[200]])

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

