### Finding multi-year birds for SSF analysis
### Hester Bronnvik
### 10.08.2022

library(move)
library(amt)
library(lubridate)
library(sf)
library(tidyverse)

# set wd and load login data
setwd("C:/Users/heste/Desktop/HB_storks")
load("loginStored.rdata")

# collect the names of the stork studies to use
studies <- c(24442409, 212096177, 76367850, 21231406, 1176017658, 173641633)

# load required functions
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
  
  
tracks <- lapply(studies, function(x){
  
  # get the names of the animals to download
  birds <- getMovebankAnimals(study =  x, login = loginStored) %>% 
    # reduce the variables
    dplyr::select(individual_id, timestamp_start, timestamp_end) %>% 
    filter(timestamp_start != "") %>% 
    # calculate the time each bird transmitted for
    mutate(duration = difftime(timestamp_end, timestamp_start, units = "days", tz = "UTC")) %>% 
    # select only birds that transmitted more than 300 days of data
    filter(duration > 300)
  birds <- unique(birds$individual_id)
  
  # download the data
  locations <- lapply(birds, function(y){
    
    print(y)
    
    ind_locs <- getMovebankLocationData(study = x, animalName =  y, sensorID = "GPS", login = loginStored) %>% 
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
      mutate(daily_dist = distHaversine(c(head(location.long,1), head(location.lat, 1)), c(tail(location.long,1), tail(location.lat, 1))),
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
  
  save(locations, file = paste0("thinned_migration_multi-year/", x, ".Rdata"))
  
  print(paste0("Downloaded data from study ", which(studies == x), ", ID ", x, "."))
  
  return(locations)
  
}) %>% reduce(rbind)

# read in the saved data and combine into a single df
track_files <- list.files("thinned_migration_multi-year/", full.names = T)

tracks <- lapply(track_files, function(x){
  load(x)
  if("acceleration.raw.x" %in% names(locations)){
    locations <- locations %>% dplyr::select(-"acceleration.raw.x", -"acceleration.raw.y", -"acceleration.raw.z", -
                                               "barometric.height", -"battery.charge.percent", -"battery.charging.current", -
                                               "external.temperature", -"gps.hdop", -"gps.time.to.fix", -
                                               "height.above.msl", -"light.level", -"ornitela.transmission.protocol", -
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
  mutate(migration = sorting$migration[sorting$yr_phase == yr_phase & sorting$individual.id == individual.id]) %>% 
  ungroup()

# plot the tracks
countries <- c("Germany", "Switzerland", "France", "Spain", "Morocco", "Italy", "Portugal", "Belgium", "Netherlands",
               "Tunisia", "Algeria", "Liechtenstein", "Luxembourg", "Mauritania", "Mali", "Libya", "Austria", "Slovenia",
               "Czech Republic", "Western Sahara", "Senegal", "Gambia", "Guinea Bissau", "Guinea", "Mali", "Sierra Leone",
               "Liberia", "Ivory Coast", "Togo", "UK", "Ghana", "Burkina Faso", "Benin", "Niger", "Nigeria", "Chad", "Cameroon", "Denmark")


ggplot(tracks, aes(x = location.long, y = location.lat, fill = phase, color = phase)) + 
  borders("world", fill = "gray50") +
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
    mutate(gap = difftime(lead(t_), t_), 
           change = ifelse(round(gap) < 45, 0, 1)) %>% 
    # group the days and the gaps
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
      mutate(heading = NCEP.loxodrome.na(y1_, y2_, x1_, x2_), 
             timestamp = paste0(t1_, "000"))  %>%
      rename(location.long = x2_, 
             location.lat = y2_)
  } else {ind <- data.frame()}
  
  return(ind)
  
}) %>% reduce(rbind)

length(unique(gsub("\\_.*", "", burst_tracks$track_id)))

burst_tracks %>% 
  filter(case_ == 1) %>%
  ggplot(aes(sl_, fill = factor(track_id))) + ggtitle("Observed") +
  geom_density(alpha = 0.4) +
  labs(x = "Step length (degrees)", y = "Density") +
  theme_classic() +
  theme(legend.position="none", axis.text = element_text(color = "black"))
burst_tracks %>% 
  filter(case_ == 1) %>%
  ggplot(aes(ta_*180/pi, fill = factor(track_id))) + ggtitle("Observed") +
  geom_density(alpha = 0.4) +
  labs(x = "Turning angle (degrees)", y = "Density") +
  theme_classic() +
  theme(legend.position="none", axis.text = element_text(color = "black"))

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

