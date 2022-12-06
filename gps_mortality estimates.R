### Using recursions to estimate mortality dates of White Storks
### Hester Bronnvik
### 30.11.2022

# The goal is simply to draw a radius around the first point in the data and 
# ask when the bird arrives in this radius for good. If that date is a threshold number of days
# from the last location or from now, the bird is lost.

library(recurse)
library(move)
library(mapview)
library(sf)
library(lubridate)
library(tidyverse)

r_thresh <- 5 # number of days that the animal is inside the radius to count as dead
studies <- c(24442409, 212096177, 76367850, 21231406, 1176017658, 173641633, 908232414)
load("C:/Users/hbronnvik/Documents/loginStored.rdata")

# the death dates of the birds for which ACC estimates functioned
done <- readRDS("C:/Users/hbronnvik/Documents/storkSSFs/estimated_death_dates_2022-11-12.rds") %>% 
  # drop_na(death_date) %>% 
  select(local_identifier) %>% 
  deframe() 

# birds that migrated (50 km segmentation) named found_birds
load("C:/Users/hbronnvik/Documents/storkSSFs/clean_data/reduced50_15m_2022-10-27.RData")
migratory <- found_birds %>% 
  ungroup() %>% 
  select(individual.local.identifier) %>% 
  rename(local_identifier = individual.local.identifier) %>% 
  deframe() %>% 
  unique()

# the data for the birds in all studies
info <- lapply(studies, function(x){
  i <- getMovebankAnimals(study = x , login = loginStored) %>%
    filter(sensor_type_id == 2365683) %>% 
    mutate(study = x)
  return(i)
}) %>% reduce(rbind)

# the names of the birds that have not yet been estimated
info <- info %>% 
  filter(!local_identifier %in% done & local_identifier %in% migratory)

# the recurse approach: works well for birds that lie undisturbed at their sites of death until 
# end deployment, but fails on birds that move post-mortem or die suddenly
# thus, added a filter for number of days in one location and large gaps
death_estimates_gps <- lapply(1:nrow(info), function(x){
  print(info$local_identifier[x], quote = F)
  # the data for one bird's last 6 months. This means only 1 of breeding or fledging can be present.
  mv <- getMovebankData(info$study[x], info$local_identifier[x], loginStored, T,
                        timestamp_start = paste0(gsub("[[:punct:]]| ", "",
                                                      (as.POSIXct(sub("\\.000", "", info$timestamp_end[x]), tz = "UTC") - 6*30*24*3600)), "000"))
  # select only locations at a lower than 5 minute resolution (reduce memory burden on recursions)
  mv$lag <- c(NA, move::timeLag(mv, units = "mins"))
  mv <- mv[which(mv$lag > 5),]
  # check for large gaps in transmission
  dd <- data.frame(dates = mv$timestamp) %>% # difference between this and last transmission
    mutate(change = difftime(dates, lag(dates), units = "days"))
  # are the data gappy?
  if(max(na.omit(dd$change)) > 60){
    lost_loc <- data.frame(location_lat = as.numeric(mv$location_lat[timestamps(mv) == dd$dates[which(dd$change > 60)-1]]),
                           location_long = as.numeric(mv$location_long[timestamps(mv) == dd$dates[which(dd$change > 60)-1]]))
    found_loc <- data.frame(location_lat = as.numeric(mv$location_lat[timestamps(mv) == dd$dates[which(dd$change > 60)]]),
                            location_long = as.numeric(mv$location_long[timestamps(mv) == dd$dates[which(dd$change > 60)]]))
    # take the last location before the gap if the tag was moved afterwards
    if(distVincentyEllipsoid(lost_loc, found_loc) > 100) {
      loss <- as.numeric(difftime(tail(mv$timestamp, 1), dd$dates[which(dd$change > 60)-1], units = "days")) > r_thresh
      dod <- data.frame(local_identifier = mv@idData$local_identifier, dod = dd$dates[which(dd$change > 60)-1], loss = loss, comment = "gap")
    # but take the recurse approach if the animal did not move during the gap
      }else{rec <- getRecursions(mv, 0.003, timeunits = "days")
    # after the animal starts moving (discard time in/on a nest)
    minimum <- rec$revisitStats$entranceTime[which.min(rec$revisitStats$timeInside)]
    # add a threshold for amount of residence time that can be alive
    revists <- rec$revisitStats %>% 
      filter(entranceTime > minimum & timeInside > 7)
    # save the date of entrance to the radius, or if there was none, no evidence of death 
    if(nrow(revists) > 0){
      dod <- revists$entranceTime[which.max(revists$timeInside)]
      loss <- as.numeric(difftime(tail(mv$timestamp, 1), dod, units = "days")) > r_thresh
      dod <- data.frame(local_identifier = mv@idData$local_identifier, dod = dod, loss = loss, comment = "death")
    }else(dod <- data.frame(local_identifier = mv@idData$local_identifier, dod = NA, loss = F, comment = "no death"))}    
  # run the recurse approach on animals that do not have gappy data
  }else{rec <- getRecursions(mv, 0.003, timeunits = "days")
  minimum <- rec$revisitStats$entranceTime[which.min(rec$revisitStats$timeInside)]
  revists <- rec$revisitStats %>% 
    filter(entranceTime > minimum & timeInside > 7)
  if(nrow(revists) > 0){
    dod <- revists$entranceTime[which.max(revists$timeInside)]
    loss <- as.numeric(difftime(tail(mv$timestamp, 1), dod, units = "days")) > r_thresh
    dod <- data.frame(local_identifier = mv@idData$local_identifier, dod = dod, loss = loss, comment = "death")
  }else(dod <- data.frame(local_identifier = mv@idData$local_identifier, dod = NA, loss = F, comment = "no death"))}
  return(dod)
  gc()
}) %>% reduce(rbind)

ggplot(rec$revisitStats, aes(entranceTime, timeInside)) +
  geom_point() + 
  geom_segment(aes(x=entranceTime, xend=entranceTime, y=0, yend=timeInside)) +
  geom_point(data = revists[which.max(revists$timeInside),], color = "red") + 
  geom_segment(data = revists[which.max(revists$timeInside),], 
               aes(x=entranceTime, xend=entranceTime, y=0, yend=timeInside),
               color = "red") +
  theme_classic()

# the ground speed approach
build <- lapply(build_birds, function(x){
  # the data for one bird
  mv <- x#getMovebankData(info$study[x], info$local_identifier[x], loginStored, T)
  df <- mv %>% 
    as.data.frame() %>% 
    arrange(timestamp) %>% 
    mutate(new_event = ifelse(round(ground_speed) <= 10, F, T),
           # take the cumulative sum to act as a unique ID for each burst identified in line 37
           cumu_check_for_event = cumsum(new_event)) %>% 
    group_by(cumu_check_for_event) %>% 
    # for the bursts, calculate the time difference between the last and first locations
    mutate(event_length = ifelse(n()>1, difftime(tail(timestamp,1), head(timestamp, 1), units = "secs"), NA),
           # add an ID to each burst
           event_ID = cur_group_id())
  
  check <- df[df$event_ID == unique(df$event_ID)[length(unique(df$event_ID))],]
  wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")
  check <- st_as_sf(check, coords = c("location_long", "location_lat"), crs = wgs)
  mapview(check, zcol = "ground_speed")
  
  dod <- check$timestamp[1]
  loss <- as.numeric(difftime(tail(mv$timestamp, 1), dod, units = "days")) > r_thresh
  dod <- data.frame(local_identifier = mv@idData$local_identifier, dod = dod, loss = loss)
  return(dod)
}) %>% reduce(rbind)

### Visually determine the death dates of example birds for a "truth" to compare estimation to
mapview(build_birds[[10]], zcol = "date")

vis_dod <- data.frame(ID = c(build_birds[[1]]@idData$local_identifier,
                             build_birds[[2]]@idData$local_identifier,
                             build_birds[[3]]@idData$local_identifier,
                             build_birds[[4]]@idData$local_identifier,
                             build_birds[[5]]@idData$local_identifier,
                             build_birds[[6]]@idData$local_identifier,
                             build_birds[[7]]@idData$local_identifier,
                             build_birds[[8]]@idData$local_identifier,
                             build_birds[[9]]@idData$local_identifier,
                             build_birds[[10]]@idData$local_identifier),
                      dod = c("2014-12-01",
                              "2015-04-17",
                              NA, # still alive
                              "2016-02-27", # disappeared
                              "2015-10-15",
                              "2016-06-10", # Sudan bird transported after disappearance
                              "2020-04-24",
                              "2017-01-20",
                              "2016-08-05",
                              "2016-12-09"), # disappeared
                      comment = c("dead", "dead", "alive", "disappeared", "dead", 
                                  "disturbed", "dead", "dead", "dead", "disappeared")) 

check <- build_birds[[16]]
check <- check[timestamps(check) > "2021-03-07"]
mapView(check, zcol = "date")

# check the next ten birds as a test of the approach
vis_dod2 <- data.frame(ID = c(build_birds[[11]]@idData$local_identifier,
                             build_birds[[12]]@idData$local_identifier,
                             build_birds[[13]]@idData$local_identifier,
                             build_birds[[14]]@idData$local_identifier,
                             build_birds[[15]]@idData$local_identifier,
                             build_birds[[16]]@idData$local_identifier,
                             build_birds[[17]]@idData$local_identifier,
                             build_birds[[18]]@idData$local_identifier,
                             build_birds[[19]]@idData$local_identifier,
                             build_birds[[20]]@idData$local_identifier),
                      dod = c(NA, # still alive
                              NA, # still alive
                              "2020-03-24", # disappeared, moved along a highway
                              NA,
                              "2020-10-15",
                              "2021-03-14", # lost near buildings, reappeared in a home a month later
                              "2014-09-03", # on migration
                              "2015-09-22",
                              "2015-08-03", # lost on a highway near Zaragoza
                              NA), # alive?
                      comment = c("alive", "alive", "disturbed", "alive", "dead", 
                                  "disturbed", "disappeared", "dead", "dead", "alive")) 
vis_dod <- rbind(vis_dod, vis_dod2)

# the new estimator is in agreement with visualization on all 20 birds (+/- 1 day)
# given NAs for birds that disappeared or that are still alive

checkers <- death_estimates_gps %>% 
  filter(comment == "gap" & local_identifier != "Sophie + / DER AT898 (eobs 2664)") %>% 
  select(local_identifier) %>% 
  deframe()

build_birds <- lapply(41:nrow(info), function(x){
  # x <- which(info$local_identifier == ind)
  # the data for one bird's last six months
  mv <- getMovebankData(info$study[x], info$local_identifier[x], loginStored, T, 
                        timestamp_start = paste0(gsub("[[:punct:]]| ", "", 
                                                      (as.POSIXct(sub("\\.000", "", info$timestamp_end[x]), tz = "UTC") - 8*30*24*3600)), "000"))
  mv$lag <- c(NA, move::timeLag(mv, units = "mins"))
  # select only locations at a higher than 5 minute resolution
  mv <- mv[which(mv$lag > 5),]
  mv$date <- as.character(date(mv$timestamp))
  return(mv)
})

mapview(build_birds[[13]], zcol = "date")
build_birds[[13]]@idData$local_identifier
mv <- build_birds[[11]]
mv <- mv[timestamps(mv) > "2021-12-18"]
mapView(mv, zcol = "date")

# "Kiki + / DER AN872 (eobs 3648)" died by a highway, disappeared for 4 months, and reappeared in the same spot
# "Eamy + / DER AL585 (eobs3063)" vanished moving north from Gibraltar, appeared 4 months later in France, and then a month later in Spain
# "Petra + A9100 (eobs 8016)" lay in the same place for the whole 6 months

# metadata for corrections to the death dates made on the basis of visual examination of maps (2022-12-06 HB)
dod_corrections <- data.frame(local_identifier = c("Isolde + / DER AU641 (eobs 2760)",
                                                   "Mattis +/ DER AW839 (eobs 3041)",
                                                   "Mira / DER AU650 (eobs 3027)",
                                                   "Borni II + / DER AX351 (eobs 3076)",
                                                   "Nicole + / DER AX555 (eobs 3020)",
                                                   "Maxi2 + / DER AZ919 (e-obs 6586)",
                                                   "Niclas / DER AU053 (eobs 3341)",
                                                   "Tobi + / DER AU052 (eobs 3340)",
                                                   "Redrunner + / DER AU057 (eobs 3339)",
                                                   "Sierit-chick + / DER AZ347 (eobs 4565)",
                                                   "Petra + A9100 (eobs 8016)"), 
                              dod = c("2014-08-28 18:46:08",
                                      "2016-07-07 15:06:54",
                                      NA,
                                      "2018-09-18 18:50:47",
                                      "2018-10-30 09:49:18",
                                      "2020-11-08 08:35:49",
                                      NA,
                                      NA,
                                      "2019-03-27 23:59:58",
                                      "2016-09-12 12:59:03",
                                      "2021-03-21 15:30:08"), 
                              loss = c(T, T, T, T, T, T, T, T, T, T, T), 
                              comment = c("death: misestimation due to rapid end_deployment and long nesting bouts",
                                          "death: misestimation due to residence in the last location for 7 months",
                                          "disappearance: misestimation due to long stay on a Cordoba landfill",
                                          "death: misestimation due to rapid tag recovery",
                                          "death: misestimation due to rapid tag recovery",
                                          "death: misestimation due to 21 day gap",
                                          "disappearance: misestimation due to breeding (?)",
                                          "disappearance: misestimation due to long stay on a Kenitra landfill",
                                          "death: misestimation due to long stay on a Fez mechouar",
                                          "death: misestimation due to long residence in/on the nest",
                                          "death: misestimation due to residence in the last location for a year and GPS reflection"))

# add the corrections to the estimates
death_estimates_gps <- death_estimates_gps %>% 
  filter(!local_identifier %in% dod_corrections) %>% 
  rbind(dod_corrections)
