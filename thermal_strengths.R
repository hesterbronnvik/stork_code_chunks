### using the onboard data to calculate wind strengths and thermal strengths
### originally written by Dr. Andrea Flack

### hier musst Du die locations anpassen

source('getTrackSegments_updated.R', echo=TRUE)
source('thermallingFeaturesFunction.R')
source('getWindEstimates_update.R', echo=TRUE)
source('MoveApps_download.R', echo=TRUE)
source('SplitGPSBursts.R', echo=TRUE)
#source("getWindEstimate_update.R", echo=TRUE)

require(move)
library(sf)
require(moveWindSpeed)
library(lubridate)

require(parallel)
library(tidyverse)
require(foreach)

# 
setwd("C:/Users/aflack/Desktop/HB_storks")
load("loginStored.RData")

## input name of Movebank study
study = "LifeTrack White Stork SW Germany CASCB"
## download names of animals in study
animals <- getMovebankAnimals(study,login = loginStored)
## remove animals with acceleration sensor - just keep sensor_type 653
animals = animals[-which(animals$sensor_type_id!=653),]
animal_names <- animals$local_identifier

## download first animal of animals variable
a = download.Movebank(study=study, animals=animals[c(1),"local_identifier"],username=username, password=password,thin = FALSE,thin_numb=0,minarg=FALSE)
a <- getMovebankData(study = study, animalName = animal_names[c(37,38,39)], login = loginStored, removeDuplicatedTimestamps = T)
#plot(a)

## segment migrations

migration_threshold <- 1e+05

# downsample to one fix per 15 minutes (significantly speeds the process)
a2 <- a %>% 
  as.data.frame() %>% 
  group_by(local_identifier) %>% 
  mutate(dt_15min = round_date(timestamp, "15 minutes")) %>% 
  group_by(dt_15min) %>% 
  slice(1) %>%  
  ungroup() %>% 
  # add columns by which to categorize the data
  group_by(date(timestamp)) %>% 
  # the Haversine distance between first and last locations of the day
  mutate(daily_dist = distHaversine(c(head(location_long,1), head(location_lat, 1)), c(tail(location_long,1), tail(location_lat, 1))),
         # the rhumbline bearing of between first and last locations of the day
         daily_direction = bearingRhumb(c(head(location_long,1), head(location_lat, 1)), c(tail(location_long,1), tail(location_lat, 1))),
         # finally a binary category denoting migratory or not, it is unnecessary but simplifies code
         migratory = ifelse(daily_dist > migration_threshold, 1, 0),
         compass_direction = ifelse(daily_direction > 90 & daily_direction < 270, "southward", "northward"), 
         phase = ifelse(migratory == 1 & compass_direction == "southward", "autumn_migration", 
                        # migrating north 
                        ifelse(migratory == 1 & compass_direction == "northward", "spring_migration", NA))) %>% 
  ungroup()

migration_dates <- a2 %>% 
  filter(migratory == 1) %>% 
  mutate(journey = ifelse(compass_direction == "southward", "autumn_migration", "spring_migration")) %>% 
  group_by(year(timestamp), journey) %>% 
  summarize(migration_start = head(timestamp, 1),
            migration_end = tail(timestamp, 1)) %>% 
  ungroup() %>% 
  rename("year" = 1) %>% 
  arrange(migration_start)

# define breeding, wintering, and stopover locations
a3 <- a2 %>% 
  group_by(year(timestamp)) %>% 
  mutate(phase = ifelse(!is.na(phase), phase,
                        # earlier than the first fall migration is when it is in or around its natal nest
                        ifelse(timestamp < migration_dates[1,3], "fledging",
                               # fall stopover (between start and end of fall migration in a given year)
                               ifelse(timestamp > migration_dates[migration_dates$year == year(head(timestamp, 1)) & migration_dates$journey == "autumn_migration",3] & timestamp < migration_dates[migration_dates$year == year(head(timestamp, 1)) & migration_dates$journey == "autumn_migration",4], "autumn_stopover",
                                      # spring stopover (between start and end of spring migration in a given year) 
                                      ifelse(timestamp > migration_dates[migration_dates$year == year(head(timestamp, 1)) & migration_dates$journey == "spring_migration",3] & timestamp < migration_dates[migration_dates$year == year(head(timestamp, 1)) & migration_dates$journey == "spring_migration",4], "spring_stopover", 
                                             # breeding (between end of spring and start of fall migration or before start of fall in a given year)
                                             ifelse(timestamp > migration_dates[migration_dates$year == year(head(timestamp, 1)) & migration_dates$journey == "spring_migration",4] & timestamp < migration_dates[migration_dates$year == year(head(timestamp, 1)) & migration_dates$journey == "autumn_migration",3], "breeding",  
                                                    # wintering (between end of fall and start of spring migration in a given year)
                                                    ifelse(timestamp > migration_dates[migration_dates$year == year(head(timestamp, 1)) & migration_dates$journey == "autumn_migration",4], "wintering", 
                                                           ifelse(timestamp < migration_dates[migration_dates$year == year(head(timestamp, 1)) & migration_dates$journey == "spring_migration", 3], "wintering", "no ID"))))))))

a3 <- data.frame()
for (i in unique(year(a2$timestamp))) {
  temp <- a2[which(a2$migratory == 0 & year(a2$timestamp) == i),]
  
  working_dates <- migration_dates[migration_dates$year == i,]
  
  # before the autumn
  temp$phase[temp$timestamp < working_dates$migration_start[working_dates$journey == "autumn_migration"]] <- "breeding"
  
  # during autumn
  temp$phase[temp$timestamp > working_dates$migration_start[working_dates$journey == "autumn_migration"] &
               temp$timestamp < working_dates$migration_end[working_dates$journey == "autumn_migration"]] <- "autumn_stopover"
  # winter
  ifelse("spring_migration" %in% working_dates$journey, 
         temp$phase[temp$timestamp > working_dates$migration_end[working_dates$journey == "autumn_migration"] & 
                    temp$timestamp < working_dates$migration_start[working_dates$journey == "spring_migration"]] <- "wintering",
         temp$phase[temp$timestamp > working_dates$migration_end[working_dates$journey == "autumn_migration"]] <- "wintering")  
  
  # spring
  ifelse("spring_migration" %in% working_dates$journey,
         temp$phase[temp$timestamp > working_dates$migration_start[working_dates$journey == "spring_migration"] & 
                    temp$timestamp < working_dates$migration_end[working_dates$journey == "spring_migration"]] <- "spring_stopover",
         temp$phase <- temp$phase)
  
  # breeding
  ifelse("autumn_migration" %in% working_dates$journey, 
         temp$phase[temp$timestamp > migration_dates$migration_end[migration_dates$journey == "spring_migration"] & 
                    temp$timestamp < migration_dates$migration_start[migration_dates$journey == "autumn_migration"]] <- "breeding",
         temp$phase[temp$timestamp > migration_dates$migration_end[migration_dates$journey == "spring_migration"]] <- "breeding")  

  a3 <- rbind(a3, temp)
}


check <- a3[is.na(a3$phase),]
plot(check$location_long, check$location_lat)

#timestamp > range(a2$timestamp[a2$phase == "autumn_migration"])[1] & timestamp < range(a2$timestamp[a2$phase == "autumn_migration"])[2]
# library(mapview)
# wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")
a2sf <- st_as_sf(a2[is.na(a2$phase),], coords = c("location_long", "location_lat"), crs = wgs)
# mapview(a2sf, zcol = "phase")

## remove data that are not in Bursts
SplitGPSBursts(a,MaxTimeDiff = 10, MinBurstLength = 120)

## 07.05.22 HB attempt to remove non-burst data
#a$time_lag <- unlist(lapply(timeLag(a, units="secs"),  c, NA))

l <- a %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  rename(ID = rowname) %>% 
  mutate(burstID = 1,
         time_lag = difftime(timestamp, lag(timestamp), units = "secs"))
l <- l %>% 
  mutate(burstID = ifelse(is.na(time_lag), 1, 
                          ifelse(time_lag <= 1, dplyr::lag(burstID), dplyr::lag(burstID) + 1)))

for (i in 2:nrow(l)) {
  l$burstID[i] <- ifelse(l$time_lag[i] <= 1, l$burstID[i-1], l$burstID[i-1]+1)
}

## calculate wind estimates
ll <- getWindEstimates(l,windowSize=29)
#x <- lapply(data,getWindEstimates,windowSize=29)

# saveRDS(a,file="animal1.rds")
