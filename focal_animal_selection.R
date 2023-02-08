### Bird selection
### 2023-01-20
### Hester Bronnvik

### The goal here is to select the birds to be used for route selection analyses on the basis of 
### having completed migrations.
library(move)
library(lubridate)
library(tidyverse)

s_thresh <- 6*7 # the number of days a bird has to have survived since its last migratory day to have completed migration
studies <- c(24442409, 212096177, 76367850, 21231406, 1176017658, 173641633)#, 908232414)

# clean_files <- list.files("C:/Users/hbronnvik/Documents/storkSSFs/clean_data/", pattern = "embc", full.names = T)
visDODs <- readRDS("C:/Users/hbronnvik/Documents/storkSSFs/visual_ACC_deaths_2023-02-01.rds")
load("C:/Users/hbronnvik/Documents/loginStored.rdata")

# take the threshold number of days before death
# are they migratory?
# filter and save the final IDs

# cleaned_files <- list.files("C:/Users/hbronnvik/Documents/storkSSFs/classified_gps/", full.names = T)
# cleaned_files <- cleaned_files[sapply(cleaned_files, file.size) > 10000]
# migration_locations <- lapply(1:length(cleaned_files), function(x){
#   ml <- readRDS(cleaned_files[x])
#   ml <- ml %>% 
#     dplyr::select(individual.id, individual.local.identifier, timestamp, location.lat, location.long, phase, track_id, daily_dist, daily_direction, gr_speed)
#   return(ml)
# }) %>% reduce(rbind)

# the full data
# migration_locations <- readRDS("C:/Users/hbronnvik/Documents/storkSSFs/classified_gps/full_migration_locations_40km_2020_2023-02-07.rds")
# # down-sample to 15 minute locations
# migration_locations <- migration_locations%>% 
#   mutate(seq_15 = round_date(timestamp, unit = "15 minutes")) %>% 
#   group_by(individual.id, seq_15) %>% 
#   slice(1) %>% 
#   ungroup()
# saveRDS(migration_locations, file = "C:/Users/hbronnvik/Documents/storkSSFs/classified_gps/full_migration_locations_40km15min_2020_2023-02-07.rds")

# the data down-sampled to 15 minutes, with migration defined as moving more than 40 km/day, 
# pre-breeding as Dec. to the next July, and post-breeding as Aug. to Nov.
migration_locations <- readRDS("C:/Users/hbronnvik/Documents/storkSSFs/classified_gps/full_migration_locations_40km15min_2020_2023-02-07.rds")

# make sure that removing low file sizes removed missing tracks
migration_locations <- migration_locations %>% 
  drop_na(track_id)
# define birds that took eastern routes as ones that are ever east of 15.7 longitude (East Germany)
eastern_birds <- unique(migration_locations$individual.id[migration_locations$location.long > 15.7])
# remove the eastern birds
migration_locations <- migration_locations %>% 
  filter(!individual.id %in% eastern_birds)
# get the track IDs from the migratory birds
rs_ids <- migration_locations %>% 
  group_by(individual.id, track_id) %>% 
  slice(1) %>% 
  dplyr::select(individual.id, track_id) %>% 
  ungroup()
# the total number of migrations attempted and completed by each animal in each season 
meta <- migration_locations %>%
  mutate(season = ifelse(grepl("post", phase), "post-breeding", "pre-breeding")) %>% 
  group_by(individual.id, season) %>% 
  count(track_id) %>% 
  mutate(journey_number = row_number()) %>% 
  ungroup() %>% 
  group_by(individual.id) %>%
  mutate(total_journeys = n()) %>% 
  ungroup()
# add the number of journeys
migration_locations <- migration_locations %>% 
  rowwise() %>% 
  mutate(journey_number = meta$journey_number[which(meta$track_id == track_id)]) %>% 
  ungroup()

lost <- migration_locations %>% 
  group_by(individual.id) %>% 
  slice(n()) %>% 
  ungroup() %>% 
  dplyr::select(individual.id, individual.local.identifier) %>% 
  filter(!individual.id %in% visDODs$individual_id)

# info <- lapply(studies, function(x){
#   info <- getMovebankAnimals(x, loginStored) %>% 
#     filter(sensor_type_id == 653 & individual_id %in% unique(lost$individual.id))
#   return(info)
# }) %>% reduce(rbind)

found <- data.frame(local_identifier = lost$individual.local.identifier,
                    dod = c("2015-08-08 02:23:00", "2016-06-14 08:56:00", NA),
                    loss = c(T, T, NA),
                    comment = c("classified by acc", "classified by acc", "classified by acc"),
                    individual_id = c(80671851, 1173982416, 1176059034),
                    study = c(76367850, 24442409, 1176017658))
visDODs <- rbind(visDODs, found)

# take the last time stamps from the full GPS data for each animal
info <- lapply(studies, function(x){
  i <- getMovebankAnimals(study = x , login = loginStored) %>%
    filter(sensor_type_id == 653 & individual_id %in% unique(migration_locations$individual.id)) %>% 
    mutate(study = x)
  return(i)
}) %>% reduce(rbind) %>% 
  mutate(timestamp_end = as.POSIXct(timestamp_end, tz = "UTC"))

# find the success/failure of each migration track
migration_locations <- lapply(split(migration_locations, migration_locations$individual.id), function(x){
  # take the last migratory route
  final_track <- x %>% 
    arrange(timestamp) %>% 
    filter(track_id == track_id[n()])
  # if the last GPS time point is greater than the time the data were downloaded, use the download date
  t_end <- ifelse(info$timestamp_end[info$individual_id == unique(x$individual.id)] > "2023-01-28", 
                  "2023-01-28 00:00:00", info$timestamp_end[info$individual_id == unique(x$individual.id)])
  # take either the confirmed DOD or the last transmitted time stamp
  loss_time <- visDODs %>% 
    filter(individual_id == unique(x$individual.id)) %>% 
    mutate(dod = as.POSIXct(ifelse(is.na(dod), t_end, dod), tz = "UTC", origin = "1970-01-01")) %>% 
    dplyr::select(dod) %>% 
    deframe()
  # compare the DOD to the end of the migration
  loss <- max(final_track$timestamp) > loss_time - days(s_thresh)
  # add a column containing the outcome of the migratory track
  x <- x %>% 
    group_by(track_id) %>% 
    mutate(track_status = ifelse(loss == T & track_id == unique(final_track$track_id), "incomplete", "complete")) %>% 
    ungroup()
  return(x)
}) %>% reduce(rbind)

complete_ml <- migration_locations %>% 
  filter(track_status == "complete")

# the total number of migrations attempted and completed by each animal in each season 
meta <- migration_locations %>%
  mutate(season = ifelse(grepl("post", track_id), "post-breeding", "pre-breeding")) %>% 
  group_by(individual.id, season) %>% 
  count(track_id) %>% 
  mutate(journey_number = row_number()) %>% 
  ungroup() %>% 
  group_by(individual.id) %>%
  mutate(total_journeys = n()) %>% 
  ungroup() %>% 
  rowwise() %>% 
  mutate(track_status = unique(migration_locations$track_status[which(migration_locations$track_id == track_id)])) %>% 
  ungroup()

# determine the final tracks to use
rs_ids <- meta %>% 
  filter(track_status == "complete") %>% 
  group_by(individual.id, track_id) %>% 
  slice(1) %>% 
  dplyr::select(individual.id, track_id) %>% 
  arrange(individual.id)

# saveRDS(meta, file = "C:/Users/hbronnvik/Documents/storkSSFs/rs_ids_40km_2023-02-08.rds")
