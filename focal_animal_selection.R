### Bird selection
### 2023-01-20
### Hester Bronnvik

### The goal here is to select the birds to be used for route selection analyses on the basis of 
### having completed migrations.
library(move)
library(lubridate)
library(tidyverse)

s <- 21 # the number of days a bird has to have survived since its last migratory day to have completed migration
studies <- c(24442409, 212096177, 76367850, 21231406, 1176017658, 173641633)#, 908232414)

# clean_files <- list.files("C:/Users/hbronnvik/Documents/storkSSFs/clean_data/", pattern = "embc", full.names = T)
visDODs <- readRDS("C:/Users/hbronnvik/Documents/storkSSFs/visual_ACC_deaths_2023-02-01.rds")
load("C:/Users/hbronnvik/Documents/loginStored.rdata")

# migration_locations <- lapply(clean_files, function(x){
#   locs <- readRDS(x)
#   locs$acceleration.raw.x <- NULL
#   locs$acceleration.raw.y <- NULL
#   locs$acceleration.raw.z <- NULL
#   locs$barometric.height <- NULL
#   locs$battery.charge.percent <- NULL
#   locs$battery.charging.current <- NULL
#   locs$external.temperature <- NULL
#   locs$gps.hdop <- NULL
#   locs$gps.time.to.fix <- NULL
#   locs$height.above.msl <- NULL
#   locs$light.level <- NULL
#   locs$manually.marked.outlier <- NULL
#   locs$ornitela.transmission.protocol <- NULL
#   locs$tag.voltage <- NULL
#   locs$manually.marked.valid <- NULL
#   return(locs)
# }) %>% reduce(rbind)
# 
# # determine the identities of the nestling birds (remove any care center adults)
# info <- lapply(studies, function(x){
#   print(x)
#   birds <- getMovebankReferenceTable(study = x, login = loginStored) %>%
#     drop_na(animal_id) %>%
#     filter(sensor_type_id == 653)
#   if("animal_life_stage" %in% colnames(birds)){
#     chicks <- birds %>% 
#       filter(grepl("juv|chick|nestling", animal_life_stage, ignore.case = T) & grepl("release", animal_comments, ignore.case = T) == F) %>% 
#       select(animal_id, study_id)
#     juv <- birds %>% 
#       filter(animal_life_stage == "" & grepl("release|adult", animal_comments, ignore.case = T) == F) %>% 
#       select(animal_id, study_id)
#     chicks <- rbind(chicks, juv)
#   }else{
#     chicks <- birds %>% 
#       filter(!grepl("release|adult", birds$animal_comments, ignore.case = T)) %>% 
#       select(animal_id, study_id)
#   }
#   return(chicks)
# }) %>% reduce(rbind)
# 
# # remove the information from non-migration, this should remove all birds that did not migrate
# migration_locations <- migration_locations %>% 
#   filter(individual.id %in% info$animal_id & !str_detect(track_id, "no_migration"))

# take the threshold number of days before death
# are they migratory?
# filter and save the final IDs

migration_locations <- readRDS("C:/Users/hbronnvik/Documents/storkSSFs/clean_data/full_migration_locations_embc3_15min_2023-01-31.rds")

# define birds that took eastern routes as ones that are ever east of 15.7 longitude (East Germany)
eastern_birds <- unique(migration_locations$individual.id[migration_locations$location.long > 15.7])
# remove the eastern birds
migration_locations <- migration_locations %>% 
  filter(!individual.id %in% eastern_birds)
# get the track IDs from the migratory birds
rs_ids <- readRDS("C:/Users/hbronnvik/Documents/storkSSFs/rs_ids_embc_2023-02-01.rds") %>% 
  filter(individual.id %in% unique(migration_locations$individual.id))
# add the number of journeys
locs <- migration_locations %>% 
  rowwise() %>% 
  mutate(journey_number = rs_ids$journey_number[which(rs_ids$track_id == track_id)]) %>% 
  ungroup()
# the classification of tracks split the data based on year,
# but storks do not care about years, some migrate in December, others in July
# here, we correct the track_id column to account for this.
# First, we re-run EMbC (which is very time-consuming), then we group journeys based
# on time lag and allow the most common track ID to set the others
test <- split(locs, locs$individual.id)[1]
corrected_locs <- lapply(split(locs, locs$individual.id), function(x){
  # prep for EMbC
  ind <- x %>% 
    select(timestamp, location.long, location.lat) %>% 
    as.data.frame()
  # run EMbC
  print(paste0("Running EMbC out for individual ", unique(x$individual.id), ". "), quote = F)
  bc <- EMbC::stbc(ind, smth = 24)
  
  ind <- x %>% 
    mutate(embc_category = bc@A) %>%
    filter(embc_category == 3) %>%
    as.data.frame()
  
  season <- lapply(1:nrow(ind), function(y){
    s <- str_split(ind$track_id, "_")[[y]][3]
  }) %>% unlist()
  
  # find days with six weeks' non-migration on either side
  solo <- ind %>% 
    group_by(date(timestamp)) %>% 
    slice(1) %>% 
    ungroup() %>% 
    mutate(gap_lag = difftime(timestamp, lag(timestamp), units = "weeks"),
           gap_lead = difftime(timestamp, lead(timestamp), units = "weeks"),
           solo = ifelse(abs(gap_lag) >= 6 & abs(gap_lead) >= 6, T, F)) %>% 
    filter(solo == T) %>% 
    select(timestamp) %>% 
    deframe()
  print(paste0("Correcting sticky ends for individual ", unique(x$individual.id), ". "), quote = F)
  # account for those long-distance movements that are not migration (solo)
  if(length(solo) > 0){
    ind <- ind %>%
      mutate(season = season) %>% 
      filter(!date(timestamp) %in% date(solo)) %>% 
      mutate(gap_lag = difftime(timestamp, lag(timestamp), units = "weeks"),
             gap_lead = difftime(timestamp, lead(timestamp), units = "weeks"),
             gap_lag = ifelse(is.na(gap_lag), 0, gap_lag),
             new_season = ifelse(gap_lag >= 6, T, F),
             check = cumsum(new_season)) %>% 
      group_by(check) %>% 
      mutate(track_id = names(which(table(track_id) == max(table(track_id)))))
  }else{
    ind <- ind %>% 
      mutate(season = season,
             gap_lag = difftime(timestamp, lag(timestamp), units = "weeks"),
             gap_lead = difftime(timestamp, lead(timestamp), units = "weeks"),
             gap_lag = ifelse(is.na(gap_lag), 0, gap_lag),
             new_season = ifelse(gap_lag >= 6, T, F),
             check = cumsum(new_season)) %>% 
      group_by(check) %>% 
      mutate(track_id = names(which(table(track_id) == max(table(track_id)))))
  }
  return(ind)
}) %>% reduce(rbind)
# saveRDS(corrected_locs, file = "C:/Users/hbronnvik/Documents/storkSSFs/corrected_locs_2023-02-03.rds")

build <- ind %>%
  mutate(season = season,
         date = date(timestamp)) %>% 
  filter(!date(timestamp) %in% date(solo)) %>% 
  group_by(date) %>% 
  slice(1) %>% 
  ungroup() %>% 
  mutate(gap_lag = difftime(timestamp, lag(timestamp), units = "weeks"),
         gap_lead = difftime(lead(timestamp), timestamp, units = "weeks"),
         gap_lag = ifelse(is.na(gap_lag), 0, gap_lag),
         gap_lead = ifelse(is.na(gap_lead), 0, gap_lead),
         new_season = ifelse(gap_lag >= 5, T, F),
         check = cumsum(new_season)) %>% 
  group_by(track_id, check) %>% 
  mutate(unit_lag = head(gap_lag, 1),
         unit_lead = tail(gap_lead, 1),
         unit_time = n(),
         unit_travel = abs(tail(location.lat, 1) - head(location.lat,1))) %>% 
  ungroup() %>% 
  mutate(bad = ifelse(unit_lag > 3 & unit_lead > 3 & unit_travel < 1 & unit_time <= 2, T, F)) %>% 
  filter(bad == F)

ggplot(build, aes(location.long, location.lat, color = track_id)) +
  geom_path() +
  geom_point(data = build, aes(color = as.character(month(timestamp)))) +
  scale_color_viridis(discrete = T) +
  theme_classic()

# check <- migration_locations %>% 
#   filter(individual.id == "1176031140") %>% 
#   mutate(hour = round_date(timestamp, "hour")) %>% 
#   group_by(hour) %>% 
#   slice(1) %>% 
#   ungroup() %>% 
#   select(timestamp, location.long, location.lat) %>% 
#   # mutate(cat = bc@A) %>% 
#   as.data.frame()
# 
# bc <- EMbC::stbc(check, smth = 24)
# 
# check <- migration_locations %>% 
#   filter(individual.id == "1176031140") %>% 
#   mutate(hour = round_date(timestamp, "hour")) %>% 
#   group_by(hour) %>% 
#   slice(1) %>% 
#   ungroup() %>% 
#   mutate(cat = bc@A) %>%
#   filter(cat == 3) %>%
#   as.data.frame()
# 
# season <- lapply(1:nrow(check), function(x){
#   s <- str_split(check$track_id, "_")[[x]][3]
# }) %>% unlist()
# 
# # find days with six weeks' non-migration on either side
# solo <- check %>% 
#   group_by(date(timestamp)) %>% 
#   slice(1) %>% 
#   ungroup() %>% 
#   mutate(gap_lag = difftime(timestamp, lag(timestamp), units = "weeks"),
#          gap_lead = difftime(timestamp, lead(timestamp), units = "weeks"),
#          solo = ifelse(abs(gap_lag) >= 6 & abs(gap_lead) >= 6, T, F)) %>% 
#   filter(solo == T) %>% 
#   select(timestamp) %>% 
#   deframe()
# 
# build <- check %>% 
#   # filter(!date(timestamp) %in% date(solo)) %>% 
#   mutate(season = season,
#          gap_lag = difftime(timestamp, lag(timestamp), units = "weeks"),
#          gap_lead = difftime(timestamp, lead(timestamp), units = "weeks"),
#          gap_lag = ifelse(is.na(gap_lag), 0, gap_lag),
#          new_season = ifelse(gap_lag >= 6, T, F),
#          check = cumsum(new_season)) %>% 
#   group_by(check) %>% 
#   mutate(tid = names(which(table(track_id) == max(table(track_id)))))
# find any animals that slipped through the cracks for estimation of DOD
lost <- migration_locations %>% 
  group_by(individual.id) %>% 
  slice(n()) %>% 
  ungroup() %>% 
  dplyr::select(individual.id, individual.local.identifier) %>% 
  filter(!individual.id %in% visDODs$individual_id)

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
    filter(track_id == track_id[n()])
  # take either the confirmed DOD or the last transmitted time stamp
  loss_time <- visDODs %>% 
    filter(individual_id == unique(x$individual.id)) %>% 
    mutate(dod = as.POSIXct(ifelse(is.na(dod), info$timestamp_end[info$individual_id == individual_id], dod), tz = "UTC", origin = "1970-01-01")) %>% 
    select(dod) %>% 
    deframe()
  # compare the DOD to the end of the migration
  loss <- max(final_track$timestamp) > loss_time - days(s)
  # add a column containing the outcome of the migratory track
  x <- x %>% 
    group_by(track_id) %>% 
    mutate(track_status = ifelse(loss == T & track_id == unique(final_track$track_id), "incomplete", "complete")) %>% 
    ungroup()
  return(x)
}) %>% reduce(rbind)


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
  group_by(individual.id) %>% 
  slice(1) %>% 
  dplyr::select(individual.id, track_id) %>% 
  arrange(individual.id)

# saveRDS(meta, file = "C:/Users/hbronnvik/Documents/storkSSFs/rs_ids_embc_2023-02-01.rds")
