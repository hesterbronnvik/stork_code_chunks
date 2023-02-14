### data location to find the names of birds that should be included (ie. have completed more than one migration)
### Hester Brønnvik
### 29.09.2022

library(lubridate)
library(geosphere)
library(EMbC)
library(move)
library(stringr)
library(tidyverse)
d_thresh <- 40000 # meters
w_thresh <- 6 # weeks
s_thresh <- 40 # m/s
l_thresh <- 3 # degrees latitude

# required information
load("C:/Users/hbronnvik/Documents/loginStored.rdata")
studies <- c(24442409, 212096177, 76367850, 21231406, 1176017658, 173641633)


# determine the identities of the nestling birds (remove any care center adults)
info <- lapply(studies, function(x){
  print(x)
  birds <- getMovebankReferenceTable(study = x, login = loginStored) %>%
    drop_na(animal_id) %>%
    filter(sensor_type_id == 653)
  if("animal_life_stage" %in% colnames(birds)){
    chicks <- birds %>% 
      filter(grepl("juv|chick|nestling", animal_life_stage, ignore.case = T) & grepl("release", animal_comments, ignore.case = T) == F) %>% 
      dplyr::select(animal_id, study_id, animal_local_identifier)
    juv <- birds %>% 
      filter(animal_life_stage == "" & grepl("release|adult", animal_comments, ignore.case = T) == F) %>% 
      dplyr::select(animal_id, study_id, animal_local_identifier)
    chicks <- rbind(chicks, juv)
  }else{
    chicks <- birds %>% 
      filter(!grepl("release|adult", birds$animal_comments, ignore.case = T)) %>% 
      dplyr::select(animal_id, study_id, animal_local_identifier)
  }
  return(chicks)
}) %>% reduce(rbind)

# select the nestling data files
full_files <- list.files("C:/Users/hbronnvik/Documents/storkSSFs/full_data", full.names = T)
todo <- str_sub(full_files, 50, -16)
todo <- todo[todo %in% info$animal_id]

full_files <- lapply(1:length(full_files), function(x){
  name <- str_sub(full_files[x], 50, -16)
  if(name %in% todo){
    return(full_files[x])
  }
}) %>% unlist()

# remove errors
start_time <- Sys.time()
clean_locations <- lapply(full_files, function(x){
  # load the data
  ind <- readRDS(x) # "C:/Users/hbronnvik/Documents/storkSSFs/full_data/78031713_2023-01-28.rds"
  # clean the data
  locs_df <- ind %>% 
    drop_na(location.long) %>% 
    mutate(index = row_number())
  # remove duplicated locations because they prevent accurate calculations of distance and speed
  doubles <- locs_df[which(locs_df$timestamp %in% locs_df[duplicated(locs_df$timestamp), "timestamp"]),] %>% 
    filter(is.na(height.above.ellipsoid))
  
  locs_df <- locs_df %>% 
    filter(!index %in% doubles$index) 
  
  # warn if a duplicated timestamp contains information other than location (not usually the case)
  if(nrow(locs_df[which(locs_df$timestamp %in% locs_df[duplicated(locs_df$timestamp),
                                                       "timestamp"]),]) > 0){print("Duplicates containing HAE and DOP values exist.", quote = F)}
  
  # remake doubles in the event of duplicates that hold values
  doubles <- locs_df[which(locs_df$timestamp %in% locs_df[duplicated(locs_df$timestamp), "timestamp"]),]%>% 
    mutate(event = round_date(timestamp, "minute"))
  
  if(nrow(doubles) > 0){
    # if there is more than one instance of duplicates with information
    check <- lapply(unique(round_date(doubles$timestamp, "minute")), function(q){
      # take the duplicates within one minute
      dd <- doubles %>% 
        filter(event == q)
      # determine the last location before the duplicates (point of interest)
      poi <- locs_df %>% 
        filter(timestamp < dd$timestamp[1]) %>% 
        slice(n())
      # find the distance from each duplicate to the poi, even for true points this may increase
      cc <- lapply(1:nrow(dd), function(p){
        d <- distVincentyEllipsoid(c(poi$location.long, poi$location.lat), c(dd$location.long[p], dd$location.lat[p]))
        return(d)
      }) %>% unlist()
      # calculate the ground speeds from each duplicate to the poi
      dd <- dd %>% 
        mutate(dist_from_unique = cc, 
               time_since_unique = as.numeric(difftime(dd$timestamp[1], poi$timestamp, units = "sec")),
               speed_after_unique = dist_from_unique/time_since_unique) %>% 
        filter(speed_after_unique > s_thresh)
      return(dd)
    }) %>% reduce(rbind)
    # filter out the ground speeds higher than reasonable
    locs_df <- locs_df %>% 
      filter(!index %in% check$index)
  }
  
  # down-sample to 15 minutes
  locs_df <- locs_df %>%
    mutate(td = as.numeric(difftime(timestamp, lag(timestamp), units = "secs"))) %>%
    filter(td >= 300) %>%
    mutate(seq15 = round_date(timestamp, unit = "15 minutes")) %>%
    group_by(seq15) %>%
    slice(1) %>%
    ungroup() %>%
    dplyr::select(-seq15, -td) %>% 
    # calculate ground speeds
    mutate(distance = distVincentyEllipsoid(cbind(location.long, location.lat), cbind(lag(location.long), lag(location.lat))),
           timediff = as.numeric(difftime(timestamp, lag(timestamp), units = "secs")),
           ground_speed_15 = distance/timediff) %>% 
    # remove absurd speeds to clean outliers
    filter(ground_speed_15 < 50)
  
  # add daily metrics
  ind <- locs_df %>% 
    filter(!is.na(location.lat)) %>% 
    mutate(date = date(timestamp)) %>%
    group_by(date) %>% 
    # the Haversine distance between first and last locations of the day
    mutate(daily_dist = distVincentyEllipsoid(c(head(location.long,1), head(location.lat, 1)), c(tail(location.long,1), tail(location.lat, 1))),
           # the rhumbline bearing between first and last locations of the day
           daily_direction = bearingRhumb(c(head(location.long,1), head(location.lat, 1)), c(tail(location.long,1), tail(location.lat, 1))),
           compass_direction = ifelse(daily_direction > 90 & daily_direction < 270, "southward", "northward")) %>% 
    ungroup() %>% 
    arrange(timestamp)

  # ggplot(ind, aes(location.long, location.lat, color = phase)) +
  #   borders(xlim = c(-5, 5), ylim = c(30, 48)) +
  #   geom_point() +
  #   theme_classic()
  print(ind$individual.local.identifier[1])
  return(ind)
})
# DOP duplucates:
# "Snöfrid + A6Y49 (eobs 7971)"
# "Wiesenmann (A6R70, eobs 8027)"
# "Maximilian (E0343, e-obs 8032)"

# saveRDS(clean_locations, file = "C:/Users/hbronnvik/Documents/storkSSFs/clean_locations_2023-02-10.rds")

clean_locations <- clean_locations[lapply(clean_locations, nrow) > 1]

# take each phase and find long-distance travel days
migration_locations <- lapply(clean_locations, function(y){
  m <- lapply(split(y, y$phase), function(x){
    print(paste0(unique(x$individual.local.identifier), " ", unique(x$phase)))
    
    # allow rbinding in spite of column differences
    if(!"manually.marked.valid" %in% colnames(x)){
      x <- x %>% 
        mutate(manually.marked.valid = NA)
    }
    if(!"acceleration.raw.x" %in% colnames(x)){
      x <- x %>% 
        mutate(acceleration.raw.x = NA,            
               acceleration.raw.y = NA,           
               acceleration.raw.z = NA,          
               barometric.height = NA,            
               battery.charge.percent = NA,         
               battery.charging.current = NA,      
               external.temperature = NA,          
               gps.hdop = NA,                       
               gps.time.to.fix = NA,              
               height.above.msl = NA,             
               light.level = NA,                   
               ornitela.transmission.protocol = NA,
               tag.voltage = NA)
    }
    # Each bout consists only of > 40km travel days
    # Each bout consists of more than one day
    # Each bout consists of locations not separated by 6 weeks or more
    
    # 1 bout:
    # the bout
    
    # 2 bouts:
    # the bout with more displacement
    
    # 3+ bouts:
    # The first bout with unit lag and unit lead greater than the threshold
    # or (if there are stopovers), the first with lead < threshold and lag > threshold
    df <- x %>% 
      mutate(seq_15 = round_date(timestamp, unit = "15 minutes"),
             displacement = max(location.lat) - min(location.lat)) %>% 
      group_by(seq_15) %>% 
      slice(1) %>% 
      ungroup() %>% 
      filter(daily_dist > d_thresh & displacement > l_thresh) %>% 
      # is there a large gap between points?
      mutate(gap_lag = difftime(timestamp, lag(timestamp), units = "weeks"),
             gap_lead = difftime(lead(timestamp), timestamp, units = "weeks"),
             gap_lag = ifelse(is.na(gap_lag), 0, gap_lag),
             gap_lead = ifelse(is.na(gap_lead), 0, gap_lead),
             new_season = ifelse(gap_lag >= w_thresh, T, F),
             check = cumsum(new_season)) %>% 
      dplyr::select(timestamp, gap_lag, gap_lead, new_season, check, daily_dist, daily_direction, date, location.lat, location.long) %>% 
      group_by(check) %>% 
      mutate(unit_lag = head(gap_lag, 1),
             unit_lead = tail(gap_lead, 1),
             days = length(unique(date(timestamp)))) %>% 
      ungroup()
    # if there is migration:
    if(length(unique(df$date)) > 1){
      # if there are multiple bouts,
      if(length(unique(df$check)) > 2){
        # if those bouts are separated by stopovers,
        if(sum(unique(df$unit_lag) < w_thresh) > 1){
          # the first time a bout of migration is within six weeks of another
          on <- df$timestamp[which(df$unit_lead < w_thresh)[1]]
          # the last time a bout is within six weeks of another
          off <- df$timestamp[which(df$unit_lag > w_thresh)[length(which(df$unit_lag > w_thresh))]]
        }else{# if those bouts are separated by breeding/wintering
          # the bout with both unit lag and unit lead greater than w_thresh
          on <- df$timestamp[which(df$unit_lead > w_thresh & df$unit_lag > w_thresh)[1]]
          off <- df$timestamp[which(df$unit_lead > w_thresh & df$unit_lag > w_thresh)[nrow(df[which(df$unit_lead > w_thresh & df$unit_lag > w_thresh),])]]
        }
      }else{
        # if there are only one or two bouts, use the one with the longest displacement
        d0 <- df %>% 
          filter(check == unique(check)[1]) %>% 
          summarize(d = max(location.lat) - min(location.lat)) %>% 
          deframe()
        d1 <- df %>% 
          filter(check == unique(check)[length(unique(check))]) %>% 
          summarize(d = max(location.lat) - min(location.lat)) %>% 
          deframe()
        if(d1 > d0){
          on <- df %>% 
            filter(check == unique(check)[length(unique(check))]) %>% 
            slice(1) %>% 
            dplyr::select(timestamp) %>% 
            deframe()
          off <- off <- df %>% 
            filter(check == unique(check)[length(unique(check))]) %>% 
            slice(n()) %>% 
            dplyr::select(timestamp) %>% 
            deframe()
        }else{
          on <- df$timestamp[which(df$check == unique(df$check)[1])][1]
          off <- df$timestamp[which(df$check == unique(df$check)[1])][nrow(df[which(df$check == unique(df$check)[1]),])]
        }
      }
      # ggplot(df %>% group_by(date) %>% slice(1) %>% ungroup(), aes(x = timestamp, y = daily_dist, color = check)) +
      #   geom_segment(aes(x = timestamp, xend = timestamp, y = 0, yend = daily_dist), color="grey") +
      #   geom_point() +
      #   theme_classic() +
      #   theme(
      #     panel.grid.major.x = element_blank(),
      #     panel.border = element_blank(),
      #     axis.ticks.x = element_blank()
      #   ) +
      #   xlab("Time") +
      #   ylab("Daily distance (m)")
      
      migration <- x %>% 
        filter(between(timestamp, on, off)) %>% 
        mutate(track_id = paste(unique(x$individual.id), unique(x$phase), sep = "_")) 
    }else{
      migration <- x %>% 
        slice(1) %>% 
        mutate(track_id = NA)
    }
    migration <- migration %>% 
      mutate(td = as.numeric(difftime(timestamp, lag(timestamp), units = "secs"))) %>% 
      filter(td >= 300) %>%
      dplyr::select(individual.id, individual.local.identifier, timestamp, location.lat, location.long, phase, track_id, daily_dist, daily_direction, ground.speed) %>% 
      mutate(seq15 = round_date(timestamp, unit = "15 minutes")) %>% 
      group_by(seq15) %>% 
      slice(1) %>% 
      ungroup() %>% 
      dplyr::select(-seq15)
    saveRDS(migration, file = paste0("C:/Users/hbronnvik/Documents/storkSSFs/classified_gps/migration_locations_40km15min_", paste(unique(x$individual.id), unique(x$phase), sep = "_"), "_", Sys.Date(), ".rds"))
    return(migration)
  })# %>% reduce(rbind)
  return(m)
})# %>% reduce(rbind) # using pre-cleaned data, 8.95 min
Sys.time() - start_time
saveRDS(migration_locations, file = "C:/Users/hbronnvik/Documents/storkSSFs/migration_locations_40km15min_20230210.rds")

ggplot(migration_locations[[81]][[1]] %>% drop_na(track_id), aes(location.long, location.lat, color = track_id)) +
  borders(xlim = c(-5, 5), ylim = c(30, 50)) +
  geom_point() +
  theme_classic() +
  facet_wrap(~track_id)

x <- readRDS("C:/Users/hbronnvik/Documents/storkSSFs/classified_gps/migration_locations_40km15min_1578616761_pre-breeding_2021_2023-02-07.rds")






ind_full <- readRDS("C:/Users/hbronnvik/Documents/storkSSFs/full_data/173199826_2023-01-28.rds")

ind_fa_17 <- ind %>% 
  filter(phase == "post-breeding_2017") %>% 
  mutate(td = as.numeric(difftime(timestamp, lag(timestamp), units = "secs"))) %>% 
  filter(td >= 300)%>%
  mutate(seq15 = round_date(timestamp, unit = "15 minutes")) %>% 
  group_by(seq15) %>% 
  slice(1) %>% 
  ungroup() %>% 
  dplyr::select(-seq15) %>% 
  group_by(date(timestamp)) %>% 
  slice(1) %>% 
  ungroup()

ggplot(ind_fa_17, aes(location.long, location.lat, color = daily_dist > 40000)) +
  borders(xlim = c(-10, 5), ylim = c(30, 50)) +
  geom_point() +
  theme_classic()



ggplot(complete_ml %>% filter(individual.id %in% unique(individual.id)[1:(length(unique(individual.id))/2)]), aes(location.long, location.lat, color = individual.id)) +
  borders(xlim = c(-5, 5), ylim = c(30, 50)) +
  geom_point() +
  theme_classic() +
  facet_wrap(~track_id)


n <- lapply(clean_locations, function(x){
  df <- data.frame(unique(x$individual.id))
  return(df)
}) %>% reduce(rbind)

y <- clean_locations[[190]]
x <- split(y, y$phase)[[1]]

clean_locations <- readRDS("C:/Users/hbronnvik/Documents/storkSSFs/clean_locations_2023-02-10.rds")

all_long_dist <- lapply(clean_locations, function(y){
  m <- lapply(split(y, y$phase), function(x){
    print(paste0(unique(x$individual.local.identifier), " ", unique(x$phase)))
    
    # allow rbinding in spite of column differences
    if(!"manually.marked.valid" %in% colnames(x)){
      x <- x %>% 
        mutate(manually.marked.valid = NA)
    }
    if(!"acceleration.raw.x" %in% colnames(x)){
      x <- x %>% 
        mutate(acceleration.raw.x = NA,            
               acceleration.raw.y = NA,           
               acceleration.raw.z = NA,          
               barometric.height = NA,            
               battery.charge.percent = NA,         
               battery.charging.current = NA,      
               external.temperature = NA,          
               gps.hdop = NA,                       
               gps.time.to.fix = NA,              
               height.above.msl = NA,             
               light.level = NA,                   
               ornitela.transmission.protocol = NA,
               tag.voltage = NA)
    }
    
    df <- x %>% 
      mutate(td = as.numeric(difftime(timestamp, lag(timestamp), units = "secs"))) %>% 
      filter(td >= 300) %>% 
      mutate(seq_15 = round_date(timestamp, unit = "15 minutes"),
             displacement = max(location.lat) - min(location.lat)) %>% 
      group_by(seq_15) %>% 
      slice(1) %>% 
      ungroup() %>% 
      filter(daily_dist > d_thresh & displacement > l_thresh)
    if(length(unique(df$date)) > 1){
      on <- df$timestamp[1]
      off <- df$timestamp[nrow(df)]
      
      migration <- x %>% 
        filter(between(timestamp, on, off))%>% 
        mutate(td = as.numeric(difftime(timestamp, lag(timestamp), units = "secs"))) %>% 
        filter(td >= 300) %>% 
        mutate(track_id = paste(unique(x$individual.id), unique(x$phase), sep = "_"))  %>%
        dplyr::select(individual.id, individual.local.identifier, timestamp, location.lat, location.long, phase, track_id, daily_dist, daily_direction, ground.speed) %>% 
        mutate(seq15 = round_date(timestamp, unit = "15 minutes")) %>% 
        group_by(seq15) %>% 
        slice(1) %>% 
        ungroup() %>% 
        dplyr::select(-seq15)
    }else{
      migration <- x %>% 
        slice(1) %>% 
        mutate(track_id = NA) %>% 
        dplyr::select(individual.id, individual.local.identifier, timestamp, location.lat, location.long, phase, track_id, daily_dist, daily_direction, ground.speed)}

    
  }) %>% reduce(rbind)
  return(m)
}) %>% reduce(rbind) %>% 
  drop_na(track_id)

ggplot(migration %>% filter(daily_dist > 40000), aes(location.long, location.lat)) +
  borders(xlim = c(-10, 5), ylim = c(30, 50)) +
  geom_point(color = "black") +
  theme_classic()

ggplot(all_long_dist %>% filter(daily_dist > 40000), aes(location.long, location.lat)) +
  borders(xlim = c(-10, 5), ylim = c(30, 50)) +
  geom_point(color = "black") +
  theme_classic() +
  facet_wrap(~track_id)

ld_dates_end <- all_long_dist %>% 
  group_by(track_id, date(timestamp)) %>% 
  slice(1) %>% 
  ungroup() %>% 
  group_by(track_id) %>% 
  mutate(n = n(),
         id = cur_group_id()) %>% 
  filter(n > 7) %>% 
  slice((n()-7):n()) %>% 
  mutate(timed = as.numeric(difftime(timestamp, lag(timestamp), units = "weeks")))

ld_dates_start <- all_long_dist %>% 
  group_by(track_id, date(timestamp)) %>% 
  slice(1) %>% 
  ungroup() %>% 
  group_by(track_id) %>% 
  mutate(n = n(),
         id = cur_group_id()) %>% 
  filter(n > 7) %>% 
  slice(1:7) %>% 
  mutate(timed = as.numeric(difftime(timestamp, lag(timestamp), units = "weeks")))

ld_dates <- rbind(ld_dates_end, ld_dates_start)

# broken tracks
check <- all_long_dist %>% 
  filter(track_id %in% unique(ld_dates$track_id[which(ld_dates$timed > 0.5)]))

eastern_birds <- unique(check$individual.id[check$location.long > 15.7])
# remove the eastern birds
check <- check %>% 
  filter(!individual.id %in% eastern_birds)
ggplot(check %>% filter(daily_dist > 40000), aes(location.long, location.lat)) +
  borders(xlim = c(-10, 5), ylim = c(30, 50)) +
  geom_point() +
  theme_classic() +
  facet_wrap(~track_id)

ggplot(track %>% filter(year(timestamp) == 2021 & month(timestamp) %in% c(1:7)), aes(location.long, location.lat, color = month(timestamp))) +
  borders(xlim = c(-10, 5), ylim = c(30, 50)) +
  geom_point() +
  theme_classic()

unique(check$track_id)
# "1178290125_pre-breeding_2021" # has a large transmission gap, but appears to be valid
# "1578604607_post-breeding_2021" # same as above, same area 
# "80436240_post-breeding_2015"  # same as above, same area


# saveRDS(all_long_dist, file = "C:/Users/hbronnvik/Documents/storkSSFs/classified_gps/full_all_long_dist_40km15min_2023-02-09.rds")

# the data down-sampled to 15 minutes, with migration defined as moving more than 40 km/day, 
# pre-breeding as Dec. to the next July, and post-breeding as Aug. to Nov.
# all_long_dist <- readRDS("C:/Users/hbronnvik/Documents/storkSSFs/classified_gps/full_all_long_dist_40km15min_2023-02-09.rds")

# define birds that took eastern routes as ones that are ever east of 15.7 longitude (East Germany)
eastern_birds <- unique(all_long_dist$individual.id[all_long_dist$location.long > 15.7])
# remove the eastern birds
all_long_dist <- all_long_dist %>% 
  filter(!individual.id %in% eastern_birds)
# get the track IDs from the migratory birds
rs_ids <- all_long_dist %>% 
  group_by(individual.id, track_id) %>% 
  slice(1) %>% 
  dplyr::select(individual.id, track_id) %>% 
  ungroup()
# the total number of migrations attempted and completed by each animal in each season 
meta <- all_long_dist %>%
  mutate(season = ifelse(grepl("post", phase), "post-breeding", "pre-breeding")) %>% 
  group_by(individual.id, season) %>% 
  count(track_id) %>% 
  mutate(journey_number = row_number()) %>% 
  ungroup() %>% 
  group_by(individual.id) %>%
  mutate(total_journeys = n()) %>% 
  ungroup()
# add the number of journeys
all_long_dist <- all_long_dist %>% 
  rowwise() %>% 
  mutate(journey_number = meta$journey_number[which(meta$track_id == track_id)]) %>% 
  ungroup()


# take the last time stamps from the full GPS data for each animal
info <- lapply(studies, function(x){
  i <- getMovebankAnimals(study = x , login = loginStored) %>%
    filter(sensor_type_id == 653 & individual_id %in% unique(migration_locations$individual.id)) %>% 
    mutate(study = x)
  return(i)
}) %>% reduce(rbind) %>% 
  mutate(timestamp_end = as.POSIXct(timestamp_end, tz = "UTC"))

# find the success/failure of each migration track
all_long_dist <- lapply(split(all_long_dist, all_long_dist$individual.id), function(x){
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

complete_ml <- all_long_dist %>% 
  filter(track_status == "complete")

ggplot(complete_ml %>% filter(daily_dist > 40000 & grepl("post-", track_id)), aes(location.long, location.lat)) +
  borders(xlim = c(-10, 5), ylim = c(30, 50)) +
  geom_point(color = "black") +
  theme_classic() +
  facet_wrap(~track_id)







## take the clean_locations at 15 min intervals and with ground_speed below 50 m/s

# did the animal migrate?
migratory <- lapply(clean_locations, function(x){
  m <- ifelse(max(x$daily_dist) > d_thresh, T, F)
  w <- ifelse(max(x$location.long)  > 15.7, F, T)
  mt <- data.frame(individual.id = unique(x$individual.id),
                   migrated = m,
                   western = w)
  return(mt)
}) %>% reduce(rbind) %>% 
  filter(migrated == T & western == T)

# are there long-distance movements in nov/jan?
ld <- lapply(clean_locations, function(x){
  t <- x %>% 
    filter(month(timestamp) %in% c(11, 12, 1))
  l <- ifelse(max(t$daily_dist) > d_thresh, T, F)
  mt <- data.frame(individual.id = unique(x$individual.id),
                   migrated = l)
  return(mt)
}) %>% reduce(rbind) %>% 
  filter(migrated == T)


sub <- rbind(clean_locations[[1]],
             clean_locations[[2]],
             clean_locations[[3]],
             clean_locations[[4]],
             clean_locations[[5]]) %>% 
  as.data.frame()

library(EMbC)
cl <- lapply(clean_locations, function(x){
  df <- x %>% 
    dplyr::select(timestamp, location.long, location.lat, individual.id, individual.local.identifier, ground_speed_15, daily_dist)
}) %>% reduce(rbind) %>% 
  filter(individual.id %in% migratory$individual.id) %>% 
  as.data.frame()

bc <- stbc(cl %>% filter(individual.id %in% unique(individual.id)[7:15]))
sctr(bc)
view(bc)







## long distance re-locations are short in time
## they are in off seasons (11-2)
## they may be lateral (225-315, 45-135)

# find long-distance movements that may not be migration and remove them if need be
ld2 <- lapply(clean_locations, function(x){
  # x <- clean_locations[[100]]
  ld <- x %>% 
    filter(daily_dist > d_thresh) %>% 
    group_by(date(timestamp)) %>% 
    slice(1) %>% 
    ungroup() %>% 
    mutate(gap = as.numeric(difftime(timestamp, lag(timestamp), units = "days")),
           lateral = ifelse(between(daily_direction, 225, 315) | between(daily_direction, 45, 135), T, F),
           potential_ldd = ifelse(gap > 7 & lateral == T & month(timestamp) %in% c(1, 2, 11, 12), T, F)) %>% 
    ungroup() %>% 
    mutate(pldd = ifelse(length(unique(ld$potential_ldd)) == 2, T, F)) %>% 
    dplyr::select(timestamp, location.long, location.lat, individual.id, individual.local.identifier, ground_speed_15, daily_dist, pldd, gap, potential_ldd)
  
  
  # png(paste0("C:/Users/hbronnvik/Documents/storkSSFs/m_locs_check/", unique(x$individual.id), ".png"))
  print(ggplot(ld, aes(location.long, location.lat, color = timestamp)) +
          borders(xlim = c(-10,10), ylim = c(30, 50)) +
          geom_path() +
          geom_point(data = ld %>% filter(potential_ldd == T), cex = 1, color = "red") +
          labs(title = ld$individual.id[1]) +
          theme_classic())
  # dev.off()
  
  if(nrow(ld) > 0){
    return(ld)
  }else{
    ld <- x %>% 
      slice(1) %>% 
      mutate(gap = NA,
             lateral = NA,
             potential_ldd = NA,
             pldd = NA) %>% 
      dplyr::select(timestamp, location.long, location.lat, individual.id, individual.local.identifier, ground_speed_15, daily_dist, pldd, gap, potential_ldd)
    return(ld)
  }
}) %>% reduce(rbind) %>% 
  filter(pldd == T & individual.id %in% migratory$individual.id)

f <- x %>% 
  mutate(date = date(timestamp),
         week = ceiling(as.numeric(difftime(date, date[1], units = "weeks")))) %>% 
  group_by(week) %>% 
  mutate(weekly_dist = distVincentyEllipsoid(c(head(location.long,1), head(location.lat, 1)), c(tail(location.long,1), tail(location.lat, 1))),
         weekly_direction = bearingRhumb(c(head(location.long,1), head(location.lat, 1)), c(tail(location.long,1), tail(location.lat, 1))),
         lateral = ifelse(between(weekly_direction, 225, 315) | between(weekly_direction, 45, 135), T, F),
         transit = abs(head(location.lat, 1) - tail(location.lat, 1))) %>% 
  ungroup()

ggplot(f %>% filter(daily_dist > 40000), aes(location.long, location.lat, color = transit)) +
  borders(xlim = c(-10,10), ylim = c(30, 50)) +
  geom_point() +
  theme_classic()

# manual corrections to remove programmatically identified long distance re-locations
corrections <- data.frame(date = c("2021-11-08 02:00:23", "2021-11-10 02:00:23", "2022-11-04 02:00:23",
                                   "2020-11-09 02:00:23", "2020-11-10 02:00:23", "2020-11-12 02:00:23",
                                   "2020-11-01 02:00:23", "2020-11-01 02:00:23", "2020-11-04 02:00:23", "2020-11-07 02:00:23", "2020-11-09 02:00:23", "2020-11-17 02:00:23", "2020-11-26 02:00:23", "2020-11-27 02:00:23", "2020-12-18 02:00:23", "2020-12-19 02:00:23", "2020-01-20 02:00:23", "2020-02-23 02:00:23", "2021-03-01 02:00:23", "2021-03-02 02:00:23", "2021-11-18 02:00:23", "2021-11-19 02:00:23", "2022-02-08 02:00:23", "2022-03-03 02:00:23", "2022-03-04 02:00:23",
                                   "2022-10-01 02:00:23","2022-10-02 02:00:23","2022-10-03 02:00:41","2022-10-15 02:00:23","2022-10-16 02:00:23","2022-10-26 02:00:23","2022-10-30 02:00:23","2022-11-04 02:00:23","2022-11-06 02:00:23","2022-11-08 02:00:23","2022-11-09 02:00:23","2022-11-11 02:00:23","2022-11-12 02:00:23",
                                   "2020-09-15 02:00:23", "2020-09-16 02:00:23", "2020-09-20 02:00:23", "2020-09-21 02:00:23",
                                   "2020-09-25 02:00:23", "2020-09-26 02:00:23", "2020-09-28 02:00:23",
                                   "2020-09-29 02:00:23", "2020-10-01 02:00:23", "2020-10-02 02:00:23",
                                   "2020-10-07 02:00:23", "2020-10-13 02:00:23", "2020-10-14 02:00:23",
                                   "2020-10-15 06:00:09", "2020-11-06 02:00:23", "2020-11-07 02:00:23",
                                   "2020-12-28 10:30:23", "2020-12-29 10:30:23", "2020-12-31 11:10:14",
                                   "2021-01-01 02:00:23", "2021-01-02 02:00:23", "2021-01-03 02:00:23",
                                   "2021-01-04 02:00:23", "2021-01-05 02:00:23",
                                   "2020-10-04 02:00:23",
                                   "2020-10-23 10:15:41", "2020-12-04 02:00:23",
                                   "2020-11-02 02:00:23","2020-11-10 02:00:23","2021-02-12 02:00:22","2021-02-13 02:00:23",
                                   "2021-06-11 02:00:23",
                                   "2021-10-10 02:00:41","2021-10-30 02:00:23","2021-11-05 02:00:23","2021-11-06 02:00:23","2021-11-15 02:00:23","2021-12-13 02:00:23",
                                   "2020-10-22 02:00:23", "2020-11-21 02:00:23", "2021-10-31 02:00:23", 
                                   "2020-10-13 02:00:23", "2020-10-14 02:00:23", "2020-10-15 02:00:23",
                                   "2020-10-19 02:00:23", "2020-10-25 02:00:23", "2020-11-11 02:00:23",
                                   "2020-11-12 02:00:23", "2020-11-24 02:00:23", "2020-11-25 02:00:23",
                                   "2020-11-26 02:00:23", "2020-11-27 02:00:23", "2020-12-10 02:00:23",
                                   "2020-12-20 02:00:23", "2020-12-21 02:00:23", "2021-01-19 02:00:23",
                                   "2021-02-07 02:00:41", "2021-02-08 02:00:23", "2021-02-10 02:00:23",
                                   "2021-02-13 02:00:23",
                                   "2020-10-12 02:00:23", "2020-10-13 02:00:23",
                                   "2020-10-18 02:00:23", "2020-12-05 02:00:23", "2020-12-15 02:00:23",
                                   "2020-12-16 02:00:23", "2020-12-17 02:00:23", "2020-12-18 02:00:23",
                                   "2020-12-20 02:00:41", "2020-12-21 02:00:23", "2020-12-27 02:00:23",
                                   "2020-12-30 02:00:23", "2020-12-31 02:00:23", "2021-01-01 02:00:23",
                                   "2021-01-11 02:00:23", "2021-01-14 02:00:23", "2021-01-15 02:00:23",
                                   "2021-01-16 02:00:23", "2021-12-07 02:00:23",
                                   "2022-11-29 02:00:23","2023-01-27 02:00:23",
                                   "2020-11-08 02:00:39", "2020-11-09 02:00:23", "2020-11-11 02:00:23", "2020-11-12 02:00:23", "2020-11-13 02:00:23","2021-01-17 02:00:23", "2021-01-21 02:00:23", "2021-01-22 02:00:23",
                                   "2021-01-23 02:00:23", "2021-01-24 02:00:23", "2021-01-25 02:00:23",
                                   "2021-02-01 02:00:23", "2021-02-02 02:00:23", "2021-02-03 02:00:23",
                                   "2021-02-05 02:00:23", "2021-02-07 02:00:23", "2021-02-10 02:00:23",
                                   "2021-02-13 02:00:23", "2021-02-24 02:00:23", "2021-03-02 02:00:23",
                                   "2021-03-09 02:00:23", "2021-03-10 02:00:23", "2021-03-12 02:00:23",
                                   "2021-03-13 02:00:23",
                                   "2020-09-30 02:00:41", "2020-10-01 02:00:23", "2020-10-02 02:00:23",
                                   "2020-10-03 02:00:23", "2020-10-12 02:00:23", "2020-10-14 02:00:23",
                                   "2020-10-16 02:00:23", "2020-10-21 02:00:23", "2020-10-22 02:00:23",
                                   "2020-11-02 02:00:23", "2020-12-18 02:00:23", "2020-12-19 02:00:23",
                                   "2020-12-20 02:00:23", "2020-12-21 02:00:23", "2021-01-04 02:00:23",
                                   "2021-01-06 02:00:23", "2021-01-11 02:00:23", "2021-01-15 02:00:23",
                                   "2021-01-16 02:00:23", "2021-01-21 02:00:23", "2021-01-22 02:00:23",
                                   "2021-01-23 02:00:23", "2021-01-24 02:00:23", "2021-01-25 02:00:23",
                                   "2021-02-24 02:00:23", "2021-02-25 02:00:23", "2021-03-05 02:00:23",
                                   "2021-03-06 02:00:23",
                                   "2021-11-01 02:00:23 UTC", "2021-11-08 02:00:23 UTC",
                                   "2021-11-19 02:00:23 UTC", "2021-11-20 02:00:23 UTC", "2021-11-21 02:00:23 UTC",
                                   "2021-12-04 02:00:23 UTC", "2021-12-13 02:00:23 UTC", "2021-12-28 02:00:23 UTC",
                                   "2022-01-03 02:00:23 UTC", "2022-01-15 02:00:23 UTC", "2022-01-18 02:00:23 UTC",
                                   "2022-01-19 02:00:23 UTC", "2022-01-22 02:00:23 UTC", "2022-01-23 02:00:23 UTC",
                                   "2022-01-25 02:00:39 UTC", "2022-01-27 02:00:23 UTC", "2022-01-29 02:00:23 UTC",
                                   "2022-01-31 02:00:23 UTC", "2022-02-01 02:00:23 UTC", "2022-02-03 02:00:23 UTC",
                                   "2022-02-04 02:00:23 UTC", "2022-02-06 02:00:23 UTC", "2022-02-07 02:00:23 UTC",
                                   "2022-02-13 02:00:23 UTC", "2022-02-14 02:00:23 UTC", "2022-03-03 02:00:23 UTC",
                                   "2022-03-04 02:00:23 UTC", "2022-03-05 02:00:23 UTC", "2022-05-29 02:00:23 UTC",
                                   "2022-06-01 02:00:23 UTC", "2022-06-03 02:00:23 UTC",
                                   "2022-06-06 02:00:23 UTC", "2022-06-15 02:00:23 UTC", "2022-06-17 02:00:23 UTC",
                                   "2022-06-19 02:00:23 UTC", "2022-07-03 02:00:23 UTC", "2022-07-04 02:00:23 UTC",
                                   "2022-07-05 02:00:23 UTC", "2022-07-11 02:00:23 UTC", "2022-07-13 02:00:23 UTC",
                                   "2022-07-16 02:00:23 UTC", "2022-07-20 02:00:23 UTC", "2022-07-21 02:00:23 UTC",
                                   "2022-08-07 02:00:23 UTC", "2023-01-21 02:00:23", "2020-10-11 02:00:23", "2020-10-17 02:00:23", "2020-10-20 02:00:23", "2020-10-21 02:00:23", "2020-10-22 02:00:23", "2020-10-23 02:00:23", "2020-11-09 02:00:23", "2021-01-08 10:30:23", "2021-01-27 10:30:23",
                                   "2021-08-25 02:00:23 UTC",
                                   "2021-08-26 02:00:23 UTC", "2021-08-30 02:00:23 UTC", "2021-09-03 02:00:23 UTC",
                                   "2021-09-06 02:00:23 UTC", "2021-09-08 02:00:23 UTC", "2021-09-14 02:00:23 UTC",
                                   "2021-09-15 02:00:23 UTC", "2021-09-16 02:00:23 UTC", "2021-09-18 02:00:23 UTC",
                                   "2021-09-19 02:00:23 UTC", "2021-10-18 02:00:23 UTC", "2021-10-30 02:00:23 UTC",
                                   "2021-10-31 02:00:23 UTC", "2021-11-01 02:00:23 UTC", "2021-11-02 02:00:23 UTC",
                                   "2021-11-03 02:00:23 UTC", "2021-11-05 02:00:23 UTC", "2021-11-06 02:00:23 UTC",
                                   "2021-11-07 02:00:23 UTC", "2021-12-23 02:00:23 UTC", "2021-12-26 02:00:23 UTC",
                                   "2021-12-31 02:00:23 UTC", "2022-01-02 02:00:23 UTC", "2022-01-03 02:00:23 UTC",
                                   "2022-01-05 02:00:23 UTC", "2022-01-26 02:00:23 UTC", "2022-01-31 02:00:23 UTC",
                                   "2022-02-02 02:00:23 UTC", "2020-10-23 02:00:23",
                                   "2021-06-01 02:00:23","2021-06-03 02:00:23","2021-06-11 02:00:23","2021-06-14 02:00:23","2021-06-15 02:00:23",
                                   "2022-10-01 02:00:23 UTC", "2022-10-03 02:00:23 UTC",
                                   "2022-10-12 02:00:23 UTC", "2022-10-14 02:00:23 UTC", "2022-10-23 02:00:23 UTC",
                                   "2022-11-11 02:00:23 UTC", "2022-11-13 02:00:23 UTC", "2022-11-20 02:00:41 UTC",
                                   "2020-10-01 02:00:23 UTC", "2020-10-03 02:00:23 UTC", "2020-10-23 02:00:23 UTC",
                                   "2020-10-24 02:00:23 UTC", "2020-11-25 02:00:23 UTC", "2020-11-26 02:00:23 UTC",
                                   "2021-10-20 02:00:23", "2021-11-02 02:00:23", "2022-01-22 02:00:23 UTC", "2022-01-23 02:00:23 UTC", "2022-02-14 02:00:23 UTC",
                                   "2022-02-15 02:00:23 UTC", "2022-02-20 02:00:23 UTC", "2022-02-23 02:00:23 UTC",
                                   "2022-10-12 02:00:23 UTC", "2022-10-15 02:00:23 UTC",
                                   "2022-10-16 02:00:23 UTC", "2022-10-24 02:00:23 UTC", "2022-10-30 02:00:23 UTC",
                                   "2023-01-26 02:00:23 UTC",
                                   "2021-10-01 02:00:23", "2021-10-02 02:00:23", "2022-03-05 02:00:23",
                                   "2021-10-23 02:00:41",
                                   "2021-10-21 02:00:23", "2022-01-16 02:00:23",
                                   "2021-10-10 02:00:23", "2021-10-12 02:00:23", "2021-12-08 02:00:23",
                                   "2021-11-28 02:00:23", "2022-02-03 02:00:23",
                                   "2021-10-31 02:00:23 UTC", "2021-11-01 02:00:23 UTC", "2021-12-11 03:40:08 UTC", "2021-12-12 02:00:23 UTC", "2021-12-13 02:00:23 UTC", "2021-12-17 02:00:23 UTC",
                                   "2021-10-07 02:00:23", "2021-11-04 02:00:23", "2022-06-25 02:00:23", "2022-07-02 02:00:23", "2022-07-03 02:00:23", "2022-07-31 02:00:23", "2022-11-19 02:00:23",
                                   "2021-10-30 02:00:23 UTC",
                                   "2021-10-31 02:00:23 UTC", "2021-11-10 02:00:41 UTC", "2021-11-11 02:00:23 UTC",
                                   "2021-11-12 02:00:23 UTC", "2021-11-23 02:00:23 UTC", "2021-11-25 02:00:23 UTC",
                                   "2021-12-27 02:00:23 UTC", "2022-01-04 02:00:23 UTC", "2022-01-07 02:00:23 UTC",
                                   "2022-01-08 02:00:23 UTC", "2022-01-09 02:00:41 UTC", "2022-01-15 02:00:23 UTC",
                                   "2022-01-19 02:00:23 UTC", "2022-01-20 02:00:23 UTC", "2022-01-22 02:00:23 UTC",
                                   "2022-01-25 02:00:23 UTC", "2022-01-27 02:00:23 UTC", "2022-01-28 02:00:23 UTC",
                                   "2022-01-30 02:00:23 UTC", "2022-02-05 02:00:23 UTC", "2022-02-06 02:00:23 UTC",
                                   "2022-02-07 02:00:23 UTC", "2022-02-10 02:00:23 UTC", "2022-02-11 02:00:23 UTC",
                                   "2022-02-12 02:00:23 UTC", "2022-02-14 02:00:23 UTC", "2022-02-15 02:00:23 UTC",
                                   "2022-02-18 02:00:23 UTC", "2022-02-19 02:00:23 UTC",
                                   "2022-10-01 02:00:23 UTC",
                                   "2022-10-21 02:00:23 UTC", "2022-10-23 02:00:23 UTC", "2022-11-10 02:00:23 UTC",
                                   "2022-11-14 06:20:08 UTC", "2022-11-17 02:00:23 UTC", "2022-11-21 03:30:07 UTC",
                                   "2022-11-22 02:00:23 UTC", "2022-11-23 02:00:23 UTC", "2022-12-05 02:00:23 UTC",
                                   "2022-12-24 02:00:23 UTC", "2022-12-25 02:00:23 UTC", "2022-12-26 02:00:23 UTC",
                                   "2023-01-02 02:00:23 UTC", "2023-01-04 02:00:23 UTC", "2023-01-05 02:00:23 UTC",
                                   "2023-01-06 02:00:23 UTC", "2023-01-07 02:00:23 UTC", "2023-01-11 02:00:23 UTC",
                                   "2023-01-25 02:00:23 UTC", "2023-01-26 02:00:23 UTC"),
                          individual.id =  c(1176031140, 1176031140, 1176031140,
                                             1176034659, 1176034659, 1176034659,
                                             rep(1176037652, times = 32),
                                             rep(1176038499, times = 24),
                                             1176046609,
                                             1176047129, 1176047129,
                                             1176048179, 1176048179, 1176048179, 1176048179,
                                             1176050456,
                                             rep(1176050456, times = 6),
                                             1176054402, 1176054402, 1176054402, 
                                             rep(1176056064, times = 19),
                                             rep(1176056064, times = 19),
                                             1176058760, 1176058760,
                                             rep(1176065189, times = 24),
                                             rep(1178289602, times = 28),
                                             rep(1178289602, times = 45),
                                             rep(1178290125, times = 9),
                                             rep(1178290125, times = 29), 
                                             1183382998, 1183382998, 1183382998, 1183382998, 1183382998, 1183382998, rep(1183382998, times = 8),
                                             rep(1187055990, times = 6), 
                                             1573099352, 1573099352, rep(1573103673, times = 12),
                                             1573153892, 1573153892, 1573153892,
                                             1573164387,
                                             1576763743, 1576763743,
                                             1576770126, 1576770126, 1576770126,
                                             1576777861, 1576777861,
                                             rep(1576778502, times = 6),
                                             rep(1576785208, times = 7),
                                             rep(1576785841, times = 51)))



1176037652
between(timestamp, "2020-11-01 02:00:23", "2021-03-19 02:00:23")

1176038499
timestamp > "2020-09-06 02:00:23"


# 1176060432 is an example of no ldds and perfect 40km accuracy

# example bird 1178289602

ex <- ld %>% 
  filter(individual.id == 1178289602) %>% 
  filter(!date(timestamp) %in% date(corrections$date[corrections$individual.id == 1178289602]))

ggplot(ex, aes(location.long, location.lat, color = timestamp)) +
  borders(xlim = c(-10,10), ylim = c(20, 50)) +
  geom_point() +
  geom_point(data = ex_ldds, color = "red") +
  theme_classic()

ex_ldds <- ld %>% 
  filter(individual.id == 1178289602 & date(timestamp) %in% date(corrections$date[corrections$individual.id == 1178289602]))

multi_ex <- ld %>% 
  filter(individual.id %in% corrections$individual.id)

multi_ex_ldds <- lapply(split(multi_ex, multi_ex$individual.id), function(x){
  corr <- x %>% 
    filter(date(timestamp) %in% date(corrections$date[corrections$individual.id == unique(x$individual.id)]))
}) %>% reduce(rbind)

ggplot(multi_ex, aes(location.long, location.lat, color = timestamp)) +
  borders(xlim = c(-10,10), ylim = c(20, 50)) +
  geom_hline(yintercept = min(multi_ex$location.lat)) +
  geom_point() +
  geom_point(data = multi_ex_ldds, color = "red") +
  theme_classic() +
  facet_wrap(~individual.id)

multi_ex$individual.id_yr <- paste0(multi_ex$individual.id, "_", year(multi_ex$timestamp))

lapply(split(multi_ex, multi_ex$individual.id_yr), function(x){
  fall_x <- x %>% 
    filter(month(timestamp) %in% c(7:11))
  
  if(nrow(fall_x) > 0){
    png(paste0("C:/Users/hbronnvik/Documents/storkSSFs/m_locs_check/fall_", unique(x$individual.id_yr), ".png"))
    print(ggplot(fall_x, aes(location.long, location.lat, color = timestamp)) +
            borders(xlim = c(-10,10), ylim = c(20, 50)) +
            geom_hline(yintercept = min(fall_x$location.lat)) +
            geom_point() +
            geom_point(data = multi_ex_ldds %>% filter(individual.id == unique(fall_x$individual.id)), color = "red", pch = 1) +
            labs(title = unique(fall_x$individual.id)) +
            theme_classic())
    dev.off()
  }
  
  spring_x <- x %>% 
    filter(month(timestamp) %in% c(12:7))
  
  if(nrow(spring_x) > 0){
    png(paste0("C:/Users/hbronnvik/Documents/storkSSFs/m_locs_check/spring_", unique(x$individual.id_yr), ".png"))
    print(ggplot(spring_x, aes(location.long, location.lat, color = timestamp)) +
            borders(xlim = c(-10,10), ylim = c(20, 50)) +
            geom_hline(yintercept = min(fall_x$location.lat)) +
            geom_point() +
            geom_point(data = multi_ex_ldds %>% filter(individual.id == unique(spring_x$individual.id)), color = "red", pch = 1) +
            labs(title = unique(spring_x$individual.id)) +
            theme_classic())
    dev.off()
  }
})

ex <- ex %>% 
  mutate(season = ifelse(month(timestamp) %in% c(7:11), "post", "pre"),
         track_id = paste(season, year(timestamp), individual.id, sep = "_"))

mins <- ex %>% 
  filter(season == "post") %>% 
  group_by(track_id) %>% 
  summarize(lat = min(location.lat),
            long = location.long[which.min(location.lat)],
            ts = timestamp[which.min(location.lat)])

lasts <- ex %>% 
  filter(season == "post") %>% 
  group_by(track_id) %>% 
  filter(timestamp < mins$ts[mins$track_id == unique(track_id)]) %>% 
  slice(n())

d_last <- ex %>% 
  filter(season == "post") %>% 
  group_by(track_id) %>% 
  filter(daily_dist > 100000) %>% 
  slice(n())

sliced_ex <- ex %>% 
  filter(season == "post") %>% 
  group_by(track_id) %>% 
  filter(timestamp <= d_last$timestamp[d_last$track_id == unique(track_id)])

last_displaced <- ex %>% 
  mutate(week = ceiling(as.numeric(difftime(timestamp, head(timestamp, 1), units = "weeks")))) %>% 
  group_by(week) %>% 
  mutate(week_dist = abs(head(location.lat, 1) - tail(location.lat, 1)),
         ld_week = ifelse(abs(head(location.lat, 1) - tail(location.lat, 1)) > 1, T, F)) %>% 
  ungroup() %>% 
  filter(ld_week == T) %>% 
  group_by(track_id) %>% 
  slice(n())
  

ggplot(ex, aes(location.long, location.lat, color = track_id)) +
  borders(xlim = c(-10,10), ylim = c(20, 50))  +
  geom_path() +
  # geom_point(data = mins, aes(long, lat), cex = 3) +
  geom_point(data = last_displaced, cex = 3, pch = 1, color = "black") +
  theme_classic()


# the last south/north ldd of a season below 20, below 35, below 45

lowest <- lapply(clean_locations, function(x){
  lowest <- x %>% 
    mutate(season = ifelse(month(timestamp) %in% c(7:11), "post", "pre"),
           track_id = paste(season, year(timestamp), individual.id, sep = "_")) %>% 
    group_by(track_id) %>% 
    filter(location.lat == min(location.lat)) %>% 
    dplyr::select(location.long, location.lat, timestamp, track_id, season)
}) %>% reduce(rbind) %>% filter(location.long < 15.7)

rect <- data.frame(xmin = rep(-Inf, times = 3), xmax = rep(Inf, times = 3), ymin = c(5, 33, 40), ymax = c(17, 38, 44))
ggplot(lowest, aes(location.long, location.lat, color = season)) +
  borders(xlim = c(-10,10), ylim = c(20, 50))  +
  geom_point() +
  geom_rect(data = rect, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            inherit.aes = F, alpha = 0.5, color = "gray20", linewidth = 0) +
  theme_classic()

## define zones without east-west movement
# if it gets below 17 degrees -- the last day of ldd that is south and below 17
# if it gets below 38 degrees -- the last day of ldd that is south and below 38
# if it gets below 44 degrees -- the last day of ldd below 44
# intermediate birds' last ldd

sub <- c(clean_locations[12:21])
sub_df <- lapply(sub, function(x){
  df <- x %>% 
    dplyr::select(individual.id, timestamp, location.lat, location.long, daily_dist, daily_direction)
  
  df <- df %>% 
    mutate(season = ifelse(month(timestamp) %in% c(7:11), "post", "pre"),
           track_id = paste(season, ifelse(month(timestamp) == 12, year(timestamp + years(1)), year(timestamp)), individual.id, sep = "_"))
  
  migratory <- lapply(split(df, df$track_id), function(y){
    m <- max(y$daily_dist) > d_thresh
    if(m == T){
      return(y)
    }else{
      y <- y %>% 
        slice(1) %>% 
        mutate(track_id = NA)
      return(y)
    }
  }) %>% reduce(rbind)
  
  return(migratory)
}) %>% 
  reduce(rbind) %>% 
  drop_na(track_id) %>% 
  as.data.frame()

x <- sub_df %>% 
  filter(individual.id == 1176038499)

start_zones <- lapply(split(sub_df, sub_df$track_id), function(x){
  if(unique(x$season) == "post"){
    x <- x %>% 
      filter(daily_dist > d_thresh) %>% 
      slice(1)
  }else{
    start_point <- x$location.lat[1]
    if(start_point < 17){
      # this bird started in the Sahel
      # the last time it started moving north south of 17
      pre_start <- y %>% 
        filter(location.lat < 17 & daily_dist > d_thresh) %>% 
        filter(daily_direction < 225 & daily_direction > 135) %>% 
        slice(n())
      # the first time it started moving north after that
      start <- y %>% 
        filter(location.lat < 17 & daily_dist > d_thresh & timestamp > pre_start$timestamp) %>% 
        filter(daily_direction > 315 | daily_direction < 45) %>% 
        slice(1)
    }else{
      if(start_point > 17 & start_point < 36){
        # this bird started in North Africa
        # the last time it started moving north south of 36
        pre_start <- y %>% 
          filter(location.lat < 36 & daily_dist > d_thresh) %>% 
          filter(daily_direction < 225 & daily_direction > 135) %>% 
          slice(n())
        # the first time it started moving north after that
        start <- y %>% 
          filter(location.lat < 36 & daily_dist > d_thresh & timestamp > pre_start$timestamp) %>% 
          filter(daily_direction > 315 | daily_direction < 45) %>% 
          slice(1)
      }else{
        # the first time it started moving north after that
        start <- y %>% 
          filter(location.lat < 36 & daily_dist > d_thresh & timestamp > pre_start$timestamp) %>% 
          filter(daily_direction > 315 | daily_direction < 45) %>% 
          slice(1)
      }
    }
  }
  return(start)
}) %>% reduce(rbind)


# each bird
zones <- lapply(split(sub_df, sub_df$individual.id), function(x){
  x <- x %>% 
    mutate(season = ifelse(month(timestamp) %in% c(7:11), "post", "pre"),
           track_id = paste(season, ifelse(month(timestamp) == 12, year(timestamp + years(1)), year(timestamp)), individual.id, sep = "_"))
  # each track
  z <- lapply(split(x, x$track_id), function(y){
    if(unique(y$season) == "post"){
      print(unique(y$track_id))
      # the point of the season when the animal reached the lowest latitude
      southernmost <- y %>% 
        filter(location.lat == min(location.lat)) %>% 
        dplyr::select(location.lat) %>% 
        deframe()
      if(southernmost < 17){
        # this bird made it to the Sahel
        # the last southwards long-distance displacement
        end <- y %>% 
          filter(daily_dist > d_thresh & between(daily_direction, 135, 225)) %>% 
          slice(n())
      }else{
        if(southernmost > 17 & southernmost < 36){
          # this bird made it to Africa
          # the last southwards ldd before a non-southwards one in Africa (allow it to go north in Spain)
          winter_north <- y %>% 
            filter(daily_direction > 270 | daily_direction < 90) %>% 
            filter(daily_dist > d_thresh & location.lat < 36) %>% 
            slice(1)
          if(nrow(winter_north) > 0){
            end <- y %>% 
              filter(daily_dist > d_thresh & between(daily_direction, 90, 270) & timestamp < winter_north$timestamp & location.lat < 36) %>% 
              slice(n())
          }else{
            end <- y %>% 
              filter(daily_dist > d_thresh & between(daily_direction, 135, 225)) %>% 
              slice(n())
          }
        }else{
          if(southernmost > 36 & southernmost < 40){
            # this bird made it to Southern or Central Spain
            # the last ldd before a northwards one in Spain (don't allow it to go north in Spain)
            winter_north <- y %>% 
              filter(daily_direction > 315 | daily_direction < 45) %>% 
              filter(daily_dist > d_thresh & location.lat < 40) %>% 
              slice(1)
            if(nrow(winter_north) > 0){
              end <- y %>% 
                filter(daily_dist > d_thresh & timestamp < winter_north$timestamp) %>% 
                slice(n())
            }else{
              end <- y %>% 
                filter(daily_dist > d_thresh) %>% 
                slice(n())
            }
          }else{
            # this bird made it to Northern Spain or did not succeed
            # the last ldd (this is a zone where east/west migration is common)
            end <- y %>% 
              filter(daily_dist > d_thresh) %>% 
              slice(n())
          }
        }
      }
  pdf(paste0("C:/Users/hbronnvik/Documents/storkSSFs/m_locs_check/", unique(y$track_id), ".pdf"), 
      width = 8.5, height = 11)
  print(ggplot(y, aes(location.long, location.lat, color = timestamp)) +
          borders(xlim = c(-10,10), ylim = c(20, 50))  +
          geom_point() +
          geom_point(data = end, color = "red") +
          # geom_point(data = winter_north, color = "orange") +
          theme_classic())
  dev.off()
    return(end)
    }else{
      print(unique(y$track_id))
      # the point of the season when the animal reached the highest latitude
      northernmost <- y %>% 
        filter(location.lat == max(location.lat))
### must define a start point because the migration all happens after that
      if(y$location.lat[1] < 17){
        # this bird started in the Sahel
        # the last time it started moving north south of 17
        pre_start <- y %>% 
          filter(location.lat < 17 & daily_dist > d_thresh) %>% 
          filter(daily_direction < 225 & daily_direction > 135) %>% 
          slice(n())
        # the first time it started moving north after that
        start <- y %>% 
          filter(location.lat < 17 & daily_dist > d_thresh & timestamp > pre_start$timestamp) %>% 
          filter(daily_direction > 315 | daily_direction < 45) %>% 
          slice(1)
        
        if(nrow(start) > 0){
          if(northernmost$location.lat < 45){
            # this bird did not return to Germany
            spring_south <- y %>% 
              filter(daily_direction < 225 & daily_direction > 135) %>% 
              filter(daily_dist > d_thresh & timestamp > start$timestamp) %>% 
              slice(n())
            if(nrow(spring_south) > 0){
              end <- y %>% 
                filter(daily_dist > d_thresh & timestamp < spring_south$timestamp) %>% 
                slice(n())
            }else{
              end <- y %>% 
                filter(daily_dist > d_thresh) %>% 
                slice(n())
            }
          }else{
            # this bird did return to Germany
            # the last marginally northward migration day (allow E/W transit through Switzerland)
            end <- y %>% 
              filter(daily_direction > 270 | daily_direction < 90) %>% 
              filter(daily_dist > d_thresh) %>% 
              slice(n())
          }
        }else{
          end <- data.frame(location.long = NA, location.lat = NA, timestamp = NA)
        }
      }else{
        if(y$location.lat[1] > 17 & y$location.lat[1] < 36){
          # this bird started in North Africa
          # the last time it started moving north south of 36
          pre_start <- y %>% 
            filter(location.lat < 36 & daily_dist > d_thresh) %>% 
            filter(daily_direction < 225 & daily_direction > 135) %>% 
            slice(n())
          # the first time it started moving north after that
          start <- y %>% 
            filter(location.lat < 36 & daily_dist > d_thresh & timestamp > pre_start$timestamp) %>% 
            filter(daily_direction > 315 | daily_direction < 45) %>% 
            slice(1)
          
          if(nrow(start) > 0){
            if(northernmost$location.lat < 45){
              # this bird did not return to Germany
              spring_south <- y %>% 
                filter(daily_direction < 225 & daily_direction > 135) %>% 
                filter(daily_dist > d_thresh & timestamp > start$timestamp) %>% 
                slice(n())
              if(nrow(spring_south) > 0){
                end <- y %>% 
                  filter(daily_dist > d_thresh & timestamp < spring_south$timestamp) %>% 
                  slice(n())
              }else{
                end <- y %>% 
                  filter(daily_dist > d_thresh) %>% 
                  slice(n())
              }
            }else{
              # this bird did return to Germany
              # the last marginally northward migration day (allow E/W transit through Switzerland)
              end <- y %>% 
                filter(daily_direction > 270 | daily_direction < 90) %>% 
                filter(daily_dist > d_thresh) %>% 
                slice(n())
            }
          }else{
            end <- data.frame(location.long = NA, location.lat = NA, timestamp = NA)
          }
        }else{
          # this bird started in Spain
          # the last time it started moving north south of 44
          pre_start <- y %>% 
            filter(location.lat < 44 & daily_dist > d_thresh) %>% 
            filter(daily_direction < 225 & daily_direction > 135) %>% 
            slice(n())
          # the first time it started moving north after that
          start <- y %>% 
            filter(location.lat < 44 & daily_dist > d_thresh & timestamp > pre_start$timestamp) %>% 
            filter(daily_direction > 315 | daily_direction < 45) %>% 
            slice(1)
          
          if(nrow(start) > 0){
            if(northernmost$location.lat < 45){
              # this bird did not return to Germany
              spring_south <- y %>% 
                filter(daily_direction < 225 & daily_direction > 135) %>% 
                filter(daily_dist > d_thresh & timestamp > start$timestamp) %>% 
                slice(n())
              if(nrow(spring_south) > 0){
                end <- y %>% 
                  filter(daily_dist > d_thresh & timestamp < spring_south$timestamp) %>% 
                  slice(n())
              }else{
                end <- y %>% 
                  filter(daily_dist > d_thresh) %>% 
                  slice(n())
              }
            }else{
              # this bird did return to Germany
              # the last marginally northward migration day (allow E/W transit through Switzerland)
              end <- y %>% 
                filter(daily_direction > 270 | daily_direction < 90) %>% 
                filter(daily_dist > d_thresh) %>% 
                slice(n())
            }
          }else{
            end <- data.frame(location.long = NA, location.lat = NA, timestamp = NA)
          }
        }
      }
      
      pdf(paste0("C:/Users/hbronnvik/Documents/storkSSFs/m_locs_check/", unique(y$track_id), ".pdf"), 
          width = 8.5, height = 11)
      print(ggplot(y, aes(location.long, location.lat, color = timestamp)) +
              borders(xlim = c(-10,10), ylim = c(20, 50))  +
              geom_point() +
              geom_point(data = end, color = "red") +
              theme_classic())
      dev.off()
    }
  })
})






