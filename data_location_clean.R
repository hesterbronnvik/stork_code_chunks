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
s_thresh <- 30 # m/s

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
      dplyr::select(animal_id, study_id)
    juv <- birds %>% 
      filter(animal_life_stage == "" & grepl("release|adult", animal_comments, ignore.case = T) == F) %>% 
      dplyr::select(animal_id, study_id)
    chicks <- rbind(chicks, juv)
  }else{
    chicks <- birds %>% 
      filter(!grepl("release|adult", birds$animal_comments, ignore.case = T)) %>% 
      dplyr::select(animal_id, study_id)
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

# remove errors and group by season
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
  
  # calculate ground speeds for the full data now that faulty locations are removed
  dist <- lapply(1:nrow(locs_df), function(x){
    d <- distVincentyEllipsoid(c(locs_df$location.long[x], locs_df$location.lat[x]), c(locs_df$location.long[x+1], locs_df$location.lat[x+1]))
    return(d)
  }) %>% unlist()
  
  locs_df <- locs_df %>% 
    mutate(distance = dist,
           timediff = as.numeric(difftime(timestamp, lag(timestamp), units = "secs")),
           gr_speed = distance/timediff) %>%
    # now that the duplicates have been dealt with, remove speed outliers
    filter(gr_speed < s_thresh)
  
  # split the data according to time range
  ind <- locs_df %>% 
    filter(!is.na(location.lat)) %>% 
    mutate(phase = ifelse(month(timestamp) %in% c(1:7), paste0("pre-breeding_", year(timestamp)), 
                          ifelse(month(timestamp) == 12, paste0("pre-breeding_", year(timestamp + years(1))),
                                 paste0("post-breeding_", year(timestamp)))))%>% 
    mutate(date = date(timestamp)) %>%
    group_by(date) %>% 
    # the Haversine distance between first and last locations of the day
    mutate(daily_dist = distVincentyEllipsoid(c(head(location.long,1), head(location.lat, 1)), c(tail(location.long,1), tail(location.lat, 1))),
           # the rhumbline bearing between first and last locations of the day
           daily_direction = bearingRhumb(c(head(location.long,1), head(location.lat, 1)), c(tail(location.long,1), tail(location.lat, 1))),
           # finally a binary category denoting migratory or not, it is unnecessary but simplifies code
           compass_direction = ifelse(daily_direction > 90 & daily_direction < 270, "southward", "northward")) %>% 
    ungroup() %>% 
    arrange(timestamp)
  # 
  # ggplot(ind, aes(location.long, location.lat, color = phase)) +
  #   borders(xlim = c(-5, 5), ylim = c(30, 48)) +
  #   geom_point() +
  #   theme_classic()
  print(ind$individual.local.identifier[1])
  return(ind)
})
# "Snöfrid + A6Y49 (eobs 7971)"
# saveRDS(clean_locations, file = "clean+locations_2023-02-07.rds")
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
    
    df <- x %>% 
      mutate(seq_15 = round_date(timestamp, unit = "15 minutes")) %>% 
      group_by(seq_15) %>% 
      slice(1) %>% 
      ungroup() %>% 
      filter(daily_dist > d_thresh) %>% 
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
             unit_lead = tail(gap_lead, 1)) %>% 
      ungroup()
    # if there is migration:
    if(length(unique(df$date)) > 1){
      # if there are multiple bouts separated by stopovers,
      if(length(unique(df$check)) > 2){
        # the first time a bout of migration is within six weeks of another
        on <- df$timestamp[which(df$unit_lead < w_thresh)[1]]
        # the last time a bout is within six weeks of another
        off <- df$timestamp[which(df$unit_lag > w_thresh)[length(which(df$unit_lag > w_thresh))]]
      }else{
        # if there are only one or two bouts, use the one with the longest duration (most non-consecutive days)
        dur0 <- length(unique(df$date[df$check == 0]))
        dur1 <- length(unique(df$date[df$check == 1]))
        if(dur1 > dur0){
          on <- df$timestamp[which(df$check == 1)][1]
          off <- df$timestamp[which(df$check == 1)][nrow(df[which(df$check == 1),])]
        }else{
          on <- df$timestamp[which(df$check == 0)][1]
          off <- df$timestamp[which(df$check == 0)][nrow(df[which(df$check == 0),])]
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
    saveRDS(migration, file = paste0("C:/Users/hbronnvik/Documents/storkSSFs/classified_gps/migration_locations_40km15min_", paste(unique(x$individual.id), unique(x$phase), sep = "_"), "_", Sys.Date(), ".rds"))
    return(migration)
  })# %>% reduce(rbind)
  return(m)
})# %>% reduce(rbind)
Sys.time() - start_time
saveRDS(migration_locations, file = "C:/Users/hbronnvik/Documents/storkSSFs/migration_locations_40km15min_20230206.rds")

ggplot(migration_locations[[81]][[3]] %>% drop_na(track_id), aes(location.long, location.lat, color = track_id)) +
  borders(xlim = c(-5, 5), ylim = c(30, 50)) +
  geom_point() +
  theme_classic() +
  facet_wrap(~track_id)

x <- readRDS("C:/Users/hbronnvik/Documents/storkSSFs/classified_gps/migration_locations_40km15min_1578616761_pre-breeding_2021_2023-02-07.rds")
