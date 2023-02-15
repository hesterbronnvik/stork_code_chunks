### data location to find the names of birds that should be included (ie. have completed more than one migration)
### Hester Brønnvik
### 29.09.2022

library(lubridate)
library(geosphere)
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
visDODs <- readRDS("C:/Users/hbronnvik/Documents/storkSSFs/visual_ACC_deaths_2023-02-01.rds")


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
}) # 13.08803 mins
# DOP duplucates:
# "Snöfrid + A6Y49 (eobs 7971)"
# "Wiesenmann (A6R70, eobs 8027)"
# "Maximilian (E0343, e-obs 8032)"

# saveRDS(clean_locations, file = "C:/Users/hbronnvik/Documents/storkSSFs/clean_locations_2023-02-10.rds")



n <- lapply(clean_locations, function(x){
  df <- data.frame(unique(x$individual.id))
  return(df)
}) %>% reduce(rbind)

clean_locations <- readRDS("C:/Users/hbronnvik/Documents/storkSSFs/clean_locations_2023-02-10.rds")


# take the clean locations, remove birds that did not migrate and add a track ID for the ones that did
migration_locations <- lapply(clean_locations, function(x){
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

# define birds that took eastern routes as ones that are ever east of 15.7 longitude (East Germany)
eastern_birds <- unique(migration_locations$individual.id[migration_locations$location.long > 15.7])
# remove the eastern birds
migration_locations <- migration_locations %>% 
  filter(!individual.id %in% eastern_birds)

# lowest <- migration_locations %>% 
#   group_by(track_id) %>% 
#   filter(location.lat == min(location.lat) | location.lat == max(location.lat))
# rect <- data.frame(xmin = rep(-Inf, times = 3), xmax = rep(Inf, times = 3), ymin = c(5, 30, 40), ymax = c(17, 36, 50))
# ggplot(lowest, aes(location.long, location.lat, color = season)) +
#   borders(xlim = c(-10,10), ylim = c(20, 50))  +
#   geom_point() +
#   geom_rect(data = rect, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
#             inherit.aes = F, alpha = 0.5, color = "gray20", linewidth = 0) +
#   theme_classic()

# find the start times of each track based on season and map zone
start_zones <- lapply(split(migration_locations, migration_locations$track_id), function(x){
  if(unique(x$season) == "post"){
    start <- x %>% 
      filter(daily_dist > d_thresh) %>% 
      slice(1)
  }else{
    start_point <- x$location.lat[1]
    if(start_point < 17){
      # this bird started in the Sahel
      # the last time it started moving north south of 17
      pre_start <- x %>% 
        filter(location.lat < 17 & daily_dist > d_thresh) %>% 
        filter(daily_direction < 225 & daily_direction > 135) %>% 
        slice(n())
      # the first time it started moving north after that
      if(nrow(pre_start)){
        start <- x %>% 
          filter(location.lat < 17 & daily_dist > d_thresh & timestamp > pre_start$timestamp) %>% 
          filter(daily_direction > 315 | daily_direction < 45) %>% 
          slice(1)
      }else{
        start <- x %>% 
          filter(daily_dist > d_thresh) %>% 
          slice(1) %>% 
          mutate(track_id = NA)
      }
    }else{
      if(start_point > 17 & start_point < 36){
        # this bird started in North Africa
        # the last time it started moving north south of 36
        pre_start <- x %>% 
          filter(location.lat < 36 & daily_dist > d_thresh) %>% 
          filter(daily_direction < 225 & daily_direction > 135) %>% 
          slice(n())
        # the first time it started moving north after that
        if(nrow(pre_start)){
        start <- x %>% 
          filter(location.lat < 36 & daily_dist > d_thresh & timestamp > pre_start$timestamp) %>% 
          filter(daily_direction > 315 | daily_direction < 45) %>% 
          slice(1)
        }else{
          start <- x %>% 
            filter(daily_dist > d_thresh) %>% 
            slice(1) %>% 
            mutate(track_id = NA)
        }
      }else{
        # this bird started it's pre-breeding migration in Spain
        start <- x %>% 
          filter(daily_dist > d_thresh) %>% 
          filter(daily_direction > 270 | daily_direction < 90) %>% 
          slice(1)
      }
    }
  }
  return(start)
}) %>% reduce(rbind) %>% 
  drop_na(track_id)

# only use tracks with a start point
migration_locations <- migration_locations %>% 
  filter(track_id %in% start_zones$track_id)

# different parts of the flyway are more or less east/west or north/south
# split the definitions so that in north/south zones (Gibraltar, France) we look for turn points
# whereas in east/west zones (Spain, Sahel) we look for north/south directional movement
end_zones <- lapply(split(migration_locations, migration_locations$track_id), function(y){
  # use different definitions of end for different seasons and map zones
    if(unique(y$season) == "post"){
      print(unique(y$track_id))
      # the point of the season when the animal reached the lowest latitude
      southernmost <- y %>% 
        filter(location.lat == min(location.lat)) %>% 
        dplyr::select(location.lat) %>% 
        slice(1) %>% # in case it is at that point twice
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
  # pdf(paste0("C:/Users/hbronnvik/Documents/storkSSFs/m_locs_check/", unique(y$track_id), ".pdf"), 
  #     width = 8.5, height = 11)
  # print(ggplot(y, aes(location.long, location.lat, color = timestamp)) +
  #         borders(xlim = c(-10,10), ylim = c(20, 50))  +
  #         geom_point() +
  #         geom_point(data = end, color = "red") +
  #         # geom_point(data = winter_north, color = "orange") +
  #         theme_classic())
  # dev.off()
    return(end)
    }else{
      print(unique(y$track_id))
      # the point of the season when the animal reached the highest latitude
      northernmost <- y %>% 
        filter(location.lat == max(location.lat))
        # the first time it started migrating north 
        start <- start_zones %>% 
          filter(track_id == unique(y$track_id))

            # the last marginally northward migration day (allow E/W transit through Spain & Switzerland)
            end <- y %>% 
              filter(daily_direction > 270 | daily_direction < 90 & timestamp > start$timestamp) %>% 
              filter(daily_dist > d_thresh) %>% 
              slice(n())
            return(end)
      }
      
      # pdf(paste0("C:/Users/hbronnvik/Documents/storkSSFs/m_locs_check/", unique(y$track_id), ".pdf"), 
      #     width = 8.5, height = 11)
      # print(ggplot(y, aes(location.long, location.lat, color = timestamp)) +
      #         borders(xlim = c(-10,10), ylim = c(20, 50))  +
      #         geom_point() +
      #         geom_point(data = end, color = "red") +
      #         theme_classic())
      # dev.off()
}) %>% reduce(rbind)

start_zones <- start_zones %>% 
  mutate(point = "start") %>% 
  filter(track_id %in% end_zones$track_id)
end_zones <- end_zones %>% 
  mutate(point = "end")
start_end <- rbind(start_zones, end_zones)

# use only tracks where the animal actually moved north or south by at least l_thresh
start_end <- start_end %>% 
  group_by(track_id)

migration_locations <- migration_locations %>% 
  filter(track_id %in% start_end$track_id)

ml <- lapply(split(migration_locations, migration_locations$track_id), function(x){
  se <- start_end %>% 
    filter(track_id == unique(x$track_id))
  x <- x %>% 
    filter(between(timestamp, se$timestamp[se$point == "start"], se$timestamp[se$point == "end"]))
  return(x)
}) %>% reduce(rbind)

# ggplot(ml, aes(location.long, location.lat, color = track_id)) +
#   borders(xlim = c(-10,10), ylim = c(20, 50))  +
#   geom_point(cex = 0.5, alpha = 0.5) +
#   scale_color_viridis(option = "D", discrete = T) +
#   theme_classic() +
#   theme(legend.position = "none")

## Gappy data -- use or discard?
# 23460528
# pre_2020_891793609
# pre_2022_78031713
# pre_2022_219402518
# pre_2022_1576782003

## just weird
# pre_2017_78031713

# get the track IDs from the migratory birds
rs_ids <- migration_locations %>% 
  group_by(individual.id, track_id) %>% 
  slice(1) %>% 
  dplyr::select(individual.id, track_id) %>% 
  ungroup()
# the total number of migrations attempted and completed by each animal in each season 
meta <- migration_locations %>%
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

# find any animals that slipped through the cracks for visual estimates of DOD
lost <- migration_locations %>% 
  group_by(individual.id) %>% 
  slice(n()) %>% 
  ungroup() %>% 
  dplyr::select(individual.id) %>% 
  filter(!individual.id %in% visDODs$individual_id)

# info <- lapply(studies, function(x){
#   info <- getMovebankAnimals(x, loginStored) %>%
#     filter(sensor_type_id == 653 & individual_id %in% unique(lost$individual.id)) %>% 
#     mutate(study = x)
#   return(info)
# }) %>% reduce(rbind)

found <- data.frame(local_identifier = c("Manfred + / DER A2M42 (e-obs 6995)",   "Pfaffenhofen2 + A2M98 (eobs 7346)",   
                                         "Emma + / DER AV742 (eobs 4367)",       "Victor II + / DER AX381 (eobs 4345)", 
                                         "Sommerwind + / DER AW242 (eobs 4357)", "Nils + / DER AU058 (eobs 3335)",      
                                         "Pieps + / DER AL586 (KN2021)",         "Chris + A8V35 (eobs 8050)",           
                                         "Simba + ABB59 (eobs 9660)"),
                    dod = c("2019-08-19 03:02:00", NA, "2015-08-08 02:58:00", NA, "2012-07-25 09:18:00", "2015-02-04 02:10:00", NA, NA, "2022-12-24 15:30:18"),
                    loss = c(T, NA, T, NA, T, T, NA, NA, NA),
                    comment = c(rep("classified by acc", times = 6), "classified_by gps", "classified by acc", "classified by acc"),
                    individual_id = c(891813394, 1173982416, 80671851, 173198760, 173667185, 24564548, 909029799, 1176059034, 2135830100),
                    study = c(24442409, 24442409, 76367850, 76367850, 76367850, 21231406, 21231406, 1176017658, 1176017658))
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
    filter(track_id == track_id[n()]) %>% 
    mutate(track_displacement = abs(location.lat[1] - location.lat[n()]))
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
  full_dist <- unique(final_track$track_displacement) > l_thresh
  # add a column containing the outcome of the migratory track
  x <- x %>% 
    group_by(track_id) %>% 
    mutate(track_status = ifelse(loss == T & track_id == unique(final_track$track_id) & full_dist  == T, "incomplete", "complete")) %>% 
    ungroup()
  return(x)
}) %>% reduce(rbind)

complete_ml <- migration_locations %>% 
  filter(track_status == "complete")


# lowest <- complete_ml %>%
#   group_by(track_id) %>%
#   filter(location.lat == min(location.lat) | location.lat == max(location.lat))
# rect <- data.frame(xmin = rep(-Inf, times = 3), xmax = rep(Inf, times = 3), ymin = c(5, 30, 40), ymax = c(17, 36, 50))
# ggplot(lowest, aes(location.long, location.lat, color = season)) +
#   borders(xlim = c(-10,10), ylim = c(20, 50))  +
#   geom_point() +
#   geom_rect(data = rect, aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
#             inherit.aes = F, alpha = 0.5, color = "gray20", linewidth = 0) +
#   theme_classic()

saveRDS(migration_locations, file = paste0("C:/Users/hbronnvik/Documents/storkSSFs/migration_locations_40kmCompass_", Sys.Date(), ".rds"))
