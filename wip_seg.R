### data location to find the names of birds that should be included (ie. have completed more than one migration)
### Hester Brønnvik
### 29.09.2022

library(lubridate)
library(geosphere)
library(move)
library(stringr)
library(tidyverse)
library(data.table)
library(mapview)
d_thresh <- 40000 # meters
w_thresh <- 6 # weeks
s_thresh <- 30 # days
l_thresh <- 3 # degrees latitude
g_thresh <- 7 # days (allowable gap in transmission)

# required information
load("C:/Users/hbronnvik/Documents/loginStored.rdata")
studies <- c(24442409, 212096177, 76367850, 21231406, 1176017658, 173641633)
visDODs <- readRDS("C:/Users/hbronnvik/Documents/storkSSFs/visual_ACC_deaths_2023-04-11.rds")


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

# find data since the last collection using date
# update_info <- lapply(studies, function(x){
#   df <- getMovebankAnimals(x, loginStored) %>% 
#     filter(number_of_events > 0) %>% 
#     mutate(timestamp_end = as.POSIXct(sub("\\.000", "", timestamp_end), tz = "UTC"),
#            study_id = x) %>% 
#     filter(sensor_type_id == 653 & timestamp_end > as.POSIXct("2023-01-01 01:00:00", tz = "UTC"))
# }) %>% 
#   reduce(rbind) %>% 
#   filter(individual_id %in% info$animal_id)
# 
# update <- lapply(1:nrow(update_info), function(x){
#   ff <- getMovebankLocationData(update_info$study_id[x], 653, update_info$individual_id[x], loginStored)
#   saveRDS(ff, file = paste0("C:/Users/hbronnvik/Documents/storkSSFs/full_data/", unique(ff$individual.id), "_", Sys.Date(), ".rds"))
#   ff
# })

# download data
lapply(1:nrow(info), function(x){
  df <- getMovebankLocationData(info$study_id[x], sensorID = 653, 
                                animalName = info$animal_id[x], login = loginStored)
  saveRDS(df, file = paste0("C:/Users/hbronnvik/Documents/storkSSFs/full_data/", info$animal_id[x], "_24052023.rds"))
})

# select the nestling data files
full_files <- list.files("C:/Users/hbronnvik/Documents/storkSSFs/full_data", pattern = "_24052023.rds", full.names = T)
# todo <- sapply(full_files, function(x){
#   td <- gsub("data/", "", str_split(full_files[[1]], "_")[[1]][2])
#   td
# })
# todo <- todo[todo %in% info$animal_id]
# # todo <- todo[todo %in% update_info$individual_id]
# 
# full_files <- lapply(1:length(full_files), function(x){
#   name <- str_sub(full_files[x], 50, -16)
#   if(name %in% todo){
#     return(full_files[x])
#   }
# }) %>% unlist()

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
}) # 13.011 mins
Sys.time() - start_time
# DOP duplucates:
# "Snöfrid + A6Y49 (eobs 7971)"
# "Wiesenmann (A6R70, eobs 8027)"
# "Maximilian (E0343, e-obs 8032)"

# saveRDS(clean_locations, file = paste0("C:/Users/hbronnvik/Documents/storkSSFs/clean_locations_", Sys.Date(),".rds"))

# compress the updated data
# updated_cl <- rbindlist(clean_locations, fill = T)
# # access the old data
# clean_locations <- readRDS("C:/Users/hbronnvik/Documents/storkSSFs/clean_locations_2023-02-15.rds")
# # get the indices in the old data of the updated data
# updated <- lapply(1:length(clean_locations), function(x){
#   id <- unique(clean_locations[[x]]$individual.id)
#   if(id %in% unique(updated_cl$individual.id)){return(x)}
# }) %>% unlist()
# # remove the data from the updated individuals
# clean_locations <- clean_locations[-updated]
# updated_cl <- split(updated_cl, updated_cl$individual.id)
# # add the updates
# clean_locations <- c(clean_locations, updated_cl)
# saveRDS(clean_locations, file = paste0("C:/Users/hbronnvik/Documents/storkSSFs/clean_locations_", Sys.Date(),".rds"))

# cl <- lapply(clean_locations, function(x){
#   cl <- x %>%
#     filter(timestamp > "2022-10-01 01:00:00")
#   if(max(cl$location.lat) > 47){
#     return(cl)
#   }
# })
# cl <- cl[lapply(cl,length)>0] 
# sp <- cl[[6]]
# coordinates(sp) <- ~location.long+location.lat
# proj4string(sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
# mapview(sp, zcol = "individual.id")
# lapply(clean_locations, function(ind){
#   png(paste0("C:/Users/hbronnvik/Documents/storkSSFs/sp_tracks/", unique(ind$individual.id), ".png"), width = 11, height = 8.5, units = "in", res = 100)
#   print(ggplot(ind %>% filter(timestamp > as.POSIXct("2022-12-01 01:00:00", tz = "UTC")), aes(location.long, location.lat, color = timestamp)) +
#           borders("world", xlim = c(-20, 20), ylim = c(20, 50)) +
#           geom_point() +
#           labs(y = "Latitude", x = "Longitude") +
#           theme_classic() +
#           theme(legend.position = "none",
#                 text = element_text(size = 20),
#                 axis.text = element_text(color = "black")))
#   dev.off()
# })
# ff <- rbindlist(clean_locations, fill = T)
# x <- ff %>%
#   filter(individual.id == 173199625) %>%
#   as.data.frame()

# sp <- a_data[1:1000,]
# sp <- x %>%
#   filter(year(timestamp) %in% c(2019)) %>%
#   mutate(yr = as.factor(year(timestamp)), mo = as.factor(month(timestamp)))
# sp <- sp %>% mutate(yr_mo = paste0(yr, mo),
#                     ts = as.character(timestamp),
#                     dd = as.factor(daily_dist))
# coordinates(sp) <- ~location.long+location.lat
# proj4string(sp) <- "EPSG:4326"
# lines <- SpatialLines(list(Lines(list(Line(coordinates(sp))), ID = "a")))
# mapView(lines, zcol = "trackID")
# 
# ggplot(x %>% filter(year(timestamp) == 2016) %>%  mutate(yr = as.factor(year(timestamp))), aes(location.long, location.lat)) +
#   borders("world", xlim = c(-10, 20), ylim = c(0, 60)) +
#   geom_point() +
#   theme_classic() +
#   facet_wrap(~trackID)
# 
# library(gganimate)
# library(viridis)
# start <- x %>% 
#   filter(daily_dist > d_thresh) %>% 
#   slice(1)
# plotting <- x %>% 
#   filter(timestamp > start$timestamp) %>% 
#   arrange(timestamp)
# map_with_animation <- ggplot(plotting, aes(location.long, location.lat, color = trackID)) +
#   borders("world", xlim = c(-10, 10), ylim = c(0, 60), colour = "black") +
#   geom_point(alpha = 0.5, size = 2) +
#   # scale_color_viridis(option = "A", discrete = T) +
#   theme_classic() +
#   theme(legend.position = "none",
#         text = element_text(size = 20, color = "black"),
#         axis.line = element_line(color = "black"),
#         axis.text = element_text(size = 20, color = "black")) +
#   ggtitle('Year: {frame_time}') +
#   labs(x = "Longitude", y = "Latitude") +
#   transition_time(timestamp) +
#   shadow_mark()
# animate(map_with_animation, nframes = 75, height = 1500, width = 1000)

clean_locations <- readRDS("C:/Users/hbronnvik/Documents/storkSSFs/clean_locations_2023-04-10.rds")

ld_locs <- lapply(clean_locations, function(x){tryCatch({
  print(x$individual.id[1])
  # print(x$individual.local.identifier[1])
  df <- x %>% 
    dplyr::select(individual.id, timestamp, location.lat, location.long, daily_dist, daily_direction, ground_speed_15) %>% 
    filter(daily_dist > d_thresh) %>%
    arrange(timestamp) %>%
    mutate(timeLag = as.numeric(timestamp - lag(timestamp), units = "secs"),
           # insert an enormous time lag for the first location to remove the NA
           timeLag = ifelse(is.na(timeLag), 1e5, timeLag),
           newBurst = ifelse(round(timeLag) <= weeks(1), F, T),
           newCluster = ifelse(round(timeLag) <= weeks(6), F, T),
           # take the cumulative sum to act as a unique ID for each burst
           cumu_check_for_burst = cumsum(newBurst),
           cumu_check_for_clust = cumsum(newCluster)) %>% 
    group_by(cumu_check_for_burst) %>% 
    # for the bursts, calculate the time difference between the last and first locations
    mutate(burstLength = length(unique(date(timestamp))),
           # add an ID to each burst
           burstID = as.character(cur_group_id())) %>%
    # remove the bursts that do not meet user-set criteria
    # filter(burstLength > MinBurstLength) %>% 
    ungroup() %>% 
    group_by(cumu_check_for_clust) %>% 
    # for the bursts, calculate the time difference between the last and first locations
    mutate(clusterLength = length(unique(date(timestamp))),
           # add an ID to each burst
           clusterID = as.character(cur_group_id())) %>% 
    ungroup() %>% 
    # clean up the sorting columns
    dplyr::select(-"newBurst", -"cumu_check_for_burst", -"newCluster", -"cumu_check_for_clust")
  
  # sp <- x
  # sp$check <- ifelse(sp$timestamp %in% df$timestamp, "yes", "no")
  # coordinates(sp) <- ~location.long+location.lat
  # proj4string(sp) <- "EPSG:4326"
  # mapView(sp, zcol = "check")

  # fpt <- lapply(split(df, df$burstID), function(f){
  #   bd <- unique(f$burstID)
  #   f <- move(f$location.long, f$location.lat, f$timestamp, f, "EPSG:4326")
  #   f <- spTransform(f, center = T)
  #   fpt <- fpt(as(f, "ltraj"), r = 20000, units = "days")
  #   res <- data.frame(FPT_days = fpt[[1]]$r1, date = as.POSIXct(attr(fpt[[1]], "date")))
  #   d <- res %>% 
  #     filter(FPT_days < 1) %>% 
  #     slice(n()) %>% 
  #     mutate(burstID = bd)
  # }) %>% reduce(rbind)
  
  # embc <- lapply(split(df, df$clusterID), function(y){
  #   start <- y %>% 
  #     arrange(timestamp) %>% 
  #     filter(daily_dist > d_thresh) %>% 
  #     slice(1) %>% 
  #     dplyr::select(timestamp) %>% 
  #     deframe()
  #   end <- y %>% 
  #     arrange(timestamp) %>% 
  #     filter(daily_dist > d_thresh) %>% 
  #     slice(n()) %>% 
  #     dplyr::select(timestamp) %>% 
  #     deframe()
  #   expth <- x %>% 
  #     filter(between(timestamp, start, end)) %>% 
  #     dplyr::select(timestamp, location.long, location.lat, individual.id, daily_dist, daily_direction, ground_speed_15, burstID) %>% 
  #     as.data.frame()
  #   mybcp <- stbc(expth, info=-1)#, smth = 18)
  #   # sctr(mybcp)
  #   # view(mybcp)
  #   expth <- expth %>% 
  #     mutate(embc = as.factor(mybcp@A))
  #   sp <- expth %>% mutate(burstID = unique(y$clusterID))#%>% filter(burstID %in% unique(burstID[which(embc == 3)]))
  #   sp
  # 
  # }) %>% rbindlist()
  # sp <- embc %>%
  #   filter(burstID %in% unique(burstID[which(embc == 3)]))
  # # sp$speed <- ifelse(sp$daily_dist > d_thresh, "high", "low")
  # coordinates(sp) <- ~location.long+location.lat
  # proj4string(sp) <- "EPSG:4326"
  # mapView(sp, zcol = "burstID")

  class_df <- df %>%  #embc %>% 
    group_by(burstID) %>% 
    mutate(burst_angle = ifelse(location.lat[1] > location.lat[n()], "south", "north"),
           burst_season = ifelse(month(timestamp[n()]) %in% c(8:11), "late", 
                                 ifelse(month(timestamp[n()]) %in% c(1:6), "early", "change")),
           burst_class = ifelse(burst_season == "early", "spring", 
                                ifelse(burst_season == "late", "fall",
                                       ifelse(burst_angle == "north" & burst_season == "change", "spring", 
                                              ifelse(burst_angle == "south" & burst_season == "change", "fall", NA)))),
           burst_year = ifelse(month(timestamp[n()]) == 12 & burst_class == "spring", year(years(1) + timestamp[n()]), year(timestamp[n()])),
           trackID = paste(unique(individual.id), burst_class, burst_year, sep = "_")) %>% 
    ungroup()#%>% 
    # mutate(timeLag = difftime(timestamp, lag(timestamp), units = "weeks"),
    #        newTrack = ifelse(round(timeLag) <= weeks(6), F, T)) %>% 
    # drop_na(timeLag) %>% 
    # mutate(cumu_check_for_track = cumsum(newTrack)) %>% 
    # group_by(cumu_check_for_track) %>% 
    # mutate(trackLength = length(unique(date(timestamp))),
    #        track_season = ifelse(month(timestamp[n()]) %in% c(8:12), "fall", "spring"),
    #        trackID = paste(year(timestamp[n()]), track_season, as.character(cur_group_id()), sep = "_")) %>%
    # ungroup() %>% 
    # dplyr::select(-"newTrack", -"cumu_check_for_track")

  # sp <- class_df
  # coordinates(sp) <- ~location.long+location.lat
  # proj4string(sp) <- "EPSG:4326"
  # mapView(sp, zcol = "burstID")
  

  # remove ragged edges of tracks the bird survived
  # if the bird died on a ragged edge,  
  
  # identify ragged edges as bursts that do not have high speed days in them
  ragged <- class_df %>% 
    group_by(burstID) %>% 
    mutate(high_speed = max(daily_dist > 70000)) %>%
    ungroup() %>% 
    filter(high_speed == F)
  # did the bird die?
  dod <- visDODs %>% 
    filter(individual_id == unique(x$individual.id)) %>% 
    dplyr::select(dod) %>% 
    deframe()
  if(nrow(ragged) > 0){
    if(!is.na(dod)){
      # did the bird die in a ragged edge?
      cutoff <- date(dod) %in% unique(date(ragged$timestamp))
      if(cutoff == T){
        # how ragged?
        messy <- lapply(2:length(class_df$burstID), function(z){
          temp <- class_df %>%
            filter(burstID == z)
          temp_prev <- class_df %>%
            filter(burstID == z-1)
          temp$prev_long <- !unique(temp_prev$burstID) %in% unique(ragged$burstID)
          temp
        }) %>% reduce(rbind)
        # these are preceded by high speed bursts, but are low speed
        messy <- messy %>% 
          filter(prev_long == T & burstID %in% unique(ragged$burstID))
        # did the bird die in a ragged edge adjacent to a high speed burst?
        keep <- date(dod) %in% unique(date(messy$timestamp))
        if(keep == T){
          edge <- unique(class_df$burstID[date(class_df$timestamp) == date(dod)])
          # discard ragged edges that do not contain a death
          class_df <- class_df %>% 
            filter(!burstID %in% ragged$burstID | burstID == edge)
        }
      }else{
        class_df <- class_df %>% 
          filter(!burstID %in% ragged$burstID)
      }
    }else{
      class_df <- class_df %>% 
        filter(!burstID %in% ragged$burstID)
    }
  }
  
  # sp <- class_df
  # coordinates(sp) <- ~location.long+location.lat
  # proj4string(sp) <- "EPSG:4326"
  # mapView(sp, zcol = "trackID")
  
  # identify ragged edges as recurring in the radius of 20km from the center of the burst
  # library(recurse)
  # rec <- lapply(split(class_df, class_df$burstID), function(f){
  #   bd <- unique(f$burstID)
  #   # the center of a burst
  #   poi <-  geosphere::centroid(as.matrix(cbind(f$location.long, f$location.lat)))
  #   f <- move(f$location.long, f$location.lat, f$timestamp, f, "EPSG:4326")
  #   f <- spTransform(f, center = T)
  #   # visit <- getRecursions(f, 40000, timeunits = "days")
  #   visit_to <- getRecursionsAtLocations(f, as.data.frame(poi), 100000, timeunits = "days")
  #   # r <- visit$revisitStats
  #   r2 <- visit_to$revisitStats
  #   # hist(visit$revisits, breaks = 20, xlab = "Revisits (radius = 40k)", main = bd)
  #   # plot(visit, f, main = bd)
  #   # r <- data.frame(rev = visit$revisits, id = bd)
  #   r2 <- r2 %>% 
  #     mutate(id = bd)
  #   r2
  # }) %>% reduce(rbind)
  # 
  # ggplot(rec, aes(as.factor(id), timeInside)) +
  #   geom_boxplot() +
  #   theme_classic()
  # 
  # 
  # if(!is.na(dod)){
  #   # find short bursts following other short bursts
  #   t <- lapply(2:length(class_df$burstID), function(z){
  #     temp <- class_df %>%
  #       filter(burstID == z)
  #     temp_prev <- class_df %>%
  #       filter(burstID == z-1)
  #     temp$prev_long <- unique(temp_prev$burstLength) > 1
  #     temp %>% 
  #       filter(burstLength == 1 & prev_long == F)
  #   }) %>% reduce(rbind)
  #   # did the bird die in a short burst?
  #   cutoff <- date(dod) %in% unique(date(class_df$timestamp[which(class_df$burstLength ==1)]))
  #   if(cutoff == T){
  #     # if that burst was preceded by a long burst, keep it, else discard short bursts
  #     keep <- !unique(class_df$burstID[which(date(class_df$timestamp)== date(dod))]) %in% unique(t$burstID)
  #     if(keep == T){
  #       class_df <- class_df %>% 
  #         filter(!burstID %in% t$burstID | burstID == unique(class_df$burstID[which(date(class_df$timestamp)== date(dod))]))
  #     }
  #   }
  #   
  # }
  # 
  # et <- lapply(split(class_df, class_df$trackID), function(y){
  #   expth <- y %>% 
  #     # filter(between(timestamp, start, end)) %>% 
  #     dplyr::select(timestamp, location.long, location.lat, individual.id, daily_dist, daily_direction, ground_speed_15, burstID, trackID) %>% 
  #     as.data.frame()
  #   mybcp <- stbc(expth, info=-1, smth = 12)
  #   sctr(mybcp)
  #   view(mybcp)
  #   expth <- expth %>% 
  #     mutate(embc = as.factor(mybcp@A))
  #   sp <- expth %>% 
  #     filter(burstID %in% unique(burstID[which(embc == 3)]))
  #   sp
  #   
  # }) %>% rbindlist()
  # 
  # sp <- et
  # coordinates(sp) <- ~location.long+location.lat
  # proj4string(sp) <- "EPSG:4326"
  # mapView(sp, zcol = "trackID")
  # 
  # fpt <- lapply(split(test, test$burstID), function(f){
  #   bd <- unique(f$burstID)
  #   f <- move(f$location.long, f$location.lat, f$timestamp, f, "EPSG:4326")
  #   f <- spTransform(f, center = T)
  #   fpt <- fpt(as(f, "ltraj"), r = 1000, units = "hours")
  #   res <- data.frame(FPT_days = fpt[[1]]$r1, date = as.POSIXct(attr(fpt[[1]], "date")))
  #   # d <- res %>%
  #   #   filter(FPT_days < 1) %>%
  #   #   slice(n()) %>%
  #   #   mutate(burstID = bd)
  #   res
  # }) %>% reduce(rbind)
  
  tracks <- class_df %>%
    group_by(trackID) %>% 
    slice(1, n())
  tracks <- lapply(split(tracks, tracks$trackID), function(c){
    ID <- unique(c$trackID)
    track <- x %>% 
      filter(between(timestamp, c$timestamp[1], c$timestamp[2])) %>% 
      mutate(trackID = ID)
    track
  }) %>% rbindlist()
  
  # ggplot(c, aes(location.long, location.lat, color = trackID)) +
  #   borders("world", xlim = c(-10, 10), ylim = c(0, 60), colour = "black") +
  #   geom_point() +
  #   theme_classic() +
  #   theme(legend.position = "none",
  #         text = element_text(size = 20, color = "black"),
  #         axis.line = element_line(color = "black"),
  #         axis.text = element_text(size = 20, color = "black")) +
  #   facet_wrap(~trackID)
  
  x <- x %>% 
    rowwise() %>% 
    mutate(trackID = ifelse(timestamp %in% tracks$timestamp, tracks$trackID[which(tracks$timestamp == timestamp)], NA)) %>% 
    ungroup()

  # sp <- x %>% filter(trackID == "1178289602_spring_2021")
  # coordinates(sp) <- ~location.long+location.lat
  # proj4string(sp) <- "EPSG:4326"
  # mapView(sp, zcol = "burstID")
  
  # df <- df %>% 
  #   rowwise() %>% 
  #   mutate(FPT = as.POSIXct(ifelse(burstID %in% unique(fpt$burstID), 
  #                                  fpt$date[which(fpt$burstID == burstID)], NA), 
  #                           origin = "1970-01-01", tz = "UTC")) %>% 
  #   group_by(burstID) %>% 
  #   mutate(move = between(timestamp, timestamp[1], unique(FPT)),
  #          route = paste(burstID, between(timestamp, timestamp[1], unique(FPT)), sep = "_"),
  #          bd = paste(burstID, date(unique(FPT)), sep = "_")) %>% 
  #   ungroup()
  # base <- x %>% 
  #   filter(is.na(trackID))
  # plotting <- lapply(split(x, x$trackID), function(y){
  #   p <- y %>% 
  #     mutate(track_yr = year(timestamp),
  #            season = ifelse(grepl("fall", trackID), "fall", "spring"))
  #   b <- base %>% 
  #     filter(between(timestamp, min(y$timestamp)-days(60), max(y$timestamp)+days(60))) %>% 
  #     mutate(season = unique(p$season),
  #            track_yr = unique(p$track_yr))
  #   p <- p %>% 
  #     rbind(b) %>% 
  #     arrange(timestamp)
  #   p
  # }) %>% reduce(rbind)
  # 
  # 
  # ends <- plotting %>% 
  #   drop_na(trackID) %>% 
  #   group_by(trackID) %>% 
  #   slice(n()) %>% 
  #   group_by(trackID) %>% 
  #   mutate(status = unique(migration_locations$track_status[which(migration_locations$trackID == unique(trackID))]),
  #          pch = ifelse(status == "complete", 8, 4)) %>% 
  #   ungroup()
  # 
  # png(filename = paste0("C:/Users/hbronnvik/Documents/storkSSFs/diagnostic_segmentation_plots/speed_40x70/speed_", unique(df$individual.id), ".png"),
  #     width = 8, height = 11.5, units = "in", res = 1000)
  # print(ggplot(plotting, aes(location.long, location.lat, color = trackID))+
  #         borders("world", xlim = c(-10, 20), ylim = c(30, 50), colour = "black")+
  #         geom_point(cex = 0.1) +
  #         geom_point(data = ends, pch = ends$pch, color = "black", cex = 2) +
  #         theme_classic() +
  #         facet_wrap(~track_yr+season))
  # dev.off()

  # sp <- x
  # coordinates(sp) <- ~location.long+location.lat
  # proj4string(sp) <- "EPSG:4326"
  # mapView(sp, zcol = "trackID")
  # 
  # m <- df %>% 
  #   filter(move == T) %>% 
  #   group_by(burstID) %>% 
  #   mutate(season = ifelse(month(timestamp[n()]) %in% c(8:12), "fall", "spring"),
  #          trackID = paste(individual.id, season, year(timestamp[n()]), burstID, sep = "_")) %>% 
  #   ungroup()
  # m
  return(x)
  }, error = function(e){
    print(geterrmessage())
    if(geterrmessage() == "Must request at least one colour from a hue palette."){dev.off()}
    })
})

migration_locations <- ld_locs[lapply(ld_locs, length) > 1]
migration_locations <- rbindlist(migration_locations, fill = T)
migration_locations <- migration_locations %>% 
  filter(!is.na(trackID))
# get the track IDs from the migratory birds
rs_ids <- migration_locations %>% 
  group_by(individual.id, trackID) %>% 
  slice(1) %>% 
  dplyr::select(individual.id, trackID) %>% 
  ungroup()
# the total number of migrations attempted and completed by each animal in each season 
meta <- migration_locations %>%
  mutate(season = ifelse(grepl("fall", trackID), "fall", "spring")) %>% 
  group_by(individual.id, season) %>% 
  count(trackID) %>% 
  mutate(journey_number = row_number()) %>% 
  ungroup() %>% 
  group_by(individual.id) %>%
  mutate(total_journeys = n()) %>% 
  ungroup()
# add the number of journeys
migration_locations <- migration_locations %>% 
  rowwise() %>% 
  mutate(journey_number = meta$journey_number[which(meta$trackID == trackID)]) %>% 
  ungroup()

# find any animals that slipped through the cracks for visual estimates of DOD
lost <- migration_locations %>% 
  group_by(individual.id) %>% 
  slice(n()) %>% 
  ungroup() %>% 
  dplyr::select(individual.id) %>% 
  filter(!individual.id %in% visDODs$individual_id)

# define birds that took eastern routes as ones that are ever east of 16.5 longitude (East Germany)
eastern_birds <- unique(migration_locations$individual.id[migration_locations$location.long > 16.5])
# remove the eastern birds
migration_locations <- migration_locations %>% 
  filter(!individual.id %in% eastern_birds) %>% 
  group_by(trackID) %>% 
  mutate(track_displacement = abs(location.lat[1] - location.lat[n()])) %>% 
  ungroup() 

# write out the file to use for social information estimates
# this includes all western migrants regardless of success or the length or 
# transmission rates of the tracks, but holds only locations that are during migration
# 307 individuals, 776 tracks
saveRDS(migration_locations, file = "C:/Users/hbronnvik/Documents/storkSSFs/soc_info_migration_locations_2023_05_02.rds")

# now filter for the route selection analysis:
migration_locations <- migration_locations %>%
  filter(track_displacement > l_thresh)

# define the tracks of birds for which data transmission errors prevent determining the track target
# some are also birds with only one or two locations per day, which precludes the possibility
# of accurate calculations of daily displacement and therefore of migration speeds, which should have been
# removed in the data cleaning step but were missed
safe <- migration_locations %>% 
  filter(trackID == "80439771_fall_2015")
# sp <- migration_locations %>% 
#   filter(trackID == "1178289602_spring_2021") %>% 
#   mutate(ts = factor(timestamp))
# coordinates(sp) <- ~location.long+location.lat
# proj4string(sp) <- "EPSG:4326"
# mapView(sp, zcol = "ts")
ind_gaps <- c(1578604607, 23460528, 80436240, 1576790450, 80439771, 504440981)
track_gaps <- c("219392149_spring_2016", "293986322_spring_2018", "1176046609_spring_2022",
               "909029794_fall_2014", "1576782003_spring_2022", "173659608_spring_2018",
               "78031713_fall_2022", "23463463_spring_2018", "78031713_fall_2022", "80438400_spring_2019",
               "909029800_fall_2016", "1576782003_spring_2022", "23463463_fall_2017")

# 288 individuals, 690 tracks
migration_locations <- migration_locations %>% 
  filter(!individual.id %in% ind_gaps & !trackID %in% track_gaps) %>% 
  # add the one useful track from individual 1178289602
  rbind(safe)

# find data since the last collection using date
info <- lapply(studies, function(x){
  df <- getMovebankAnimals(x, loginStored) %>% 
    filter(number_of_events > 0) %>% 
    mutate(timestamp_end = as.POSIXct(sub("\\.000", "", timestamp_end), tz = "UTC"),
           study_id = x) %>% 
    filter(sensor_type_id == 653)
}) %>% 
  reduce(rbind) %>% 
  filter(individual_id %in% migration_locations$individual.id)

# find the success/failure of each migration track
migration_locations <- lapply(split(migration_locations, migration_locations$trackID), function(x){
  # take the migratory route
  track <- x %>% 
    arrange(timestamp) 
  # if the last GPS time point is greater than the time the data were downloaded, use the download date
  t_end <- ifelse(info$timestamp_end[info$individual_id == unique(x$individual.id)] > "2023-04-10", 
                  "2023-04-10 00:00:00", info$timestamp_end[info$individual_id == unique(x$individual.id)])
  # take either the confirmed DOD or the last transmitted time stamp
  loss_time <- visDODs %>% 
    filter(individual_id == unique(x$individual.id)) %>% 
    mutate(dod = as.POSIXct(ifelse(is.na(dod), t_end, dod), tz = "UTC", origin = "1970-01-01")) %>% 
    dplyr::select(dod) %>% 
    deframe()
  # compare the DOD to the end of the migration
  loss <- max(track$timestamp) > loss_time - days(s_thresh)
  # add a column containing the outcome of the migratory track
  x <- x %>% 
    mutate(track_status = ifelse(loss == T & trackID == unique(track$trackID), "incomplete", "complete"))
  return(x)
}) %>% reduce(rbind)

complete_ml <- migration_locations %>% 
  filter(track_status == "complete")

saveRDS(migration_locations, file = paste0("C:/Users/hbronnvik/Documents/storkSSFs/migration_locations_40km70km30daySpeed_", Sys.Date(),".rds"))

### ----------------------------------------------------------------------------------------------


# take the clean locations, remove birds that did not migrate and add a track ID for the ones that did
migration_locations <- lapply(clean_locations, function(x){
  df <- x %>% 
    dplyr::select(individual.id, timestamp, location.lat, location.long, daily_dist, daily_direction, ground_speed_15)
  
  df <- df %>% 
    mutate(season = ifelse(month(timestamp) %in% c(7:11), "post", "pre"),
           track_id = paste(season, ifelse(month(timestamp) == 12, year(timestamp + years(1)), year(timestamp)), individual.id, sep = "_"))
  
  migratory <- lapply(split(df, df$track_id), function(y){
    # does the bird try to migrate?
    m <- max(y$daily_dist) > d_thresh
    # are the data too gappy for us to judge?
    td <- y %>% 
      mutate(timegap = as.numeric(difftime(timestamp, lag(timestamp), units = "days")))
    td <- max(na.omit(td$timegap)) > g_thresh
    if(m == T & td == F){
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

migration_locations <- lapply(split(migration_locations, migration_locations$track_id), function(x){
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
migration_locations <- lapply(split(migration_locations, migration_locations$track_id), function(x){
  # take the migratory route
  track <- x %>% 
    arrange(timestamp) %>% 
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
  loss <- max(track$timestamp) > loss_time - days(s_thresh)
  full_dist <- unique(track$track_displacement) > l_thresh
  # add a column containing the outcome of the migratory track
  x <- x %>% 
    group_by(track_id) %>% 
    mutate(track_status = ifelse(loss == T & track_id == unique(track$track_id), "incomplete", "complete"),
           track_status = ifelse(full_dist == F, "incomplete", track_status)) %>% 
    ungroup()
  return(x)
}) %>% reduce(rbind)

complete_ml <- migration_locations %>% 
  filter(track_status == "complete")

# saveRDS(migration_locations, file = paste0("C:/Users/hbronnvik/Documents/storkSSFs/migration_locations_40kmCompass_", Sys.Date(), ".rds"))


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


# lapply(split(complete_ml, complete_ml$track_id), function(x){
#   pdf(paste0("C:/Users/hbronnvik/Documents/storkSSFs/m_locs_check/", unique(x$track_id), ".pdf"),
#       width = 8.5, height = 11)
#   print(ggplot(x, aes(location.long, location.lat, color = timestamp)) +
#           borders(xlim = c(-10,10), ylim = c(20, 50))  +
#           geom_point() +
#           theme_classic())
#   dev.off()
# })



