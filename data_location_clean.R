### cleaned up data location to find the names of birds that should be included (ie. have completed more than one migration)
### Hester Brønnvik
### 29.09.2022


# required packages
library(move)
library(lubridate)
library(tidyverse)

# required information
setwd("C:/Users/hbronnvik/Documents/stork_code_chunks")
load("C:/Users/hbronnvik/Documents/loginStored.rdata")
studies <- c(24442409, 212096177, 76367850, 21231406, 1176017658, 173641633)
m_thresh <- 50*1000 # the threshold above which movement is migratory
s_thresh <- 50 # the threshold above which speed is untrustworthy

# determine the identities of the nestling birds (remove any care center adults)
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

start_time <- Sys.time()
## each study
lost_birds <- lapply(studies, function(x){
  info <- getMovebankAnimals(x, loginStored) %>% 
    filter(sensor_type_id == 653 & individual_id %in% nestlings)
  
  ## each individual
  locs <- lapply(info$individual_id, function(w){
    ind_locs <- getMovebankLocationData(study = x, animalName = w, login = loginStored, sensorID = 653)
    # We have to clean the data before we can determine whether the bird was migrating.
    # If we classify migratory as 100km per day, an outlier could create a false positive.
    # It is best to use only trustworthy fixes to determine life, death, migratory phase, etc.
    
    ## each year
    yr_locs <- lapply(unique(year(ind_locs$timestamp)), function(y){
      
      ## Get the distance covered in each day for the full data (bursts and all)
      locs_df <- ind_locs %>% 
        filter(year(timestamp) == y) %>% 
        drop_na(location.long) %>% 
        mutate(index = row_number())
 
      ## Remove outliers
      if(nrow(locs_df) > 0){
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
              d <- pointDistance(c(poi$location.long, poi$location.lat), c(dd$location.long[p], dd$location.lat[p]), lonlat = T)
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
          d <- pointDistance(c(locs_df$location.long[x], locs_df$location.lat[x]), c(locs_df$location.long[x+1], locs_df$location.lat[x+1]), lonlat = T)
          return(d)
        }) %>% unlist()
        
        tlag <- lapply(1:nrow(locs_df), function(x){
          t <- difftime(locs_df$timestamp[x+1], locs_df$timestamp[x], units = "secs")
          return(t)
        }) %>% unlist()
        
        locs_df <- locs_df %>% 
          mutate(distance = dist,
                 timediff = tlag,
                 gr_speed = distance/tlag) %>%
          # now that the duplicates have been dealt with, remove speed outliers (using instantaneous speed from the tag)
          filter(gr_speed < s_thresh)
        
        ## Now calculate daily distances to classify migration onset
        locs_df <- locs_df %>% 
          group_by(date(timestamp)) %>% 
          mutate(daily_distance = distVincentyEllipsoid(c(head(location.long,1), head(location.lat, 1)), 
                                                        c(tail(location.long,1), tail(location.lat, 1))),
                 # the rhumbline bearing between first and last locations of the day
                 daily_direction = bearingRhumb(c(head(location.long,1), head(location.lat, 1)), 
                                                c(tail(location.long,1), tail(location.lat, 1))),
                 # the global direction of the day to differentiate spring and fall migrations
                 compass_direction = ifelse(daily_direction > 90 & daily_direction < 270, "southward", "northward")) %>% 
          ungroup()
        
        # find the first and last instances of migration in each direction
        locs_df_s <- locs_df %>% 
          filter(daily_distance >= m_thresh & compass_direction == "southward")
        
        on_s <- locs_df_s %>% 
          select(timestamp) %>% 
          slice(1)
        
        off_s <- locs_df_s %>% 
          select(timestamp) %>% 
          slice(n())
        
        locs_df_n <- locs_df %>% 
          filter(daily_distance >= m_thresh & compass_direction == "northward")
        
        on_n <- locs_df_n %>% 
          select(timestamp) %>% 
          slice(1)
        
        off_n <- locs_df_n %>% 
          select(timestamp) %>% 
          slice(n())
        
        # use these to classify migration
        if(nrow(on_s) > 0){
          locs_df_s <- locs_df %>% 
            filter(between(timestamp, on_s, off_s)) %>% 
            mutate(phase = paste0("fall_migration_", y))
        }else{locs_df_s <- locs_df %>% 
          slice(1) %>% 
          mutate(phase = "no_migration")}
        
        if(nrow(on_n) > 0){
          locs_df_n <- locs_df %>% 
            filter(between(timestamp, on_n, off_n)) %>% 
            mutate(phase = paste0("spring_migration_", y))
        }else{locs_df_n <- locs_df %>% 
          slice(1) %>% 
          mutate(phase = "no_spring_migration")}
        
        ## Finally, sub-sample to 15 minute intervals to remove burst data and reduce memory use
        locs_df <- rbind(locs_df_s, locs_df_n)  %>% 
          arrange(timestamp) %>% 
          mutate(seq15 = round_date(timestamp, "15 minutes")) %>% 
          group_by(seq15) %>% 
          slice(1)
      }else{locs_df <- ind_locs %>% 
        mutate(phase = "no_locations")}
      
      print(paste0("Completed the classification of individual ", w, " in ", y, "."), quote = F)
      return(locs_df) # finish the year
    }) %>% reduce(rbind)
    
    return(yr_locs) # finish the ID
  }) %>% reduce(rbind)
  
  print(paste0("Completed the classification of study ", x, "."), quote = F)
  # save(locs, file = paste0("C:/Users/hbronnvik/Documents/storkSSFs/clean_data/migration50_15m_", x, "_", Sys.Date(), ".RData"))
  return(locs) # finish the study
}) %>% reduce(rbind)
Sys.time() - start_time

# the files that were just created containing the location data of all adults classified 
# as migrating if moving more than m_thresh in a day and given north or south
# down sampled to 15 minute intervals
clean_files <- list.files("C:/Users/hbronnvik/Documents/storkSSFs/clean_data/", pattern = "migration", full.names = T)

found_birds <- lapply(clean_files, function(x){
  load(x)
  return(locs)
}) %>% reduce(rbind)

# remove the information from non-migration, this should remove all birds that did not migrate
found_birds <- found_birds %>% 
  filter(!str_detect(phase, "no_"))

save(found_birds, paste0("C:/Users/hbronnvik/Documents/storkSSFs/clean_data/reduced50_15m_", Sys.Date(), ".RData"))

# remove birds for which stopover cannot be differentiated from wintering
found_birds <- found_birds %>% 
  filter(year(timestamp) != year(Sys.Date()))

## identify individuals that died on migration

# find the last date of the birds' transmissions
last_known <- lapply(studies, function(x){
  info <- getMovebankAnimals(x, loginStored) %>% 
    filter(sensor_type_id == 653 & individual_id %in% nestlings)
})
# for birds that have death comments containing a date, extract that information
info[grep("[0-9]{2}[[:punct:]][0-9]{2}[[:punct:]][0-9]{4}",info$death_comments),]
# remove those birds from the search

# for remaining birds, use the last date
# find a start date a threshold of days before that
# download the ACC 
# calculate DBA
# if under a threshold, call this animal dead
# save the death date

# compare the death date to the migration
# if within a certain time of the migration, it died on migration
# remove the bird

# is the last time the animal moved 50km per day within the last week it transmitted?

# unfinished <- lapply(studies, function(x){
#   info <- getMovebankAnimals(x, loginStored) %>% 
#     filter(sensor_type_id == 653 & individual_id %in% nestlings)
#   
#   current <- info %>%
#     filter(year(timestamp_start) == 2022 & timestamp_end > Sys.Date()-2)
#   
#   return(current)
# }) %>% reduce(rbind)
# 
# # take out the birds from this year that have not yet finished their first migration
# nestlings <- nestlings[!nestlings %in% unfinished$individual_id]


# Duplicates containing HAE and DOP values exist.
# Completed the classification of individual 1576780467 in 2021.
# Duplicates containing HAE and DOP values exist.
# Completed the classification of individual 1578616761 in 2022.

