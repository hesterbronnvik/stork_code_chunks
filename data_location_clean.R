### cleaned up data location to find the names of birds that should be included (ie. have completed more than one migration)
### Hester Brønnvik
### 29.09.2022


# required packages
library(move)
library(lubridate)
library(tidyverse)

# required information
setwd("C:/Users/hbronnvik/Documents/stork_code_chunks")
load("loginStored.rdata")
studies <- c(24442409, 212096177, 76367850, 21231406, 1176017658, 173641633)

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

# did the animal ever transmit  100 km per day?
# is the last time the animal moved 100km per day within the last week it transmitted?

# the birds that are not yet done with their 2022 migration
unfinished <- lapply(studies, function(x){
  info <- getMovebankAnimals(x, loginStored) %>% 
    filter(sensor_type_id == 653 & individual_id %in% nestlings)
  
  current <- info %>%
    filter(year(timestamp_start) == 2022 & timestamp_end > Sys.Date()-2)
  
  return(current)
}) %>% reduce(rbind)

# take out the birds from this year that have not yet finished their first migration
nestlings <- nestlings[!nestlings %in% unfinished$individual_id]

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
        filter(year(timestamp) == y) 
      
      ## Remove outliers
      # calculate ground speeds
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
               gr_speed = distance/timediff)
      
      # remove duplicated locations in the same time
      duplicate_spacetime <- locs_df %>%
        dplyr::select(location.long, location.lat, individual.id, timestamp) %>%
        duplicated()
      if(sum(duplicate_spacetime) > 0){locs_df <- locs_df[!duplicate_spacetime,]}
      
      # if timestamps are duplicated in different space, remove the location with the higher move::speed
      duplicate_time <- locs_df %>%
        dplyr::select(timestamp) %>%
        duplicated()
      if(sum(duplicate_time) > 0){dt <- locs_df[which(locs_df$timestamp %in% locs_df[duplicated(locs_df$timestamp),"timestamp"]),];
      dt <- rownames_to_column(dt);
      locs_df <- locs_df[-as.numeric(dt$rowname[which.max(dt$gr_speed)]),]}
      
      # now that the duplicates have been dealt with, remove speed outliers (using instantaneous speed from the tag)
      locs_df <- locs_df %>% 
        filter(ground.speed < 50)
      
      # calculate ground speeds again now that faulty locations are removed
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
               gr_speed = distance/tlag)
      
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
      
      return(locs_df)
    })
    
    return(yr_locs)
  })
})













