### data location to find the names of birds that should be included (ie. have completed more than one migration)
### Hester Brønnvik
### 29.09.2022


# required packages
library(move)
library(lubridate)
library(stringr)
library(EMbC)
library(tidyverse)

# required information
setwd("C:/Users/hbronnvik/Documents/stork_code_chunks")
load("C:/Users/hbronnvik/Documents/loginStored.rdata")
studies <- c(24442409, 212096177, 76367850, 21231406, 1176017658, 173641633)#, 908232414)
m_thresh <- 100*1000 # the threshold above which movement is migratory (m/day)
s_thresh <- 50 # the threshold above which speed is untrustworthy (m/s)

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

info <- lapply(studies, function(x){
  info <- getMovebankAnimals(x, loginStored) %>% 
    # the gps records of animals tagged as juveniles
    filter(sensor_type_id == 653 & individual_id %in% nestlings) %>% 
    mutate(study = x)
  return(info)
}) %>% reduce(rbind)

# done <- str_sub(list.files("C:/Users/hbronnvik/Documents/storkSSFs/full_data/"), 1, -16)
# info <- info %>%
#   filter(!individual_id %in% done) %>% 
#   arrange(individual_id)

# get the data
gps <- lapply(1:nrow(info), function(x){
  df <- getMovebankLocationData(study = info$study[x], animalName = info$individual_id[x], login = loginStored)
  saveRDS(df, file = paste0("C:/Users/hbronnvik/Documents/storkSSFs/full_data/", info$individual_id[x], "_", Sys.Date()-2, ".rds"))
  print(paste0("Downloaded individual: ", info$individual_id[x], ", ", x, " of ", nrow(info), "."))
  return(df)
})

# call the data
gps_data <- list.files("C:/Users/hbronnvik/Documents/storkSSFs/full_data/", full.names = T)

start_time <- Sys.time()
## each study
clean_gps_data <- lapply(30, function(x){
    ind_locs <- readRDS(gps_data[x])
    # We have to clean the data before we can determine whether the bird was migrating.
    # If we classify migratory as based on movement parameters, an outlier could create a false positive.
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
          filter(gr_speed < s_thresh) %>% 
          group_by(date(timestamp)) %>%
          # the daily displacement (m) as a way to gauge whether individuals migrate at all
          mutate(daily_distance = distVincentyEllipsoid(c(head(location.long,1), head(location.lat, 1)), 
                                                        c(tail(location.long,1), tail(location.lat, 1)))) %>% 
          ungroup() 
        
        ## Now cluster based on speed and turning angle to classify migration (HL)
        if(max(locs_df$daily_distance) > m_thresh){
          # if it did, run EMbC to identify migratory times
          loco <- locs_df %>%
            # downsample
            mutate(aggregator = round_date(timestamp, "hour")) %>%
            group_by(aggregator) %>%
            slice(1) %>%
            ungroup() %>%
            # prep for EMbC
            select(timestamp, location.long, location.lat) %>%
            as.data.frame()
          # run EMbC
          bc <- stbc(loco, smth = 24)
          # EMbC::view(bc)
          loco <- loco %>% 
            # add on the EMbC clusters
            mutate(category_embc = bc@A,
                   year = year(timestamp))
          
          
        # find the first and last instances of migration in each season
         if(unique(loco$year) == unique(year(ind_locs$timestamp))[1]){
           # separate the first year, when fall is the first migration
              s1 <- loco %>% 
                filter(category_embc == 3) %>% 
                slice(1) %>% 
                select(timestamp) %>% 
                deframe()
              e1 <- loco %>% 
                filter(category_embc == 3) %>% 
                slice(n()) %>% 
                select(timestamp) %>% 
                deframe()
              s_df <- data.frame(individual.id = unique(locs_df$individual.id), local.identifier = unique(locs_df$individual.local.identifier),
                                 start_time = s1, end_time = e1, season = "post-breeding")
            }else{
              loco <- loco %>% 
                group_by(date(timestamp))%>%
                mutate(daily_distance = distVincentyEllipsoid(c(head(location.long,1), head(location.lat, 1)), 
                                                              c(tail(location.long,1), tail(location.lat, 1)))) %>% 
              ungroup() 
              # an arbitrary point identifying the breeding season
              nm <- loco %>%
                filter(month(timestamp) %in% 6:8) %>% 
                filter(daily_distance == min(daily_distance)) %>% 
                slice(n()) %>% 
                select(timestamp) %>% 
                deframe()
              if(length(nm) > 0){ # if the animal did survive to summer:
              # start of the pre-breeding
              s1 <- loco %>% 
                filter(category_embc == 3) %>% 
                select(timestamp) %>% 
                slice(1) %>% 
                deframe()
              # end of the pre-breeding
              e1 <- loco %>% 
                filter(timestamp < nm & category_embc == 3) %>% 
                select(timestamp) %>% 
                slice(n()) %>% 
                deframe()
              # start of post-breeding
              s2 <- loco %>% 
                filter(category_embc == 3 & timestamp > nm) %>% 
                select(timestamp) %>% 
                slice(1) %>% 
                deframe()
              # end of post-breeding
              e2 <- loco %>% 
                filter(category_embc == 3 & timestamp > nm) %>% 
                select(timestamp) %>% 
                slice(n()) %>% 
                deframe()
              s_df <- data.frame(individual.id = unique(locs_df$individual.id), local.identifier = unique(locs_df$individual.local.identifier),
                                 start_time = c(s1, s2), end_time = c(e1, e2), season = c("pre-breeding", "post-breeding"))
              }else{# if the animal did not survive to summer:
              # start of the pre-breeding
              s1 <- loco %>% 
                filter(category_embc == 3) %>% 
                slice(1) %>% 
                select(timestamp) %>% 
                deframe()
              e1 <- loco %>% 
                filter(category_embc == 3) %>% 
                slice(n()) %>% 
                select(timestamp) %>% 
                deframe()
              s_df <- data.frame(individual.id = unique(locs_df$individual.id), local.identifier = unique(locs_df$individual.local.identifier),
                                 start_time = s1, end_time = e1, season = "pre-breeding")
            }
            }
          # use these to classify migration
          locs_df <- lapply(1:nrow(s_df), function(x){
            s <- locs_df %>% 
              filter(between(timestamp, s_df$start_time[x], s_df$end_time[x])) %>% 
              mutate(track_id = paste(unique(individual.id), unique(year(s_df$start_time[x])), s_df$season[x], sep = "_"))
          }) %>% reduce(rbind)
          
          ## Finally, sub-sample to 15 minute intervals to remove burst data and reduce memory
          locs_df <- locs_df %>% 
            arrange(timestamp) %>% 
            mutate(seq15 = round_date(timestamp, "15 minutes")) %>% 
            group_by(seq15) %>% 
            slice(1) %>% 
            ungroup()
          
          print(paste0("Completed the classification of individual ", x, " in ", y, "."), quote = F)
          return(locs_df) # finish the year
        }else{
          locs_df <- locs_df %>% 
            arrange(timestamp) %>% 
            mutate(seq15 = round_date(timestamp, "15 minutes")) %>% 
            group_by(seq15) %>% 
            slice(1) %>% 
            ungroup() %>% 
            mutate(track_id = paste(unique(individual.id), unique(year(locs_df$timestamp)), "no_migration", sep = "_"))
          
          print(paste0("Completed the classification of individual ", x, " in ", y, "."), quote = F)
          return(locs_df) # finish the year
        }
      }
    }) %>% reduce(rbind)
    
  print(paste0("Completed the classification of individual ", x, " of ", length(gps_data), "."), quote = F)
  saveRDS(yr_locs, file = paste0("C:/Users/hbronnvik/Documents/storkSSFs/clean_data/migration_embc_15min_", x, "_", Sys.Date(), ".rds"))
  return(yr_locs) # finish the ID
}) %>% reduce(rbind)
Sys.time() - start_time

# check <- readRDS("C:/Users/hbronnvik/Documents/storkSSFs/clean_data/migration_embc_15min_30_2023-01-30.rds")
# 
# 
# data_sp <- check
# coordinates(data_sp) <- ~ location.long + location.lat # set coordinates
# proj4string(data_sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs") # set projection
# mapView(data_sp, zcol = "track_id", burst = F, cex = 3, color = rainbow) # plot on a map
