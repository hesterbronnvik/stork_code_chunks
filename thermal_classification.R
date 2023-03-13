

### find thermals in stork data
### compare KS/MS method to Weinzierl et al.
source("C:/Users/hbronnvik/Documents/storkSSFs/Laptop/getTrackSegments_updated.R")
source("C:/Users/hbronnvik/Documents/storkSSFs/thermallingFeaturesFunction.R")
source("C:/Users/hbronnvik/Documents/storkSSFs/getWindEstimates_update.R")
library(move)
library(dplyr)
library(tidyr)
library(purrr)
library(lubridate)
library(moveWindSpeed)
load("C:/Users/hbronnvik/Documents/storkSSFs/laptop/loginStored.RData")
source("C:/Users/hbronnvik/Documents/storkSSFs/functions.R")
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

done <- list.files("C:/Users/hbronnvik/Documents/storkBursts/full_burst_data/")
done <- sub(".rds", "", done)
info <- info %>% 
  filter(!animal_id %in% done)

start_time <- Sys.time()
full_data <- lapply(1:nrow(info), function(x){
  df <- getMovebankLocationData(info[x, 2], 653, info[x, 1], loginStored)
  saveRDS(df, file = paste0("C:/Users/hbronnvik/Documents/storkBursts/full_burst_data/", info$animal_id[x], ".rds"))
  return(df)
})
Sys.time() - start_time

full_data_files <- list.files("C:/Users/hbronnvik/Documents/storkBursts/full_burst_data/", full.names = T)
full_data <- lapply(full_data_files, function(x){
  df <- readRDS(x)
  return(df)
})

# remove errors
start_time <- Sys.time()
clean_locations <- lapply(full_data, function(x){
  # load the data
  ind <- x
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
  print(doubles$height.above.ellipsoid)
  
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
        filter(speed_after_unique < 50)
      return(dd)
    }) %>% reduce(rbind)
    locs_df <- locs_df %>% 
      filter(!index %in% check$index)
  }
  
  locs_df <- locs_df %>% 
    filter(ground.speed < 50)
  
  return(locs_df)
})
Sys.time() - start_time

clean_locations <- clean_locations[lapply(clean_locations, nrow) > 1]

# make a list of move objects
clean_locations_mv <- lapply(clean_locations, function(ind){
  print(unique(ind$individual.id))
  mv <- move(x = ind$location.long, y = ind$location.lat, time = ind$timestamp, data = ind, proj = "+proj=longlat +datum=WGS84 +no_defs")
  mv$timelag.sec <- c(NA, timeLag(mv, units="secs"))
  mv$altitude.diff <- c(NA, (mv$height.above.ellipsoid[-1] - mv$height.above.ellipsoid[-nrow(mv)]))
  mv$vert.speed <- mv$altitude.diff/mv$timelag.sec
  mv$turn.angle <- c(NA, turnAngleGc(mv), NA)
  mv$step.length <- c(NA, distance(mv))
  mv$gr.speed <- c(NA, move::speed(mv))
  return(mv)
}) 

classified_bursts <- lapply(clean_locations_mv, function(x){
  ind <- as.data.frame(x)
  ind_bursts <- split_gps_bursts_df(ind, 2, 30)
  if(nrow(ind_bursts) > 0){
    # KS/MS:
    swV <- 2 #smoothing window of 5 seconds (< min burst duration, 2 before 2 after each loc) for vertical speed for later classification
    swT <- 14 #smoothing window of 29 seconds for thermalling behavior (according to Rolf/Bart's paper)
    
    animalID <- unique(ind_bursts$individual.id)
    if(nrow(ind_bursts) > 0){
      # Remove unnecessary columns
      ind_bursts <- ind_bursts[,-grep("mag|orientation|coords|timestamps|start.timestamp|optional|import|visible|algorithm|battery|decoding|accuracy|manually|activity|checksum|acceleration",
                                      colnames(ind_bursts))]
      # Split each individual dataframe by burst ID
      burst_ls_corr <- split(ind_bursts, ind_bursts$burstID)
      # Keep only bursts with at least 40 s (30 of smoothing window will be NA) 
      burst_ls_corr_sub <- burst_ls_corr[which(sapply(burst_ls_corr, nrow) >= 40)]
      # Compute smoothed turning angle separately for each burst
      HRdf <- lapply(burst_ls_corr_sub, function(b){
        b$vertSpeed_smooth <- NA
        b$turnAngle_smooth <- NA
        for(i in (swV+1):(nrow(b)-swV)){
          b$vertSpeed_smooth[i] <- mean(b$vert.speed[(i-swV):(i+swV)], na.rm=T)}
        for(i in (swT+1):(nrow(b)-swT)){
          b$turnAngle_smooth[i] <- max(abs(cumsum(b$turn.angle[(i-swT):(i+swT)])))}
        return(b) # return df with smoothed variables
      }) %>% reduce(rbind)
      # Classify soaring only based on vertical speed
      HRdf <- HRdf[complete.cases(HRdf$vertSpeed_smooth),]
      kmeanV <- kmeans(HRdf$vertSpeed_smooth, 2)   #Get the two clusters
      soarId <- which.max(aggregate(HRdf$vertSpeed_smooth~kmeanV$cluster, FUN=mean)[,2]) # which one is the soaring one?
      soarClust <- rep("glide", length(kmeanV$cluster))
      soarClust[which(kmeanV$cluster==soarId)] <- "soar"
      HRdf$soarClust <- factor(soarClust, levels=c("soar","glide"))  
      # Now classify thermaling only based on turning angle (cumulated to a 30 s time window in previous step)
      HRdf$thermalClust <- "other"
      HRdf$thermalClust[which(HRdf$gr.speed >= 2 & HRdf$turnAngle_smooth >= 300)] <- "circular"
      HRdf$thermalClust[which(HRdf$gr.speed >= 2 & HRdf$soarClust=="soar" & HRdf$thermalClust != "circular")] <- "linear"
      HRdf$thermalClust <- factor(HRdf$thermalClust, levels=c("circular","linear","other"))
    }
    
    # Weinzierl:
    colnames(HRdf)[9] <- "height_above_ellipsoid"
    mv_bursts <- move(x = HRdf$location.long, y = HRdf$location.lat, time = HRdf$timestamp, data = HRdf, proj = "+proj=longlat +datum=WGS84 +no_defs")
    ws <- getWindEstimates(mv_bursts, isThermallingFunction = getDefaultIsThermallingFunction(300),
                           windowSize = 29)
    ws <- as.data.frame(ws)
    ind_thermals <- HRdf %>% 
      full_join(ws[, c(2, 24:40)], by = "timestamp")
    return(ind_thermals)
  }  
})

# library(rgl)
# plot3d(ws$location.long, ws$location.lat, ws$height.above.ellipsoid, 
#        col = ws$col, radius = 2, xlab = "Longitude", ylab = "Latitude", zlab = "HAE", type = "s")



