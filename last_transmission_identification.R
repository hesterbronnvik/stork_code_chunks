### The birds that do not contribute data to the study ("dead")
### Hester Brønnvik
### 27.10.2022
### hbronnvik@ab.mpg.de

# required packages
library(move) # access the data
library(moveACC) # manipulate the ACC data
library(lubridate) # manipulate dates
library(stringr) # manipulate strings
library(recurse) # approximate residence times from GPS
library(mapview) # build maps
library(sf) # manipulate spatial data
library(runner) # calculate metrics in rolling windows
library(tidyverse) # ease writing

# required information
setwd("C:/Users/hbronnvik/Documents/stork_code_chunks")
load("C:/Users/hbronnvik/Documents/loginStored.rdata")
studies <- c(24442409, 212096177, 76367850, 21231406, 1176017658, 173641633, 908232414)
d_thresh <- 6*30*24*60*60 # the number of seconds in a day threshold for defining having died
a_thresh <- 0.5 # the average activity a bird has to be below to be dead
r_thresh <- 7 # number of days that the animal is inside the radius to count as dead
g_thresh <- 60 # the number of days that an animal has to be missing to count as having died
t_thresh <- 100 # the number of meters an animal has to move during a gap to count as displaced
c_thresh <- 0.003 # the radius of the circle drawn around each GPS point for residence time estimation

# the files that were created containing the location data of all adults classified 
# as migrating if moving more than m_thresh (number in file name) in a day and given north or south
# down sampled to 15 minute intervals
clean_files <- list.files("C:/Users/hbronnvik/Documents/storkSSFs/clean_data/", pattern = "migration40", full.names = T)

found_birds <- lapply(clean_files, function(x){
  load(x)
  locs$acceleration.raw.x <- NULL
  locs$acceleration.raw.y <- NULL
  locs$acceleration.raw.z <- NULL
  locs$barometric.height <- NULL
  locs$battery.charge.percent <- NULL
  locs$battery.charging.current <- NULL
  locs$external.temperature <- NULL
  locs$gps.hdop <- NULL
  locs$gps.time.to.fix <- NULL
  locs$height.above.msl <- NULL
  locs$light.level <- NULL
  locs$manually.marked.outlier <- NULL
  locs$ornitela.transmission.protocol <- NULL
  locs$tag.voltage <- NULL
  locs$manually.marked.valid <- NULL
  return(locs)
}) %>% reduce(rbind)

# remove the information from non-migration, this should remove all birds that did not migrate
found_birds <- found_birds %>% 
  filter(!str_detect(phase, "no_"))

# save(found_birds, file = paste0("C:/Users/hbronnvik/Documents/storkSSFs/clean_data/reduced40_15m_", Sys.Date(), ".RData"))


## identify individuals that died on migration

# find the last date of the birds' transmissions
# download the ACC reference data
info <- lapply(studies, function(x){
  info <- getMovebankAnimals(x, loginStored) %>% 
    filter(sensor_type_id == 2365683 & individual_id %in% unique(found_birds$individual.id)) %>% 
    select("individual_id","local_identifier", "tag_local_identifier", "timestamp_start", "timestamp_end") %>%
    mutate(pull_start = paste0(gsub("[[:punct:]]| ", "", as.POSIXct(sub(".000", "", timestamp_end), tz = "UTC", origin = "1970-01-01")-d_thresh), "000"),
           study = x)
  return(info)
}) %>% reduce(rbind)

# library(doParallel) 
# cl <- makeCluster(detectCores()-6, type='PSOCK')
# registerDoParallel(cl)
start_time <- Sys.time()
acc_data <- lapply(1:nrow(info), function(x)tryCatch({
  df <- info[x,]
  # download the ACC data for each individual
  acc_df <- getMovebankNonLocationData(study = df$study, animalName = df$individual_id,
                                       sensorID = 2365683, #timestamp_start = info$pull_start[x],
                                       login = loginStored)
  # remove the inconsistent column to allow row binding
  if(!"manually_marked_outlier" %in% colnames(acc_df)){
    acc_df <- acc_df %>% 
      mutate(manually_marked_outlier = NA)
  }
  print(ACCtimeRange(acc_df, units="days"), quotes = F)
  return(acc_df)
}, error = function(msg){ # if an error is thrown (no ACC data), continue and record the ID
  acc_df <- data.frame(matrix(NA,
                              nrow = 1,
                              ncol = 21))
  colnames(acc_df) <- c("individual_id", "deployment_id", "tag_id", "study_id", "sensor_type_id", 
                        "individual_local_identifier", "tag_local_identifier", "individual_taxon_canonical_name", 
                        "data_decoding_software", "eobs_acceleration_axes", "eobs_acceleration_sampling_frequency_per_axis", 
                        "eobs_accelerations_raw", "eobs_key_bin_checksum", "eobs_start_timestamp", 
                        "import_marked_outlier", "timestamp", "event_id",  "visible", "study_name", 
                        "sensor_type", "manually_marked_outlier")
  acc_df$individual_id <- df$individual_id
  return(acc_df)
})) #%>% reduce(rbind) 
# registerDoSEQ()
step1.5_time <- Sys.time() - start_time
missing_acc <- acc_data[sapply(acc_data, nrow) == 1] %>% 
  reduce(rbind) %>% 
  select(individual_id) %>% 
  deframe()
acc_data <- acc_data[sapply(acc_data, nrow) > 1]
# saveRDS(acc_data, file = paste0("C:/Users/hbronnvik/Documents/storkSSFs/acc_data/acc_", length(acc_data), "ids_", Sys.Date(), ".rds"))
# acc_data <- readRDS("C:/Users/hbronnvik/Documents/storkSSFs/acc_data/acc_227ids_2022-12-08.rds") # all the data (35.9 GB)

# lapply(acc_data, function(x){
#   saveRDS(x, file = paste0("C:/Users/hbronnvik/Documents/storkSSFs/acc_data/acc_", unique(x$individual_id), Sys.Date(), ".rds"))
# })
acc_files <- list.files("C:/Users/hbronnvik/Documents/storkSSFs/acc_data/", pattern = "2022-12-12", full.names = T) # each ID's data

# if there is a difference in GPS & ACC end_timestamp, choose
# DBA?
start_time <- Sys.time()
acc_deaths <- lapply(1:length(acc_files), function(x){
  # load the data, reduce size, and remove misreads and pre-tagging/pre-fledging low DBA (within a week of deployment)
  acc <- acc_files[x] %>% 
    readRDS() %>% 
    filter(timestamp > (max(timestamp) - months(12)) & timestamp < Sys.Date() & timestamp > (as.POSIXct(info$timestamp_start[which(info$individual_id == unique(individual_id))], tz = "UTC") + days(7)))

    
  if(nrow(acc) > 1){
    # check the difference between GPS and ACC time ranges
    gps_end <- info %>% 
      filter(individual_id == unique(acc$individual_id)) %>% 
      select(timestamp_end) %>% 
      deframe()
    
    end_diff <- abs(difftime(max(acc$timestamp), as.POSIXct(sub(".000", "", gps_end), tz = "UTC"), units = "days"))

    # if there is a time difference larger than a week, skip it and determine death in a later step
    # unless the ACC has the greater amount of data, then use ACC (presumably more accurate than GPS)
    if(end_diff < 7){
      # calculate VeDBA rowwise
      DBA <- lapply(1:nrow(acc), function(x){
        accRawCol <- grep("accelerations_raw", names(acc), value=T)
        sampFreqCol <- grep("acceleration_sampling_frequency_per_axis", names(acc), value=T)
        axesCol = grep("acceleration_axes", names(acc), value=T)
        accMx <- matrix(as.integer(unlist(strsplit(as.character(acc[x, accRawCol]), " "))), ncol=3, byrow = T)
        n_samples_per_axis <- nrow(accMx)
        acc_burst_duration_s <- n_samples_per_axis/acc[x, sampFreqCol]
        avg_VeDBA <- mean(sqrt((accMx[,1]-mean(accMx[,1]))^2 + (accMx[,2]-mean(accMx[,2]))^2 + (accMx[,3]-mean(accMx[,3]))^2))
        return(avg_VeDBA)
      }) %>% reduce(rbind)
      
      acc <- acc %>% 
        mutate(VeDBA = DBA[,1])
      
      # run a rolling window across the VeDBA calculating standard deviation
      # then select the first time with a deviation of less than 1 and save it
      poi <- acc %>% 
        mutate(rolling_stdev = runner(VeDBA, sd, k = 10)) %>% 
        # discard the first 100 values because the window finds small values here as an artifact
        slice(10:n()) %>% 
        filter(rolling_stdev < a_thresh) %>% 
        select(timestamp, individual_id, individual_local_identifier) %>% 
        arrange(timestamp) %>% 
        slice(1) %>% 
        mutate(comment = "death classified by ACC")
      # ggplot(poi, aes(timestamp, rolling_stdev)) +
      #   geom_point() +
      #   geom_segment(aes(x=timestamp, xend=timestamp, y=0, yend=rolling_stdev)) +
      #   geom_point(data = poi %>% filter(rolling_stdev < a_thresh), color = "red") +
      #   labs(title = unique(acc$individual_local_identifier)) +
      #   theme_classic()
      if(nrow(poi) == 1){print(paste0(unique(acc$individual_local_identifier), " died on ", date(poi$timestamp), "."))
        # save the data individual by individual 
        write.table(poi, file = "C:/Users/hbronnvik/Documents/storkSSFs/acc_deaths_221213.csv", sep = ",",
                    append = TRUE, quote = FALSE,
                    col.names = FALSE, row.names = FALSE)
        return(poi)}else{print(paste0(unique(acc$individual_local_identifier), " did not die ", "."))
          poi <- data.frame(timestamp = NA, individual_id = unique(acc$individual_id), individual_local_identifier = unique(acc$individual_local_identifier), comment = "no death classified by ACC")
          write.table(poi, file = "C:/Users/hbronnvik/Documents/storkSSFs/acc_deaths_221213.csv", sep = ",",
                      append = TRUE, quote = FALSE,
                      col.names = FALSE, row.names = FALSE)
          return(poi)}
    }else{print(paste0(unique(acc$individual_local_identifier), " has a sensor end time disparity of ", round(as.numeric(end_diff), digits = 2), " days."))
    poi <- data.frame(timestamp = NA, individual_id = unique(acc$individual_id), individual_local_identifier = unique(acc$individual_local_identifier), comment = "sensor time disparity")
    write.table(poi, file = "C:/Users/hbronnvik/Documents/storkSSFs/acc_deaths_221213.csv", sep = ",",
                append = TRUE, quote = FALSE,
                col.names = FALSE, row.names = FALSE)
    return(poi)}
  }
}) %>% reduce(rbind) %>% 
  mutate(timestamp = as.POSIXct(timestamp, tz = "UTC", origin = "1970-01-01"))
Sys.time() - start_time
saveRDS(acc_deaths, file = paste0("C:/Users/hbronnvik/Documents/storkSSFs/acc_deaths", Sys.time(), ".rds"))

# see some estimates
  # load the data and reduce size
  acc <- acc_files[2] %>% 
    readRDS() %>% 
    filter(timestamp > (max(timestamp) - months(12)) & timestamp < Sys.Date())
  

    # check the difference between GPS and ACC time ranges
    gps_end <- info %>% 
      filter(individual_id == unique(acc$individual_id)) %>% 
      select(timestamp_end) %>% 
      deframe()
    
    end_diff <- abs(difftime(max(acc$timestamp), as.POSIXct(sub(".000", "", gps_end), tz = "UTC"), units = "days"))
    
    print(paste0(unique(acc$individual_local_identifier), " has a sensor end time disparity of ", round(as.numeric(end_diff), digits = 2), " days."))
    
      acc_sep <- lapply(1:nrow(acc), function(x){
        accRawCol <- grep("accelerations_raw", names(acc), value=T)
        sampFreqCol <- grep("acceleration_sampling_frequency_per_axis", names(acc), value=T)
        axesCol = grep("acceleration_axes", names(acc), value=T)
        accMx <- matrix(as.integer(unlist(strsplit(as.character(acc[x, accRawCol]), " "))), ncol=3, byrow = T)
        n_samples_per_axis <- nrow(accMx)
        acc_burst_duration_s <- n_samples_per_axis/acc[x, sampFreqCol]
        avg_VeDBA <- mean(sqrt((accMx[,1]-mean(accMx[,1]))^2 + (accMx[,2]-mean(accMx[,2]))^2 + (accMx[,3]-mean(accMx[,3]))^2))
        return(avg_VeDBA)
      }) %>% reduce(rbind)
      
      acc <- acc %>% 
        mutate(VeDBA = acc_sep[,1])
      
poi <- acc %>% 
  filter(VeDBA > mean(VeDBA)) %>% 
  arrange(desc(timestamp)) %>% 
  slice(1) %>% 
  mutate(loss = ifelse(timestamp > (max(acc$timestamp) - days(1)), "disappearance", "death"), 
         comment = paste0(loss, " classified by ACC")) %>% 
  select(timestamp, individual_id, individual_local_identifier, comment)
      
ggplot(acc, aes(timestamp, VeDBA)) +
  geom_point() +
  geom_segment(aes(x=timestamp, xend=timestamp, y=0, yend=VeDBA)) +
  geom_point(data = acc %>% filter(VeDBA > mean(VeDBA)), color = "red") +
  geom_vline(xintercept = poi$timestamp, color = "red") +
  labs(title = unique(acc$individual_local_identifier)) +
  theme_classic()


build <- acc %>% 
  filter(timestamp > "2019-08-20") %>% 
  arrange(timestamp) %>% 
  mutate(change = ifelse(lag(VeDBA) < mean(VeDBA) & VeDBA < mean(VeDBA), F, T))

build$change[1] <- F

build <- build %>% 
  mutate(timeLag = as.numeric(timestamp - lag(timestamp), units = "secs"),
       # take the cumulative sum to act as a unique ID for each burst identified in line 37
     cumu_check_for_event = cumsum(change)) %>% 
  group_by(cumu_check_for_event) %>% 
  # for the bursts, calculate the time difference between the last and first locations
  mutate(burstLength = ifelse(n()>1, difftime(tail(timestamp,1), head(timestamp, 1), units = "secs"), NA),
         # add an ID to each burst
         burstID = cur_group_id())

ggplot(build, aes(timestamp, VeDBA)) +
  geom_point() +
  geom_segment(aes(x=timestamp, xend=timestamp, y=0, yend=VeDBA)) +
  geom_point(data = build %>% filter(burstLength == max(burstLength) & VeDBA < mean(VeDBA)), color = "red") +
  geom_vline(xintercept = poi$timestamp, color = "red") +
  labs(title = unique(acc$individual_local_identifier)) +
  theme_classic()


check <- data.frame(index = 1:nrow(build), stdev = runner(build$VeDBA, sd, k = 100))




# the previously estimated death dates from GPS data
gpsd <- readRDS("C:/Users/hbronnvik/Documents/storkSSFs/death_estimates_gps_221206.rds")

### moveACC amplitude examination
# change acc to g using default calibrations
# transfDF <- TransformRawACC(df=acc_df, units="g")
# # FFT
# waveDF <- ACCwave(transfDF, transformedData=T)
# # exploratory plots
# wingBeatsPlot(dfw=waveDF, forclustering= c("amplitude","odbaAvg"))
# ggplot(waveDF, aes(timestamp, amplitude)) +
#   geom_point(alpha = 0.5) +
#   theme_classic()
deaths <- lapply(1:length(acc_files[1:5]), function(x)tryCatch({
  acc_df <- readRDS(acc_files[x])
  # for each individual
  print(unique(acc_df$individual_id))
  # add the tag_local_identifier because moveACC functions use it
  acc_df$tag_local_identifier <- as.integer(info$tag_local_identifier[which(info$individual_id == unique(acc_df$individual_id))])
  # reduce the number of data to save time (and remove future errors)
  acc_df <- acc_df %>% 
    filter(timestamp > (max(acc_df$timestamp) - days(6*30)) & timestamp < Sys.Date())
  # change acc to g using default calibrations
  transfDF <- TransformRawACC(df=acc_df, units="g")
  # FFT
  waveDF <- ACCwave(transfDF, transformedData=T)
  # the first day with low activity
  poi <- waveDF %>% 
    mutate(datestamp = date(timestamp)) %>% 
    mutate(grouping = (amplitude + lag(amplitude) + lead(amplitude))/3) %>% 
    filter(grouping < a_thresh) %>% 
    arrange(datestamp) %>% 
    select(datestamp) %>% 
    slice(1) %>% 
    deframe()
  
  if(length(poi) == 1){
    # within that identified day, the last burst with above average activity
    poi2 <- waveDF %>% 
      filter(date(timestamp) %in% c(poi- days(1), poi)) %>% 
      filter(amplitude > mean(amplitude)) %>% 
      select(timestamp, burstID) %>% 
      slice(n())
    ggplot(waveDF, aes(timestamp, amplitude)) +
      geom_point() +
      geom_segment(aes(x=timestamp, xend=timestamp, y=0, yend=amplitude)) +
      # geom_vline(xintercept = poi2$timestamp, color = "red") +
      theme_classic()
    death <- waveDF %>% 
      filter(burstID %in% poi2$burstID:(poi2$burstID+1)) %>% 
      mutate(dod = mean(timestamp)) %>% 
      select(dod) %>% 
      slice(1) %>% 
      deframe()
    # record the death information
    acc_death <- data.frame(individual_id = unique(acc_df$individual_id), 
                            local_identifier = unique(waveDF$individualID), death_date = death,
                            loss = T, comment = "death classified by ACC amplitude")
  }else{# record the lack of death information
    acc_death <- data.frame(individual_id = unique(acc_df$individual_id), 
                            local_identifier = unique(waveDF$individualID), death_date = NA,
                            loss = F, comment = "no death classified by ACC amplitude")}
  # save the data individual by individual because on 32 GB RAM this crashes
  # write.table(acc_death, file = "acc_deaths_221208.csv", sep = ",",
  #             append = TRUE, quote = FALSE,
  #             col.names = FALSE, row.names = FALSE)
  return(acc_death)
}, error = function(msg){print(geterrmessage())}))


wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")
gps_df <- gps_df %>% 
  mutate(date = as.character(date(timestamp))) %>% 
  filter(timestamp > (tail(timestamp, 1)-days(6*30)))
gps_sf <- st_as_sf(gps_df, coords = c("location.long", "location.lat"), crs = wgs)
mapview(gps_sf, zcol = "date")

# start_time <- Sys.time()
# deaths <- lapply(acc_data, function(x)tryCatch({
#   # for each individual
#   print(unique(x$individual_id))
#   # add the tag_local_identifier on because moveACC functions use it
#   x$tag_local_identifier <- as.integer(info$tag_local_identifier[which(info$individual_id == unique(x$individual_id))])
#   # change acc to g using default calibrations
#   transfDF <- TransformRawACC(df = x, units = "g")
#   # FFT
#   waveDF <- ACCwave(transfDF, transformedData = T)
#   
#   # find the last burst above a given threshold of activity
#   last_active_burst <- waveDF %>% 
#     filter(amplitude > a_thresh) %>%
#     arrange(timestamp) %>% 
#     slice(n()) %>% 
#     select(burstID) %>% 
#     deframe()
#   
#   dd <- waveDF %>% 
#     # take the last active burst and the first burst presumed dead
#     slice(last_active_burst:(last_active_burst + 1)) %>% 
#     # take the average timestamp between these bursts
#     mutate(timestamp = as.POSIXct(timestamp, tz = "UTC", origin = "1970-01-01"),
#            change = difftime(timestamp, lag(timestamp), units = "secs"),
#            midpoint = change/2,
#            last_time = timestamp[1]+midpoint[2]) %>% 
#     # take this midpoint as the last timestamp (it died between last active and first not)
#     select(last_time) %>% 
#     slice(1) %>% 
#     deframe()
#   
#   dates <- data.frame(individual_id = unique(x$individual_id), local_identifier = unique(x$individual_local_identifier), death_date = dd)
#   return(dates)
# }, error = function(msg){print(geterrmessage())})) #%>% reduce(rbind)
# Sys.time() - start_time

acc_deaths <- deaths[lapply(deaths, length) > 1]

acc_deaths <- lapply(acc_deaths, function(x){
  df <- data.frame(individual_id = unique(x$individual_id),
                   local_identifier = unique(x$local_identifier),
                   death_date = unique(x$death_date))
  return(df)
}) %>% reduce(rbind) #%>% 
  # mutate(death_date = as.POSIXct(ifelse(between(date(death_date), runtime -11, runtime) | death_date > Sys.Date(), NA, death_date), tz = "UTC", origin = "1970-01-01"))

saveRDS(acc_deaths, file = paste0("C:/Users/hbronnvik/Documents/storkSSFs/estimated_death_dates_", "2022-12-07", ".rds"))

# left over we have the birds with no ACC downloading or that vanished with normal activity
# use GPS data to estimate those death dates

# the death dates of the birds for which ACC estimates functioned
done <- readRDS("C:/Users/hbronnvik/Documents/storkSSFs/estimated_death_dates_2022-12-07.rds") %>% 
  select(local_identifier) %>% 
  deframe() 

# birds that migrated (50 km segmentation) named found_birds
# load("C:/Users/hbronnvik/Documents/storkSSFs/clean_data/reduced50_15m_2022-10-27.RData")
migratory <- found_birds %>% 
  ungroup() %>% 
  select(individual.local.identifier) %>% 
  rename(local_identifier = individual.local.identifier) %>% 
  deframe() %>% 
  unique()

# the data for the birds in all studies
info <- lapply(studies, function(x){
  i <- getMovebankAnimals(study = x , login = loginStored) %>%
    filter(sensor_type_id == 653) %>% 
    mutate(study = x)
  return(i)
}) %>% reduce(rbind)

# the names of the birds that have not yet been estimated
info <- info %>% 
  filter(!local_identifier %in% done & local_identifier %in% migratory)

# the recurse approach: works well for birds that lie undisturbed at their sites of death until 
# end deployment, but fails on birds that move post-mortem or die suddenly
# thus, added a filter for number of days in one location and large gaps
start_time <- Sys.time()
death_estimates_gps <- lapply(1:nrow(info), function(x){
  print(info$local_identifier[x], quote = F)
  # the data for one bird's last 6 months. This means only 1 of breeding or fledging can be present.
  mv <- getMovebankData(info$study[x], info$local_identifier[x], loginStored, T,
                        timestamp_start = paste0(gsub("[[:punct:]]| ", "",
                                                      (as.POSIXct(sub("\\.000", "", info$timestamp_end[x]), tz = "UTC") - 6*30*24*3600)), "000"))
  # select only locations at a lower than 5 minute resolution (reduce memory burden on recursions)
  mv$lag <- c(NA, move::timeLag(mv, units = "mins"))
  mv <- mv[which(mv$lag > 5),]
  # check for large gaps in transmission
  dd <- data.frame(dates = mv$timestamp) %>% # difference between this and last transmission
    mutate(change = difftime(dates, lag(dates), units = "days"))
  # are the data gappy?
  if(max(na.omit(dd$change)) > g_thresh){
    lost_loc <- data.frame(location_lat = as.numeric(mv$location_lat[timestamps(mv) == dd$dates[which(dd$change > g_thresh)-1]]),
                           location_long = as.numeric(mv$location_long[timestamps(mv) == dd$dates[which(dd$change > g_thresh)-1]]))
    found_loc <- data.frame(location_lat = as.numeric(mv$location_lat[timestamps(mv) == dd$dates[which(dd$change > g_thresh)]]),
                            location_long = as.numeric(mv$location_long[timestamps(mv) == dd$dates[which(dd$change > g_thresh)]]))
    # take the last location before the gap if the tag was moved afterwards
    if(distVincentyEllipsoid(lost_loc, found_loc) > t_thresh) {
      loss <- as.numeric(difftime(tail(mv$timestamp, 1), dd$dates[which(dd$change > g_thresh)-1], units = "days")) > r_thresh
      dod <- data.frame(local_identifier = mv@idData$local_identifier, dod = dd$dates[which(dd$change > g_thresh)-1], loss = loss, comment = "gap")
      # but take the recurse approach if the animal did not move during the gap
    }else{rec <- getRecursions(mv, c_thresh, timeunits = "days")
    # after the animal starts moving (discard time in/on a nest)
    minimum <- rec$revisitStats$entranceTime[which.min(rec$revisitStats$timeInside)]
    # add a threshold for amount of residence time that can be alive
    revists <- rec$revisitStats %>% 
      filter(entranceTime > minimum & timeInside > r_thresh)
    # save the date of entrance to the radius, or if there was none, no evidence of death 
    if(nrow(revists) > 0){
      dod <- revists$entranceTime[which.max(revists$timeInside)]
      loss <- as.numeric(difftime(tail(mv$timestamp, 1), dod, units = "days")) > r_thresh
      dod <- data.frame(local_identifier = mv@idData$local_identifier, dod = dod, loss = loss, comment = "death")
    }else(dod <- data.frame(local_identifier = mv@idData$local_identifier, dod = NA, loss = F, comment = "no death"))}    
    # run the recurse approach on animals that do not have gappy data
  }else{rec <- getRecursions(mv, c_thresh, timeunits = "days")
  minimum <- rec$revisitStats$entranceTime[which.min(rec$revisitStats$timeInside)]
  revists <- rec$revisitStats %>% 
    filter(entranceTime > minimum & timeInside > r_thresh)
  # ggplot(rec$revisitStats, aes(entranceTime, timeInside)) +
  # geom_point() + 
  # geom_segment(aes(x=entranceTime, xend=entranceTime, y=0, yend=timeInside)) +
  # geom_point(data = revists[which.max(revists$timeInside),], color = "red") + 
  # geom_segment(data = revists[which.max(revists$timeInside),], 
  #              aes(x=entranceTime, xend=entranceTime, y=0, yend=timeInside),
  #              color = "red") +
  # theme_classic()
  if(nrow(revists) > 0){
    dod <- revists$entranceTime[which.max(revists$timeInside)]
    loss <- as.numeric(difftime(tail(mv$timestamp, 1), dod, units = "days")) > r_thresh
    dod <- data.frame(local_identifier = mv@idData$local_identifier, dod = dod, loss = loss, comment = "death")
  }else(dod <- data.frame(local_identifier = mv@idData$local_identifier, dod = NA, loss = F, comment = "no death"))}
  return(dod)
  gc()
}) %>% reduce(rbind)
step3_time <- Sys.time() - start_time

full_time <- step1_time + step1.5_time + step2_time + step3_time 

# visually determine the death dates
vis_check <- lapply(1:nrow(info), function(x){
  mv <- getMovebankData(info$study[x], info$local_identifier[x], loginStored, T, 
                        timestamp_start = paste0(gsub("[[:punct:]]| ", "", 
                                                      (as.POSIXct(sub("\\.000", "", info$timestamp_end[x]), tz = "UTC") - 8*30*24*3600)), "000"))
  mv$lag <- c(NA, move::timeLag(mv, units = "mins"))
  mv <- mv[which(mv$lag > 5),]
  mv$date <- as.character(date(mv$timestamp))
  return(mv)
})

# generate a map to double-check the dates
mapview(vis_check[[13]], zcol = "date")

# metadata for corrections to the death dates made on the basis of visual examination of maps (2022-12-06 HB)
dod_corrections <- data.frame(local_identifier = c("Isolde + / DER AU641 (eobs 2760)",
                                                   "Mattis +/ DER AW839 (eobs 3041)",
                                                   "Mira / DER AU650 (eobs 3027)",
                                                   "Borni II + / DER AX351 (eobs 3076)",
                                                   "Nicole + / DER AX555 (eobs 3020)",
                                                   "Maxi2 + / DER AZ919 (e-obs 6586)",
                                                   "Niclas / DER AU053 (eobs 3341)",
                                                   "Tobi + / DER AU052 (eobs 3340)",
                                                   "Redrunner + / DER AU057 (eobs 3339)",
                                                   "Sierit-chick + / DER AZ347 (eobs 4565)",
                                                   "Petra + A9100 (eobs 8016)"), 
                              dod = c("2014-08-28 18:46:08",
                                      "2016-07-07 15:06:54",
                                      NA,
                                      "2018-09-18 18:50:47",
                                      "2018-10-30 09:49:18",
                                      "2020-11-08 08:35:49",
                                      NA,
                                      NA,
                                      "2019-03-27 23:59:58",
                                      "2016-09-12 12:59:03",
                                      "2021-03-21 15:30:08"), 
                              loss = c(T, T, T, T, T, T, T, T, T, T, T), 
                              comment = c("death: misestimation due to rapid end_deployment and long nesting bouts",
                                          "death: misestimation due to residence in the last location for 7 months",
                                          "disappearance: misestimation due to long stay on a Cordoba landfill",
                                          "death: misestimation due to rapid tag recovery",
                                          "death: misestimation due to rapid tag recovery",
                                          "death: misestimation due to 21 day gap",
                                          "disappearance: misestimation due to breeding",
                                          "disappearance: misestimation due to long stay on a Kenitra landfill",
                                          "death: misestimation due to long stay on a Fez mechouar",
                                          "death: misestimation due to long residence in/on the nest",
                                          "death: misestimation due to residence in the last location for a year and GPS reflection"))

# add the corrections to the estimates
death_estimates_gps <- death_estimates_gps %>% 
  filter(!local_identifier %in% dod_corrections) %>% 
  rbind(dod_corrections)

# saveRDS(death_estimates_gps, file = "death_estimates_gps_221206.rds")
gps_deaths <- readRDS("C:/Users/hbronnvik/Documents/storkSSFs/death_estimates_gps_221206.rds") %>% 
  rowwise() %>% 
  mutate(individual_id = info$individual_id[which(info$local_identifier == local_identifier)]) %>% 
  rename(death_date = dod)

last_transmission <- acc_deaths %>% 
  mutate(loss = ifelse(is.na(death_date), NA, T), 
         comment = ifelse(loss == T, "death", "no death")) %>% 
  rbind(gps_deaths)

