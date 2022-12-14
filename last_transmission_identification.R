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
    # select("individual_id","local_identifier", "tag_local_identifier", "timestamp_start", "timestamp_end") %>%
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
Sys.time() - start_time
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
acc_m <- read.csv("C:/Users/hbronnvik/Documents/storkSSFs/acc_deaths_221214.csv", header = T, row.names = F) %>% 
  drop_na(V2)

acc_names <- gsub("2022-12-12.rds", "", gsub("C:/Users/hbronnvik/Documents/storkSSFs/acc_data/acc_", "", acc_files))
done <- acc_m$V2
todo <- acc_names[!acc_names %in% done]
todo <- sapply(todo, function(x){grep(x, acc_files)})

library(parallel)
start_time <- Sys.time()
acc_deaths <- lapply(1:length(acc_files), function(x){
  # load the data and remove misreads and pre-tagging/pre-fledging low DBA (within a week of deployment)
  acc <- acc_files[x] %>% 
    readRDS() %>% 
    filter(timestamp < Sys.Date() & timestamp > (as.POSIXct(info$timestamp_start[which(info$individual_id == unique(individual_id))][1], tz = "UTC") + days(7))) %>% 
    dplyr::mutate(index = 1:n())
  # filter the data to reduce memory burden
  indexer <- (nrow(acc) - 35000):nrow(acc)
  acc <- acc %>% 
    filter(index %in% indexer)
    
  if(nrow(acc) > 1){
    # check the difference between GPS and ACC time ranges
    gps_end <- info %>% 
      filter(individual_id == unique(acc$individual_id)) %>% 
      dplyr::select(timestamp_end) %>% 
      slice(1) %>% 
      deframe()
    
    end_diff <- abs(difftime(max(acc$timestamp), as.POSIXct(sub(".000", "", gps_end), tz = "UTC"), units = "days"))
    # if there is a time difference larger than a week, skip it and determine death in a later step
    if(end_diff < 7){
      # set tag sensitivity
      if(as.numeric(unique(acc$tag_local_identifier))[1] < 2241){sensitive <- data.frame(TagID = as.numeric(unique(acc$tag_local_identifier)), sensitivity = "low")}
      # transform the raw ACC data into g
      acc <- TransformRawACC(df = acc, sensitivity.settings = sensitive, units = "g")
      # extract each axis' values row-wise and use a rolling window to calculate sd of each 10 bursts
      cl <- makeCluster(detectCores()-6)
      clusterExport(cl, "acc", envir = environment())
      acc_sep <- parLapply(cl, 1:nrow(acc), function(x){
        library(runner)
        accTCol <- grep("accelerationTransformed", names(acc), value=T)
        sampFreqCol <- grep("acceleration_sampling_frequency_per_axis", names(acc), value=T)
        axesCol = grep("acceleration_axes", names(acc), value=T)
        accMx <- matrix(as.numeric(unlist(strsplit(as.character(acc[x, accTCol]), " "))), ncol=3, byrow = T)
        X <- mean(na.omit(runner(accMx[,1], sd, k = 10)))
        Y <- mean(na.omit(runner(accMx[,2], sd, k = 10)))
        Z <- mean(na.omit(runner(accMx[,3], sd, k = 10)))
        V <- data.frame(X = X, Y = Y, Z = Z, timestamp = acc$timestamp[x])
        return(V)
      }) %>% reduce(rbind)
      stopCluster(cl)
      
      # select the first instance of low ACC trace sd and save
      poi <- acc_sep %>% 
        # remove night hours when the birds have low activity regardless of death
        filter(hour(timestamp) %in% 8:20) %>% 
        # find the average daytime activity for the bird
        dplyr::mutate(avgX = mean(X),
               avgY = mean(Y),
               avgZ = mean(Z)) %>% 
        group_by(date(timestamp)) %>% 
        # find the average daytime activity for the day
        dplyr::mutate(avg_dayX = mean(X),
               avg_dayY = mean(Y),
               avg_dayZ = mean(Z)) %>% 
        ungroup() %>% 
        dplyr::mutate(criterion = ifelse(X < 0.001 & lag(X) < 0.001 & Y < 0.0015 & lag(Y) < 0.0015 & Z < 0.001 & lag(Z) < 0.001, T, F)) %>% 
        # select the first time the animal was below average activity all day
        filter(avg_dayX < avgX & avg_dayY < avgY & avg_dayZ < avgZ) %>% 
        # and had almost no change in activity over a 10 burst time window
        filter(criterion == T) %>% 
        dplyr::arrange(timestamp) %>% 
        slice(1) 
      
      if(nrow(poi == 1)){
        # using that date, take the time the animal had almost no change in activity
        poi2 <- acc_sep %>% 
          filter(date(timestamp) %in% c(date(poi$timestamp), (date(poi$timestamp) - days(1)))) %>% 
          filter(X < 0.001 & Y < 0.0015 & Z < 0.001) %>% 
          dplyr::arrange(timestamp) %>% 
          slice(1)%>% 
          dplyr::mutate(individual_id = unique(acc$individual_id),
                 individual_local_identifier = unique(acc$individual_local_identifier),
                 loss = T, 
                 comment = "classified by ACC") %>% 
          dplyr::select(timestamp, individual_id, individual_local_identifier, loss, comment)
      }else{
        poi2 <- acc_sep %>% 
          filter(date(timestamp) %in% c(max(date(timestamp)), (max(date(timestamp)) - days(1)))) %>% 
          filter(X < 0.001 & Y < 0.0015 & Z < 0.001) %>% 
          dplyr::arrange(timestamp) %>% 
          slice(1)%>% 
          dplyr::mutate(individual_id = unique(acc$individual_id),
                 individual_local_identifier = unique(acc$individual_local_identifier),
                 loss = T, 
                 comment = "classified by ACC") %>% 
          dplyr::select(timestamp, individual_id, individual_local_identifier, loss, comment)}
      
      if(nrow(poi2) == 1){print(paste0(unique(acc$individual_local_identifier), " died on ", date(poi2$timestamp), "."))
        # save the data individual by individual 
        write.table(poi2, file = "C:/Users/hbronnvik/Documents/storkSSFs/acc_deaths_221214.csv", sep = ",",
                    append = TRUE, quote = FALSE,
                    col.names = FALSE, row.names = FALSE)
        return(poi2)}else{print(paste0(unique(acc$individual_local_identifier), " did not die", "."))
          poi <- data.frame(timestamp = NA, individual_id = unique(acc$individual_id), individual_local_identifier = unique(acc$individual_local_identifier), loss = F, comment = "no death classified by ACC")
          write.table(poi, file = "C:/Users/hbronnvik/Documents/storkSSFs/acc_deaths_221214.csv", sep = ",",
                      append = TRUE, quote = FALSE,
                      col.names = FALSE, row.names = FALSE)
          return(poi)}
    }else{print(paste0(unique(acc$individual_local_identifier), " has a sensor end time disparity of ", round(as.numeric(end_diff), digits = 2), " days."))
    poi <- data.frame(timestamp = NA, individual_id = unique(acc$individual_id), individual_local_identifier = unique(acc$individual_local_identifier), loss = NA, comment = "sensor time disparity")
    write.table(poi, file = "C:/Users/hbronnvik/Documents/storkSSFs/acc_deaths_221214.csv", sep = ",",
                append = TRUE, quote = FALSE,
                col.names = FALSE, row.names = FALSE)
    return(poi)}
  }
}) %>% reduce(rbind) %>% 
  dplyr::mutate(timestamp = as.POSIXct(timestamp, tz = "UTC", origin = "1970-01-01"))
Sys.time() - start_time
saveRDS(acc_deaths, file = paste0("C:/Users/hbronnvik/Documents/storkSSFs/acc_deaths", Sys.time(), ".rds"))

# see some estimates
  # load the data and reduce size
  acc <- acc_files[todo[5]] %>% 
    readRDS() %>% 
    filter(timestamp < Sys.Date()) %>% 
    filter(timestamp > (max(timestamp) - months(12))) %>% 
    mutate(index = 1:n())
  # filter the data to reduce memory burden
  indexer <- (nrow(acc) - 35000):nrow(acc)
  acc <- acc %>% 
    filter(index %in% indexer)
  
  # check the difference between GPS and ACC time ranges
  gps_end <- info %>% 
      filter(individual_id == unique(acc$individual_id)) %>% 
      select(timestamp_end) %>% 
      deframe()
    
    end_diff <- abs(difftime(max(acc$timestamp), as.POSIXct(sub(".000", "", gps_end), tz = "UTC"), units = "days"))
    
    print(paste0(unique(acc$individual_local_identifier), " has a sensor end time disparity of ", round(as.numeric(end_diff), digits = 2), " days."))
    
    
    # calculate VeDBA rowwise
    # DBA <- lapply(1:nrow(acc), function(x){
    #   accRawCol <- grep("accelerations_raw", names(acc), value=T)
    #   sampFreqCol <- grep("acceleration_sampling_frequency_per_axis", names(acc), value=T)
    #   axesCol = grep("acceleration_axes", names(acc), value=T)
    #   accMx <- matrix(as.integer(unlist(strsplit(as.character(acc[x, accRawCol]), " "))), ncol=3, byrow = T)
    #   n_samples_per_axis <- nrow(accMx)
    #   acc_burst_duration_s <- n_samples_per_axis/acc[x, sampFreqCol]
    #   avg_VeDBA <- mean(sqrt((accMx[,1]-mean(accMx[,1]))^2 + (accMx[,2]-mean(accMx[,2]))^2 + (accMx[,3]-mean(accMx[,3]))^2))
    #   return(avg_VeDBA)
    # }) %>% reduce(rbind)
    # 
    # acc <- acc %>% 
    #   mutate(VeDBA = DBA[,1])
    # transform the raw ACC data into g
    acc <- TransformRawACC(df = acc, sensitivity.settings = sensitive, units = "g")
    # extract each axis' values rowwise and use a rolling window to calculate sd of each 10 bursts
    acc_sep <- lapply(1:nrow(acc), function(x){
      accTCol <- grep("accelerationTransformed", names(acc), value=T)
      sampFreqCol <- grep("acceleration_sampling_frequency_per_axis", names(acc), value=T)
      axesCol = grep("acceleration_axes", names(acc), value=T)
      accMx <- matrix(as.numeric(unlist(strsplit(as.character(acc[x, accTCol]), " "))), ncol=3, byrow = T)
      X <- mean(na.omit(runner(accMx[,1], sd, k = 10)))
      Y <- mean(na.omit(runner(accMx[,2], sd, k = 10)))
      Z <- mean(na.omit(runner(accMx[,3], sd, k = 10)))
      V <- data.frame(X = X, Y = Y, Z = Z, timestamp = acc$timestamp[x])
      return(V)
    }) %>% reduce(rbind)
    
# select the first instance of low ACC trace sd and save
poi <- acc_sep %>% 
  # remove night hours when the birds have low activity regardless of death
  filter(hour(timestamp) %in% 8:20) %>% 
  # find the average daytime activity for the bird
  mutate(avgX = mean(X),
         avgY = mean(Y),
         avgZ = mean(Z)) %>% 
  group_by(date(timestamp)) %>% 
  # find the average daytime activity for the day
  mutate(avg_dayX = mean(X),
         avg_dayY = mean(Y),
         avg_dayZ = mean(Z)) %>% 
  ungroup() %>% 
  mutate(criterion = ifelse(X < 0.001 & lag(X) < 0.001 & Y < 0.0015 & lag(Y) < 0.0015 & Z < 0.001 & lag(Z) < 0.001, T, F)) %>% 
  # select the first time the animal was below average activity all day
  filter(avg_dayX < avgX & avg_dayY < avgY & avg_dayZ < avgZ) %>% 
  # and had almost no change in activity over a 10 burst time window
  filter(criterion == T) %>% 
  arrange(timestamp) %>% 
  slice(1) 
# using that date, take the time the animal had almost no change in activity
poi2 <- acc_sep %>% 
  filter(date(timestamp) %in% c(date(poi$timestamp), (date(poi$timestamp) - days(1)))) %>% 
  filter(X < 0.001 & Y < 0.0015 & Z < 0.001) %>% 
  arrange(timestamp) %>% 
  slice(1) %>% 
      mutate(individual_id = unique(acc$individual_id),
             individual_local_identifier = unique(acc$individual_local_identifier),
             loss = T, 
             comment = "classified by ACC") %>% 
      select(timestamp, individual_id, individual_local_identifier, loss, comment)
    
build <- acc_sep %>% 
  filter(timestamp > "2022-02-13" & timestamp < "2022-02-15") %>% 
  group_by(date(timestamp)) %>% 
  # find the average daytime activity for the day
  mutate(avg_dayX = mean(X),
         avg_dayY = mean(Y),
         avg_dayZ = mean(Z)) %>% 
  ungroup()
ggplot(build, aes(timestamp, Z)) +
  geom_line(color = "orange") +
  geom_line(data = build, aes(timestamp, Y), color = "red") +
  geom_line(data = build, aes(timestamp, X), color = "blue") +
  geom_point(data = build %>% filter(X < 0.001 & Y < 0.0015 & Z < 0.001), color = "black") +
  geom_hline(yintercept = c(poi$avgX, poi$avgY, poi$avgZ)) +
  geom_line(data = build, aes(timestamp, avg_dayX), lty = 2) +
  # geom_vline(xintercept = poi$timestamp, lty = 2) +
  labs(title = unique(acc$individual_local_identifier), y = "standard deviation (g)") +
  theme_classic()
    
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

ggplot(acc, aes(timestamp, VeDBA)) +
  geom_point() +
  geom_segment(aes(x=timestamp, xend=timestamp, y=0, yend=VeDBA)) +
  # geom_point(data = acc %>% filter(VeDBA > mean(VeDBA)), color = "red") +
  geom_vline(xintercept = poi$timestamp, color = "red") +
  labs(title = unique(acc$individual_local_identifier)) +
  theme_classic()


# the previously estimated death dates from GPS data
gpsd <- readRDS("C:/Users/hbronnvik/Documents/storkSSFs/death_estimates_gps_221206.rds")

wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")
gps_df <- gps_df %>% 
  mutate(date = as.character(date(timestamp))) %>% 
  filter(timestamp > (tail(timestamp, 1)-days(6*30)))
gps_sf <- st_as_sf(gps_df, coords = c("location.long", "location.lat"), crs = wgs)
mapview(gps_sf, zcol = "date")

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

