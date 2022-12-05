### The birds that do not contribute data to the study ("dead")
### Hester Brønnvik
### 27.10.2022
### hbronnvik@ab.mpg.de


# required packages
library(move)
library(moveACC)
library(lubridate)
library(stringr)
library(tidyverse)

# required information
setwd("C:/Users/hbronnvik/Documents/stork_code_chunks")
load("C:/Users/hbronnvik/Documents/loginStored.rdata")
studies <- c(24442409, 212096177, 76367850, 21231406, 1176017658, 173641633)
d_thresh <- 120*24*60*60 # the number of seconds in a day threshold for defining having died
a_thresh <- 0.25 # the average activity a bird has to be below to be dead

# the files that were created containing the location data of all adults classified 
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

# save(found_birds, file = paste0("C:/Users/hbronnvik/Documents/storkSSFs/clean_data/reduced50_15m_", Sys.Date(), ".RData"))


## identify individuals that died on migration

# find the last date of the birds' transmissions
# download the ACC reference data
info <- lapply(studies, function(x){
  info <- getMovebankAnimals(x, loginStored) %>% 
    filter(sensor_type_id == 2365683 & individual_id %in% unique(found_birds$individual.id)) %>% 
    select("individual_id", "timestamp_start", "timestamp_end") %>% 
    mutate(pull_start = paste0(gsub("[[:punct:]]| ", "", as.POSIXct(sub(".000", "", timestamp_end), tz = "UTC", origin = "1970-01-01")-d_thresh), "000"),
           study = x)
  
}) %>% reduce(rbind)

acc_data <- lapply(1:nrow(info), function(x)tryCatch({
  df <- info[x,]
  # download the ACC data for each individual
  acc_df <- getMovebankNonLocationData(study = df$study, animalName = 173673531,
                                       sensorID=2365683, timestamp_start = info$pull_start[x],
                                       login=loginStored)
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

# save(acc_data, file = "C:/Users/hbronnvik/Documents/storkSSFs/acc_data/acc_309ids_11.11.2022.RData")
load("C:/Users/hbronnvik/Documents/storkSSFs/acc_data/acc_94ids_31.10.2022.RData")


### moveACC amplitude examination
# change acc to g using default calibrations
# transfDF <- TransformRawACC(df=acc_df, units="g")
# # FFT
# waveDF <- ACCwave(transfDF, transformedData=T)
# # exploratory plots
# wingBeatsPlot(dfw=waveDF, forclustering= c("amplitude","odbaAvg"))
# plot(waveDF$timestamp, waveDF$amplitude, xlab = "Timestamp", ylab = "Burst amplitude")
# ggplot(waveDF, aes(timestamp, amplitude)) +
#   geom_point(alpha = 0.5) +
#   theme_classic()


deaths <- lapply(acc_data, function(x)tryCatch({
  # for each individual
  print(unique(x$individual_id))
  # change acc to g using default calibrations
  transfDF <- TransformRawACC(df = x, units = "g")
  # FFT
  waveDF <- ACCwave(transfDF, transformedData = T)
  
  # find the last burst above a given threshold of activity
  last_active_burst <- waveDF %>% 
    filter(amplitude > a_thresh) %>%
    arrange(timestamp) %>% 
    slice(n()) %>% 
    select(burstID) %>% 
    deframe()
  
  dd <- waveDF %>% 
    # take the last active burst and the first burst presumed dead
    slice(last_active_burst:(last_active_burst + 1)) %>% 
    # take the average timestamp between these bursts
    mutate(timestamp = as.POSIXct(timestamp, tz = "UTC", origin = "1970-01-01"),
           change = difftime(timestamp, lag(timestamp), units = "secs"),
           midpoint = change/2,
           last_time = timestamp[1]+midpoint[2]) %>% 
    # take this midpoint as the last timestamp (it died between last active and first not)
    select(last_time) %>% 
    slice(1) %>% 
    deframe()
  
  dates <- data.frame(individual_id = unique(x$individual_id), local_identifier = unique(x$individual_local_identifier), death_date = dd)
  return(dates)
}, error = function(msg){print(geterrmessage())})) #%>% reduce(rbind)
# 219402564
# "Tag numbers not recognized. Make sure the class of column 'tag.local.identifier' or 'tag_local_identifier' is 'integer'"
# 89349490
# [1] "<text>:1:7: unexpected '{'\n1: .code {\n          ^"

acc_deaths <- deaths[lapply(deaths, length) > 1]

acc_deaths <- lapply(acc_deaths, function(x){
  df <- data.frame(individual_id = unique(x$individual_id.individual_id),
                   local_identifier = unique(x$individual_id.individual_local_identifier),
                   death_date = unique(x$death_date))
  return(df)
}) %>% reduce(rbind) %>% 
  mutate(death_date = as.POSIXct(ifelse(between(date(death_date), runtime -11, runtime) | death_date > Sys.Date(), NA, death_date), tz = "UTC", origin = "1970-01-01"))

saveRDS(acc_deaths, file = paste0("C:/Users/hbronnvik/Documents/storkSSFs/estimated_death_dates_", Sys.Date(), ".rds"))

# left over we have the birds with no ACC downloading or that vanished with normal activity

# acc <- acc_data %>%
#   filter(individual_id == 2141594958) %>% #2141594958
#   mutate(timestamp = as.POSIXct(timestamp, tz = "UTC", origin = "1970-01-01"))
# 
# BurstSamplingScedule(acc)
# 
# PlotAccDataTIME(df=acc, bursts = c(1))
# 
# axesCol = grep("acceleration_axes", names(acc), value=T)
# if(nrow(acc) > 0 & length(unique(acc[, axesCol]))==1){
#   accDf_vedba <- acc %>%
#     select("individual_id", "timestamp", "event_id") %>%
#     mutate(n_samples_per_axis = NA,
#            acc_burst_duration_s=NA,
#            meanVeDBA=NA,
#            meanODBA=NA)
# 
#   accRawCol <- grep("accelerations_raw", names(acc), value=T)
#   sampFreqCol <- grep("acceleration_sampling_frequency_per_axis", names(acc), value=T)
# 
#   for(j in 1:nrow(acc)){
#     Naxes <- nchar(as.character(acc[j, axesCol]))
#     accMx <- matrix(as.integer(unlist(strsplit(as.character(acc[j, accRawCol]), " "))), ncol=Naxes, byrow = T)
#     n_samples_per_axis <- nrow(accMx)
#     acc_burst_duration_s <- n_samples_per_axis/acc[j, sampFreqCol]
#     if(nchar(acc[j, axesCol])<3){stop("The ACC data have fewer than 3 axes.")}
#     VeDBA <- sqrt((accMx[,1]-mean(accMx[,1]))^2 + (accMx[,2]-mean(accMx[,2]))^2 + (accMx[,3]-mean(accMx[,3]))^2)
#     ODBA <- (accMx[,1]-mean(accMx[,1])) + (accMx[,2]-mean(accMx[,2])) + (accMx[,3]-mean(accMx[,3]))
# 
#     accDf_vedba[j, c("n_samples_per_axis", "acc_burst_duration_s", "meanVeDBA", "meanODBA")] <- c(n_samples_per_axis, acc_burst_duration_s, mean(VeDBA, na.rm=T), mean(ODBA, na.rm=T))
#   }
# }
# 
# diagnostic <- accDf_vedba %>%
#   mutate(date = date(as.POSIXct(timestamp, origin = "1970-01-01", tz = "UTC"))) %>%
#   group_by(date) %>%
#   summarize(activity = mean(meanVeDBA)) %>%
#   mutate(active = ifelse(activity > a_thresh, 1, 0),
#          # add on whether the activity was high or low yesterday
#          check_for_event = lag(active) == T,
#          # and then say whether this is a two day long inactive period
#          cumu_check_for_event = ifelse(active == 0 & lag(active == 0), T, F))
# 
# 
# ggplot(accDf_vedba, aes(timestamp, meanVeDBA, color = meanVeDBA)) +
#   geom_point(alpha = 0.3) +
#   labs(x = "Date", y = "Mean VeDBA per burst") +
#   theme_classic() +
#   theme(legend.position = "none")
# 
# diagnostic
# 
# ggplot(diagnostic, aes(date, activity, color = activity)) +
#   geom_point(alpha = 0.3) +
#   labs(x = "Date", y = "Mean VeDBA per day", title = unique(acc$individual_id)) +
#   theme_classic() +
#   theme(legend.position = "none")



# one more filter for whether the last location was the last upload



# calculate DBA
# if under a threshold, call this animal dead
# save the death date

# compare the death date to the migration
# if within a certain time of the migration, it died on migration
# remove the bird

# is the last time the animal moved 50km per day within the last threshold amount of days it transmitted?

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


