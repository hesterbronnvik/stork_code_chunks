### The birds that do not contribute data to the study ("dead")
### Hester Brønnvik
### 27.10.2022


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
d_thresh <- 60*24*60*60 # the number of seconds in a day threshold for defining having died
a_thresh <- 10*2 # the average activity a bird has to be below to be dead x days 

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

# save(found_birds, file = paste0("C:/Users/hbronnvik/Documents/storkSSFs/clean_data/reduced50_15m_", Sys.Date(), ".RData"))


## identify individuals that died on migration

# find the last date of the birds' transmissions
last_known <- lapply(studies, function(x){
  info <- getMovebankAnimals(x, loginStored) %>% 
    filter(sensor_type_id == 653 & individual_id %in% found_birds$individual.id)
  return(info)
}) %>% reduce(rbind) %>% 
  select(individual_id, death_comments, timestamp_end, local_identifier)

# for birds that have death comments containing a digit date, extract that information
patterns <- c("(20)[0-9][0-9][[:punct:]]\\d{1,2}[[:punct:]]\\d{1,2}", # yyyy mm dd
              "\\d{1,2}[[:punct:]]\\d{1,2}[[:punct:]](20)[0-9][0-9]", # dd mm yy, mm dd yyyy, d m yyyy, m d yyyy
              "[0-9]{2}[[:punct:]][0-9]{2}[[:punct:]][0-9]{2}") # dd mm yy

death_dates <- last_known %>% 
  # remove km distances because they do not parse but m may be valid 
  mutate(death_comments = sub("Sept |Sept\\. ", "September ", gsub("\\d{1,2}( km|km)", "", death_comments))) %>% 
  rowwise() %>% 
  # extract dates
  mutate(char_date = gsub("\\.", "-", str_extract(death_comments, patterns[str_detect(death_comments, patterns)])[1]),
         # format the dates
         date = parse_date_time(char_date, orders = c("%d-%m-%Y", "%Y-%m-%d"), tz = "UTC"),
         # be careful of mdy and clean up two digit dates
         date = as.POSIXct(ifelse(length(grep("/", char_date)) == 1, parse_date_time(char_date, "%m/%d/%Y", tz = "UTC"), 
                                  ifelse(nchar(char_date) == 8, parse_date_time(char_date, "%d-%m-%y", tz = "UTC"), date)),
                           tz = "UTC", origin = "1970-01-01")) %>% 
  ungroup()

death_dates <- death_dates %>% 
  rowwise() %>% 
  # the dates that have months
  mutate(date = as.POSIXct(ifelse(is.na(date), parse_date_time(death_comments, "%d %b %y", tz = "UTC"), date), tz = "UTC", origin = "1970-01-01"))

# give up and do the last group manually because writing code for each typo is too onerous
deaths <- data.frame(individual_id = c(78032151, 1173989698, 218609200, 173666935, 23465786, 23466747, 24564548, 1176048179, 1573095923, 1576780467, 1576789367, 2165039729),
                     date = c("2015-09-08", "2020-09-12", "2015-08-03", "2017-06-27", "2015-05-31", "2015-01-24", "2015-02-03", "2021-04-01", "2021-09-11", "2021-09-01", "2021-09-01", "2022-08-23"))

death_dates <- death_dates %>% 
  rowwise() %>% 
  mutate(date = as.POSIXct(ifelse(individual_id %in% deaths$individual_id, deaths$date[which(deaths$individual_id == individual_id)], date), tz = "UTC", origin = "1970-01-01"))


# remove those birds from the search and continue
missing_dates <- death_dates %>% 
  filter(is.na(date))


# download the ACC 
info <- lapply(studies, function(x){
  info <- getMovebankAnimals(x, loginStored) %>% 
    filter(sensor_type_id == 2365683 & individual_id %in% missing_dates$individual_id) %>% 
    select("individual_id", "timestamp_start", "timestamp_end") %>% 
    mutate(pull_start = paste0(gsub("[[:punct:]]| ", "", as.POSIXct(sub(".000", "", timestamp_end), tz = "UTC", origin = "1970-01-01")-d_thresh), "000"),
           study = x)
  
}) %>% reduce(rbind)

# acc_df <- getMovebankNonLocationData(study = 212096177, animalName = 219397733, sensorID=2365683,  login=loginStored)
# x <- 5 # 42
acc_data <- lapply(1:nrow(info), function(x)tryCatch({
  df <- info[x,]
  # download the data for each individual
  acc_df <- getMovebankNonLocationData(study = df$study, animalName = df$individual_id,
                                       sensorID=2365683, timestamp_start = info$pull_start[x],
                                       login=loginStored)
  # remove the inconsistent column to allo rbinding
  if("manually_marked_outlier" %in% colnames(acc_df)){
    acc_df <- acc_df %>% 
      select(-"manually_marked_outlier")
  }
  print(ACCtimeRange(acc_df, units="days"), quotes = F)
  return(acc_df)
}, error = function(msg){ # if an error is thrown (no ACC data), continue and record the ID
  acc_df <- data.frame(matrix(NA,
                              nrow = 1,
                              ncol = 20))
  colnames(acc_df) <- c("individual_id", "deployment_id", "tag_id", "study_id", "sensor_type_id", 
                        "individual_local_identifier", "tag_local_identifier", "individual_taxon_canonical_name", 
                        "data_decoding_software", "eobs_acceleration_axes", "eobs_acceleration_sampling_frequency_per_axis", 
                        "eobs_accelerations_raw", "eobs_key_bin_checksum", "eobs_start_timestamp", 
                        "import_marked_outlier", "timestamp", "event_id",  "visible", "study_name", 
                        "sensor_type")
  acc_df$individual_id <- df$individual_id
  return(acc_df)
})) %>% reduce(rbind)

# save(acc_data, file = "C:/Users/hbronnvik/Documents/storkSSFs/acc_data/acc_94ids_31.10.2022.RData")
load("C:/Users/hbronnvik/Documents/storkSSFs/acc_data/acc_94ids_31.10.2022.RData")


todo <- unique(acc_data$individual_id[!is.na(acc_data$timestamp)])
todo <- todo[todo %in% missing_dates$individual_id[which(date(missing_dates$timestamp_end) < date(Sys.Date() - 3))]]
start_time <- Sys.time()
# for each individual in missing dates, calculate the DBA row wise and store it
event_DBAs <- lapply(todo, function(x)tryCatch({
  acc <- acc_data %>% 
    filter(individual_id == x)
  
  axesCol = grep("acceleration_axes", names(acc), value=T)
  if(nrow(acc) > 0 & length(unique(acc[, axesCol]))==1){
    accDf_vedba <- acc %>% 
      select("individual_id", "timestamp", "event_id") %>% 
      mutate(n_samples_per_axis = NA, 
             acc_burst_duration_s=NA, 
             meanVeDBA=NA,
             meanODBA=NA)
    
    accRawCol <- grep("accelerations_raw", names(acc), value=T)
    sampFreqCol <- grep("acceleration_sampling_frequency_per_axis", names(acc), value=T)
    
    for(j in 1:nrow(acc)){
      Naxes <- nchar(as.character(acc[j, axesCol]))
      accMx <- matrix(as.integer(unlist(strsplit(as.character(acc[j, accRawCol]), " "))), ncol=Naxes, byrow = T)
      n_samples_per_axis <- nrow(accMx)
      acc_burst_duration_s <- n_samples_per_axis/acc[j, sampFreqCol]
      if(nchar(acc[j, axesCol])<3){stop("The ACC data have fewer than 3 axes.")}
      VeDBA <- sqrt((accMx[,1]-mean(accMx[,1]))^2 + (accMx[,2]-mean(accMx[,2]))^2 + (accMx[,3]-mean(accMx[,3]))^2)
      ODBA <- (accMx[,1]-mean(accMx[,1])) + (accMx[,2]-mean(accMx[,2])) + (accMx[,3]-mean(accMx[,3]))
      
      accDf_vedba[j, c("n_samples_per_axis", "acc_burst_duration_s", "meanVeDBA", "meanODBA")] <- c(n_samples_per_axis, acc_burst_duration_s, mean(VeDBA, na.rm=T), mean(ODBA, na.rm=T))
    }
  }

  return(accDf_vedba)
}, error = function(msg){
  accDf_vedba <- data.frame(matrix(NA, nrow = 1, ncol = 6))
  colnames(accDf_vedba) <- c("timestamp", "event_id", "n_samples_per_axis", "acc_burst_duration_s", "meanVeDBA", "meanODBA" )
  accDf_vedba$individual_id <- x
  return(accDf_vedba)})
)
Sys.time() - start_time

# remove the empty elements (axes fewer than 3)
event_DBAs <- event_DBAs[lapply(event_DBAs,nrow) > 1]

# using the threshold set above, determine the inactivity of each individual
# if the animal is presumed dead, return the first day of inactivity 
found_dates <- lapply(event_DBAs, function(x){
  diagnostic <- x %>% 
    mutate(date = date(as.POSIXct(timestamp, origin = "1970-01-01", tz = "UTC"))) %>% 
    group_by(date) %>% 
    summarize(activity = mean(meanVeDBA)) %>% 
    mutate(active = ifelse(activity > 10, 1, 0),
           check_for_event = lag(active) == T,
           cumu_check_for_event = ifelse(active == 0 & lag(active == 0), T, F))
  
  if(sum(diagnostic[(nrow(diagnostic)- 1):nrow(diagnostic), "activity"]) < a_thresh){
    dod <- diagnostic %>% 
      filter(cumu_check_for_event == T) %>% 
      slice(1)
    ddf <- data.frame(individual_id = unique(x$individual_id), char_date = as.character(dod$date))
  }else{ddf <- data.frame(individual_id = unique(x$individual_id), char_date = NA)}
  return(ddf)
}) %>% reduce(rbind)

# add these estimated dates to the rest
death_dates <- death_dates %>% 
  rowwise() %>% 
  mutate(date = as.POSIXct(ifelse(individual_id %in% found_dates$individual_id, found_dates$char_date[which(found_dates$individual_id == individual_id)], date), tz = "UTC", origin = "1970-01-01"))

# move on without them
missing_dates <- death_dates %>% 
  filter(is.na(date) & date(timestamp_end) < Sys.Date()-2)


# BurstSamplingScedule(acc_df2)
# 
# PlotAccDataTIME(df=acc_df2, bursts = c(100))


ggplot(accDf_vedba, aes(timestamp, meanVeDBA, color = meanVeDBA)) +
  geom_point(alpha = 0.3) +
  labs(x = "Date", y = "Mean VeDBA per burst") +
  theme_classic() +
  theme(legend.position = "none")

diagnostic <- accDf_vedba %>% 
  mutate(date = date(as.POSIXct(timestamp, origin = "1970-01-01", tz = "UTC"))) %>% 
  group_by(date) %>% 
  summarize(activity = mean(meanVeDBA)) %>% 
  mutate(active = ifelse(activity > 10, 1, 0),
         check_for_event = lag(active) == T,
         cumu_check_for_event = ifelse(active == 0 & lag(active == 0), T, F))

diagnostic

ggplot(diagnostic, aes(date, activity, color = activity)) +
  geom_point(alpha = 0.3) +
  labs(x = "Date", y = "Mean VeDBA per burst", title = unique(acc_df$individual_id)) +
  theme_classic() +
  theme(legend.position = "none")



# one more filter for whether the last location was the last upload



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


