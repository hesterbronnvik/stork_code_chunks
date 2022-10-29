### The birds that do not contribute data to the study ("dead")
### Hester Brønnvik
### 27.10.2022


# required packages
library(move)
library(lubridate)
library(stringr)
library(tidyverse)

# required information
setwd("C:/Users/hbronnvik/Documents/stork_code_chunks")
load("C:/Users/hbronnvik/Documents/loginStored.rdata")
studies <- c(24442409, 212096177, 76367850, 21231406, 1176017658, 173641633)
d_thresh <- 60*24*60*60 # the number of seconds in a day threshold for defining having died

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

# remove birds for which stopover cannot be differentiated from wintering
found_birds <- found_birds %>% 
  filter(year(timestamp) != year(Sys.Date()))

## identify individuals that died on migration

# find the last date of the birds' transmissions
last_known <- lapply(studies, function(x){
  info <- getMovebankAnimals(x, loginStored) %>% 
    filter(sensor_type_id == 653 & individual_id %in% found_birds$individual.id)
}) %>% reduce(rbind)%>% 
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
deaths <- data.frame(individual_id = c(78032151, 1173989698, 218609200, 173666935, 23465786, 23466747, 24564548, 1176048179, 1573095923, 1576780467, 1576789367),
                     date = c("2015-09-08", "2020-09-12", "2015-08-03", "2017-06-27", "2015-05-31", "2015-01-24", "2015-02-03", "2021-04-01", "2021-09-11", "2021-09-01", "2021-09-01"))

death_dates <- death_dates %>% 
  rowwise() %>% 
  mutate(date = as.POSIXct(ifelse(individual_id %in% deaths$individual_id, deaths$date[which(deaths$individual_id == individual_id)], date), tz = "UTC", origin = "1970-01-01"))


# remove those birds from the search and continue
missing_dates <- death_dates %>% 
  filter(is.na(date))


# download the ACC 
library(moveACC)
info <- getMovebankAnimals(studies[1], loginStored) %>% 
  filter(sensor_type_id == 2365683 & individual_id %in% missing_dates$individual_id) %>% 
  select("individual_id", "timestamp_start", "timestamp_end") %>% 
  mutate(pull_start = gsub("[[:punct:]]| ", "", as.POSIXct(sub(".000", "", timestamp_end), tz = "UTC", origin = "1970-01-01") - d_thresh))
# the last location of the first bird is one year earlier than the last GPS location
# the second bird has no ACC
acc_df <- getMovebankNonLocationData(study = studies[1], animalName = info$individual_id[3],
                                     sensorID=2365683,  login=loginStored)
ACCtimeRange(acc_df, units="days")
acc_df2 <- acc_df %>% 
  filter(timestamp > (timestamp[n()] - d_thresh))
ACCtimeRange(acc_df2, units="days")
# BurstSamplingScedule(acc_df2)
# 
# PlotAccDataTIME(df=acc_df2, bursts = c(100))

# timestamp_start = info$pull_start[3],
# timestamp_end = info$timestamp_end[3],
axesCol = grep("acceleration_axes", names(acc_df2), value=T)
if(nrow(acc_df2)>0 & length(unique(acc_df2[, axesCol]))==1){
  accDf_vedba <- acc_df2 %>% 
    select("timestamp", "event_id") %>% 
    mutate(n_samples_per_axis = NA, 
           acc_burst_duration_s=NA, 
           meanVeDBA=NA,
           meanODBA=NA)
  
  accRawCol <- grep("accelerations_raw", names(acc_df2), value=T)
  sampFreqCol <- grep("acceleration_sampling_frequency_per_axis", names(acc_df2), value=T)
  
  for(j in 1:nrow(acc_df2)){
    Naxes <- nchar(as.character(acc_df2[j, axesCol]))
    accMx <- matrix(as.integer(unlist(strsplit(as.character(acc_df2[j, accRawCol]), " "))), ncol=Naxes, byrow = T)
    n_samples_per_axis <- nrow(accMx)
    acc_burst_duration_s <- n_samples_per_axis/acc_df2[j, sampFreqCol]
    if(nchar(acc_df2[j, axesCol])<3){stop("The ACC data have fewer than 3 axes.")}
    VeDBA <- sqrt((accMx[,1]-mean(accMx[,1]))^2 + (accMx[,2]-mean(accMx[,2]))^2 + (accMx[,3]-mean(accMx[,3]))^2)
    ODBA <- (accMx[,1]-mean(accMx[,1])) + (accMx[,2]-mean(accMx[,2])) + (accMx[,3]-mean(accMx[,3]))
    
    accDf_vedba[j, c("n_samples_per_axis", "acc_burst_duration_s", "meanVeDBA", "meanODBA")] <- c(n_samples_per_axis, acc_burst_duration_s, mean(VeDBA, na.rm=T), mean(ODBA, na.rm=T))
  }
}

ggplot(accDf_vedba, aes(timestamp, meanVeDBA, color = meanVeDBA)) +
  geom_point(alpha = 0.3) +
  labs(x = "Date", y = "Mean VeDBA per burst") +
  theme_classic() +
  theme(legend.position = "none")

diagnostic <- accDf_vedba %>% 
  group_by(date(timestamp)) %>% 
  summarize(activity = mean(meanVeDBA))

diagnostic

ggplot(diagnostic, aes(`date(timestamp)`, activity, color = activity)) +
  geom_point(alpha = 0.3) +
  labs(x = "Date", y = "Mean VeDBA per burst") +
  theme_classic() +
  theme(legend.position = "none")

last_known$timestamp_end[which(last_known$individual_id == unique(acc_df2$individual_id))]
last_known$death_comments[which(last_known$individual_id == unique(acc_df2$individual_id))]


# if the DBA was below the threshold for a given amount of time
# but if that time is within the last read, then it shows as not having died
# is that to be expected (not dead) or fixed?



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


