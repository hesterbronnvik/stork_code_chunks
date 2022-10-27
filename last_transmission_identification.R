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
}) %>% reduce(rbind)
# for birds that have death comments containing a date, extract that information
# extract comments with patterns of two digits, punctuation, two digits, punctuation, two digits
death_dates <- last_known[grep("[0-9]{2}[[:punct:]][0-9]{2}[[:punct:]][0-9]{2}",last_known$death_comments),c("individual_id", "death_comments")]
# extract only the digits and no words
death_dates$dd <- gsub("\\.|/", "-", gsub("[^0-9./|-]", "", death_dates$death_comments))
# add a column of just the year month day
death_dates$date <- str_extract(death_dates$dd,"(20)[0-9][0-9][[:punct:]]\\d{2}[[:punct:]]\\d{2}")
# or of the day month year
death_dates$date2 <-str_extract(death_dates$dd,"\\d{2}[[:punct:]]\\d{2}[[:punct:]](20)[0-9][0-9]") 
# combine those two formats into one date column
death_dates$date <- as.Date(ifelse(is.na(death_dates$date), as.Date(as.character(death_dates$date2), format = "%d-%m-%y"), as.Date(death_dates$date)), origin = '1970-01-01')
death_dates$date <- as.Date(ifelse(is.na(death_dates$date), as.Date(as.character(death_dates$date2), format = "%m-%d-%y"), as.Date(death_dates$date)), origin = '1970-01-01')
# catch the times that years were entered as only two digits and not four
death_dates$dd <- ifelse(is.na(death_dates$date), gsub("-20", "-2020", death_dates$dd), death_dates$dd)
death_dates$date2 <-str_extract(death_dates$dd,"\\d{2}[[:punct:]]\\d{2}[[:punct:]](20)[0-9][0-9]") 
death_dates$date <- as.Date(ifelse(is.na(death_dates$date), as.Date(as.character(death_dates$date2), format = "%d-%m-%y"), as.Date(death_dates$date)), origin = '1970-01-01')
death_dates <- death_dates[,c("individual_id", "date")]


# remove those birds from the search and continue
missing_dates <- last_known %>% 
  filter(!individual_id %in% death_dates$individual_id)

mno <- data.frame(month = c("January", "February", "March", "April", "May", "June", "July",
                            "August", "September", "October", "November", "December", "Jan", 
                            "Feb", "Mar", "Apr", "May", "Jun", "Jul","Aug", "Sep", "Oct", "Nov", "Dec", "Sept"), 
                  no = c(1:12, 1:12, 9))

md <- missing_dates %>% 
  # the ones with death comments that have no date in a format caught above
  filter(nchar(death_comments) > 0) %>% 
  select("individual_id", "death_comments", "timestamp_end") %>% 
  rowwise() %>% 
  # take out the month if one was recorded in the comment
  mutate(month = str_extract(death_comments, "January|Jan|February|Feb|March|Mar|April|Apr|May|June|July|August|Aug|September|Sept|Sep.|October|Oct|November|Nov|December|Dec"),
         m = ifelse(!is.na(month), mno$no[which(mno$month == month)], NA),
         # take out the year
         year = str_extract(death_comments, "(20)\\d{2}"),
         date = as.Date(ifelse(paste0(m, year) == paste0(month(timestamp_end), year(timestamp_end)), date(timestamp_end), NA), origin = '1970-01-01')) %>% 
  ungroup()

# add these dates to the death confirmation group and move on without them
death_dates <- rbind(death_dates, md[,c("individual_id", "date")]) %>% 
  drop_na(date)

# with only the ids for which still no date is known, extract other date formats
md <- md %>% 
  filter(!individual_id %in% death_dates$individual_id) %>% 
  mutate(date = as.Date(str_extract(death_comments,"\\d{2}[[:punct:]]\\d{1}[[:punct:]](20)[0-9][0-9]|\\d{1}[[:punct:]]\\d{1}[[:punct:]](20)[0-9][0-9]|\\d{1}[[:punct:]]\\d{2}[[:punct:]](20)[0-9][0-9]"), format = "%d.%m.%y", origin = '1970-01-01'))

# add these dates to the death confirmation group and move on without them
death_dates <- rbind(death_dates, md[,c("individual_id", "date")]) %>% 
  drop_na(date)

# now there is no year in the death comment, so only the year from the end timestamp can be used
md <- md %>% 
  filter(!individual_id %in% death_dates$individual_id) %>% 
  mutate(year = year(timestamp_end), 
         date = as.Date(ifelse(paste0(m, year) == paste0(month(timestamp_end), year(timestamp_end)), date(timestamp_end), NA), origin = '1970-01-01'))

# add these dates to the death confirmation group and move on without them
death_dates <- rbind(death_dates, md[,c("individual_id", "date")]) %>% 
  drop_na(date)

# there are only observations when the last timestamp and the comment date are not in the same month
# or when there is no day in the comment, or when there was no death but transmission stopped
# fixing these requires manual work because each one is entered differently (typos, etc.)
# where a date is simply unclear, I have chosen to investigate it again below
md <- md %>% 
  filter(!individual_id %in% death_dates$individual_id) %>% 
  mutate(date = c(NA, NA, NA, NA, NA, "2016-04-30", "2017-06-27", 
                  NA, NA, NA, "2014-09-13", NA , NA, NA, NA, NA))

# add these dates to the death confirmation group and move on without them
death_dates <- rbind(death_dates, md[,c("individual_id", "date")]) %>% 
  drop_na(date)

# for remaining birds, use the last date
missing_dates <- last_known %>% 
  filter(!individual_id %in% death_dates$individual_id)
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