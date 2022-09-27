#-------------------------------------------------------#
#                                                       #
#--------------- Function SplitGPSBursts ---------------#
#                                                       #
#-------------------------------------------------------#

# Project: stork_code_chunks
# Authors: Hester Bronnvik, Andrea Flack & Iris Bontekoe (code originally written by Dr. Flack, Iris adjusted the code and added comments, Hester readjusted it to work on move objects and for legibility)
# Date started: 14 May 2020
# Date last modified: 27th September 2022
# R version: 4.2.1
# Description: This function identifies GPS bursts & gives each burst a unique ID
# To test:
# mv <- getMovebankData(study = 24442409, animalName =  24443630, sensorID = "GPS", # 24443411
#                       login = loginStored, removeDuplicatedTimestamps = T)
# split_gps_bursts(mv, 1, 120)

split_gps_bursts <- function(data.ind, MaxTimeDiff, MinBurstLength){ # Start function SplitGPSBursts
  
  # Load packages and functions where necessary
  if(!require("lubridate")){install.packages("lubridate");library(lubridate)}
  if(!require("raster")){install.packages("raster");library(raster)}
  if(!require("data.table")){install.packages("data.table");library(data.table)}
  if(!require("tidyverse")){install.packages("tidyverse");library(tidyverse)}
  
  source("nearest.R")
  
  #---------------------------------#
  #- Preparation of the data frame -#
  #---------------------------------#
  ind_ID <- data.ind@idData$individual_id
  
  data.ind$timeLag <- c(timeLag(data.ind, units = "secs"), NA)
  data.ind$newBurst <- ifelse(round(data.ind$timeLag) <= MaxTimeDiff, F, T)
  # shift the burst column down one so that the times are from one to lead not lag to one
  data.ind$newBurst <- c(T, data.ind$newBurst[1:(nrow(data.ind)-1)])
  # check whether the previous location was part of this burst
  data.ind$check_for_event <- dplyr::lag(as.character(data.ind$newBurst), default = "") == T
  # take the cumulative sum to act as a unique ID for each burst identified in line 37
  data.ind$cumu_check_for_event <- cumsum(data.ind$newBurst)
  
  
  #--------------------------------------------------------#
  #- Assign BurstIDs to DataFrame & Reduce data to bursts -#
  #--------------------------------------------------------#
  
  data.ind <- data.ind %>% 
    # allow grouping
    as.data.frame() %>% 
    group_by(cumu_check_for_event) %>% 
    # for the bursts, calculate the time difference between the last and first locations
    mutate(burstLength = ifelse(n()>1, difftime(tail(timestamp,1), head(timestamp, 1), units = "secs"), NA),
           # add an ID to each burst
           burstID = cur_group_id()) %>%
    # remove the bursts that do not meet user-set criteria
    filter(burstLength > MinBurstLength) %>% 
    ungroup() %>% 
    # clean up the sorting columns
    dplyr::select(-"newBurst", -"check_for_event", -"cumu_check_for_event")
  
  # save it all back to a move object
  suppressWarnings(if(nrow(data.ind)>0){data <- move(x = data.ind$location_long, y = data.ind$location_lat, 
                                    time = data.ind$timestamp, proj = CRS("+proj=longlat +ellps=WGS84"), 
                                    data = data.ind);print(paste0("Data for ",
                                    ind_ID," contain only bursts and each burst has a unique identifier."), 
                                    quote = F)}else{data <- data.ind;print(paste0(ind_ID,
                                    " has no bursts meeting the defined criteria."))})
  
  
  #--------------------------------------------#
  #- Save data frame and print end-statements -#
  #--------------------------------------------#
  
  # Save data to the R-environment
  assign("data_mv", data, pos=1)
} # End function SplitGPSBursts
