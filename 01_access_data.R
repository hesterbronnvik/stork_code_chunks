### Access white stork data stored on Movebank
### Hester Bronnvik
### 2024-10-26
### hbronnvik@ab.mpg.de

library(move)
library(tidyverse)

# required information
load("/home/hbronnvik/Documents/storkSSFs/loginStored.RData") # Movebank credentials
studies <- c(24442409, 212096177, 76367850, 21231406, 1176017658, 173641633) # Movebank study IDs

# retrieve names and transmission end-dates for available GPS data (496 IDs)
info <- lapply(studies, function(x){
  print(x)
  birds <- getMovebankReferenceTable(study = x, login = loginStored) %>%
    drop_na(animal_id) %>%
    filter(sensor_type_id == 653)
  if("animal_life_stage" %in% colnames(birds)){
    chicks <- birds %>% 
      filter(grepl("juv|chick|nestling", animal_life_stage, ignore.case = T) & grepl("release", animal_comments, ignore.case = T) == F) %>% 
      dplyr::select(animal_id, study_id, animal_local_identifier, timestamp_end)
    juv <- birds %>% 
      filter(animal_life_stage == "" & grepl("release|adult", animal_comments, ignore.case = T) == F) %>% 
      dplyr::select(animal_id, study_id, animal_local_identifier, timestamp_end)
    chicks <- rbind(chicks, juv)
  }else{
    chicks <- birds %>% 
      filter(!grepl("release|adult", birds$animal_comments, ignore.case = T)) %>% 
      dplyr::select(animal_id, study_id, animal_local_identifier, timestamp_end)
  }
  return(chicks)
}) %>% reduce(rbind)

# download GPS data, 408, 489
lapply(1:nrow(info), function(x){
  print(x)
  df <- getMovebankLocationData(info$study_id[x], sensorID = 653,
                                animalName = info$animal_id[x], login = loginStored)
  saveRDS(df, file = paste0("/home/hbronnvik/Documents/WS_data/GPS_data/", info$animal_id[x], "_",gsub("-", "", Sys.Date()),"_full_GPS_data.rds"))
})

# retrieve names and transmission end-dates for available ACC data (474 IDs)
info <- lapply(studies, function(x){
  print(x)
  birds <- getMovebankReferenceTable(study = x, login = loginStored) %>%
    drop_na(animal_id) %>% 
    filter(grepl("acceleration", sensor_type_ids))
  if("animal_life_stage" %in% colnames(birds)){
    chicks <- birds %>% 
      filter(grepl("juv|chick|nestling", animal_life_stage, ignore.case = T) & grepl("release", animal_comments, ignore.case = T) == F) %>% 
      dplyr::select(animal_id, study_id, animal_local_identifier, timestamp_end)
    juv <- birds %>% 
      filter(animal_life_stage == "" & grepl("release|adult", animal_comments, ignore.case = T) == F) %>% 
      dplyr::select(animal_id, study_id, animal_local_identifier, timestamp_end)
    chicks <- rbind(chicks, juv)
  }else{
    chicks <- birds %>% 
      filter(!grepl("release|adult", birds$animal_comments, ignore.case = T)) %>% 
      dplyr::select(animal_id, study_id, animal_local_identifier, timestamp_end)
  }
  return(chicks)
}) %>% reduce(rbind)

info <- info %>% 
  group_by(animal_id) %>% 
  slice(1) %>% 
  ungroup()

# download ACC data
acc_data <- lapply(1:nrow(info), function(n){
  print(n)
  acc_df <- getMovebankNonLocationData(info$study_id[n], "Acceleration", info$animal_id[n], loginStored)
  saveRDS(acc_df, file = paste0("/home/hbronnvik/Documents/WS_data/ACC_data/", info$animal_id[n], "_",gsub("-", "", Sys.Date()),"_full_ACC_data.rds"))
})




