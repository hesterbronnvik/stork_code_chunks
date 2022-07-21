### Kernel density estimate for storks migrating across western Europe
### Hester Bronnvik
### 18.07.2022

library(move)
library(scales)
library(ks)
library(adehabitatHR)
library(lubridate)
library(tidyverse)


### retrieve the data
## this requires downloading large studies off movebank

# set wd and load login data
setwd("C:/Users/heste/Desktop/HB_storks")
load("loginStored.rdata")

# collect the names of the stork studies to use
studies <- c(21231406, 1176017658, 173641633)#24442409, 212096177, 76367850, 

start_time <- Sys.time()
locations <- lapply(studies, function(x){
  
  print(paste0("Study ", x, "."))
  # get the names of the animals to download
  birds <- getMovebankAnimals(study =  x, login = loginStored)
  # remove birds that had no data 
  bird_names <- unique(birds$individual_id[which(birds$number_of_events > 0)])
  
  individual_locs <- lapply(bird_names, function(y){
    print(paste0("Individual ", y, "."))
    # get the data themselves
    locations1 <- getMovebankLocationData(study = x, animalName =  y, sensorID = "GPS", login = loginStored) %>% 
      # make the move object manipulable
      as.data.frame(row.names = NULL) %>% 
      # remove missing locations
      drop_na(location.long) %>% 
      # create a column to ignore year
      mutate(datestamp = timestamp) %>% 
      # pare the data down to reduce burden and remove uselss columns
      dplyr::select("timestamp", "location.long", "location.lat", "individual.id", "individual.local.identifier", 
                    "ground.speed", "heading", "height.above.ellipsoid", "study.name", "datestamp") %>% 
      # add columns by which to categorize the data
      group_by(date(timestamp)) %>% 
      # the Haversine distance between first and last locations of the day
      mutate(daily_dist = distHaversine(c(head(location.long,1), head(location.lat, 1)), c(tail(location.long,1), tail(location.lat, 1))),
             # the rhumbline bearing between first and last locations of the day
             daily_direction = bearingRhumb(c(head(location.long,1), head(location.lat, 1)), c(tail(location.long,1), tail(location.lat, 1))),
             # finally a binary category denoting migratory or not, it is unnecessary but simplifies code
             migratory = ifelse(daily_dist > 1e5, 1, 0),
             compass_direction = ifelse(daily_direction > 90 & daily_direction < 270, "southward", "northward"), 
             phase = ifelse(migratory == 1 & compass_direction == "southward", "autumn_migration", 
                            # migrating north 
                            ifelse(migratory == 1 & compass_direction == "northward", "spring_migration", NA))) %>% 
      ungroup() 
    # reset the year
    year(locations1$datestamp) <- 1995
    # because the idea is to capture the locations of all birds during the migration, 
    # include data from birds not migrating, but during others' migrations (arbitrary months based on a histogram):
    locations_thin <- locations1 %>% 
      filter(month(datestamp) %in% c(3,4,5,8,9))
    
    return(locations_thin)
  }) %>% reduce(rbind)
  
  # save the data for each study
  saveRDS(individual_locs, file = paste0("C:/Users/heste/Desktop/HB_storks/study_locations/", x, ".rds"))
  return(individual_locs)
}) %>% reduce(rbind)
# the run time:
Sys.time() - start_time

# a rough idea of the months in which migration happens from movements of 47 birds on both (EW) routes
# hist(month(locations$datestamp[locations$migratory == 1])) # 3,4,5,8,9

### build a kernel estimate for each subset

study1kernels <- lapply(unique(date(locations_thin$datestamp)), function(f){
  print(f)
  data_sp <- locations_thin %>%  # store the data, then create a spatial object for plotting
    dplyr::filter(date(datestamp) == f) %>% 
    dplyr::select(location.long, location.lat)
  
  coordinates(data_sp) <- ~ location.long + location.lat
  
  #H.lscv <- Hlscv(x = data_sp)
  
  kern <- kernelUD(data_sp, h = "LSCV")
  return(kern)
})

plot(getverticeshr(study1kernels[[1]]))
points(locations_thin, pch = 16, cex = .75)

### save


dummydata <- data.frame(x = c(10, 5, 10), y = c("Wind support", "Uplift", "Uplift:Migration year"), sh = c(2, 2, 2), sl = c(-2, -2, -2))
dummydata$y <- factor(dummydata$y, levels = c("Uplift:Migration year", "Uplift", "Wind support"))

ggplot(dummydata, aes(x = x, y= y)) + 
  geom_vline(xintercept = 5, lty = 2, lwd = 20) +
  geom_errorbar(aes(xmax = x + sh, xmin = x + sl), width = 0.15, lwd = 40, color = "#481F70FF") +
  geom_point(cex = 70, color = "#481F70FF") +
  labs(x = "Coefficient", y = " ") +
  xlim(0, 15) +
  theme_classic() + 
  theme(legend.position="none", 
        axis.text = element_text(color = "black"), 
        text = element_text(size = 450),
        axis.line=element_line(size=4, color = "black"),
        axis.ticks.length=unit(8, "cm"),
        axis.ticks = element_blank(),
        axis.text.x.bottom = element_blank())


