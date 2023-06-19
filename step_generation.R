### Step generation for step selection functions using only the complete routes
### Hester Bronnvik
### hbronnvik@ab.mpg.de
### 2023-02-08

library(move)
library(CircStats)
library(circular)
library(fitdistrplus)
library(lubridate)
library(tidyverse)

hr <- 60 #minutes; determine the sub-sampling interval
tolerance <- 15 #minutes; tolerance for sub-sampling
n_alt <- 100 #number of alternative steps.
meters_proj <- CRS("+proj=moll +ellps=WGS84")#Mollweide projection (in meters) for accurate calculation of length
wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")

# the location data during migration
full_ml <- readRDS("C:/Users/hbronnvik/Documents/storkSSFs/juv_complete_migration_locations_40km70km30daySpeed_2023-06-12.rds")
# ml_compass <- readRDS("C:/Users/hbronnvik/Documents/storkSSFs/migration_locations_40kmCompass_2023-02-15.rds")

# only migratory flight locations
# filter ground speeds and daily distances (speed is scalar and a few missing data can affect it dramatically)
full_ml <- full_ml %>% 
  filter(daily_dist > 40*1000 & ground_speed_15 > 2)

check_n <- full_ml %>% 
  group_by(trackID, date(timestamp)) %>% 
  slice(1) %>% 
  ungroup() %>% 
  group_by(trackID) %>% 
  count() %>% 
  arrange(n)

ggplot(full_ml %>% filter(trackID == "173199625_fall_2016"), aes(location.long, location.lat, color = timestamp)) +
  borders(xlim = c(-10, 10), ylim = c(5, 50)) +
  geom_point() +
  theme_classic()

# 
# trk <- full_ml %>% 
#   make_track(location.long, location.lat, timestamp, id = individual.id) %>% 
#   nest(data = -"id") %>% 
#   mutate(steps = map(data, function(x) 
#     x %>% track_resample(rate = minutes(hr), tolerance = minutes(tolerance)) %>% steps_by_burst()))
# 
# sl_plot <- trk %>% select(id, steps) %>% unnest(cols = steps) %>% 
#   ggplot(aes(sl_, fill = factor(id))) +
#   geom_density(alpha = 0.4) +
#   labs(x = "Step length (degrees)", y = "Density") +
#   theme_classic() +
#   scale_x_continuous(limits = c(0, 2), expand = c(0, 0)) +
#   scale_y_continuous(limits = c(0, 25), expand = c(0, 0)) +
#   theme(legend.position = "none")
# 
# ta_plot <- trk %>% select(id, steps) %>% unnest(cols = steps) %>% 
#   ggplot(aes(ta_, fill = factor(id))) +
#   geom_density(alpha = 0.4) +
#   labs(x = "Turning angle (radians)", y = "Density") +
#   theme_classic() +
#   scale_y_continuous(limits = c(0, 15), expand = c(0, 0)) +
#   theme(legend.position = "none")
# 
# ggarrange(sl_plot, ta_plot)
# 
# tracks <- lapply(split(full_ml, full_ml$track_id), function(x){
#   trk <- x %>%
#     make_track(location.long, location.lat, timestamp, id = individual.id) %>%
#     track_resample(rate = minutes(hr), tolerance = minutes(tolerance)) %>% 
#     steps_by_burst() %>% 
#     random_steps(n_control = 500) %>% 
#     mutate(track_id = unique(x$track_id))
# })
# 
# check <- lapply(split(full_ml, full_ml$track_id), function(x){
#   if(nrow(x) < 10){return(nrow(x))}
# }) %>% unlist()
# 
# move_ls <- lapply(split(full_ml, full_ml$individual.id), function(x){
#   x <- x %>% 
#     arrange(timestamp)
#   mv <- move(x$location.long, x$location.lat, x$timestamp, x, proj = CRS("+proj=longlat +datum=WGS84 +no_defs"), animal = unique(x$individual.id))
# })

full_ml <- full_ml %>% 
  arrange(individual.id, timestamp) %>% 
  as.data.frame()

mv <- move(x = full_ml$location.long, y = full_ml$location.lat, time = full_ml$timestamp, data =full_ml, 
           proj = CRS("+proj=longlat +datum=WGS84 +no_defs"), animal = full_ml$trackID)

### Adapted from E. Nourani 2021 https://github.com/mahle68/global_seascape_public/blob/main/step_generation.R
#import required functions

NCEP.loxodrome.na <- function (lat1, lat2, lon1, lon2) {
  deg2rad <- pi/180
  acot <- function(x) {
    return(atan(1/x))
  }
  lat1 <- deg2rad * lat1
  lat2 <- deg2rad * lat2
  lon1 <- deg2rad * lon1
  lon2 <- deg2rad * lon2
  deltaLon <- lon2 - lon1
  pi4 <- pi/4
  Sig1 <- log(tan(pi4 + lat1/2))
  Sig2 <- log(tan(pi4 + lat2/2))
  deltaSig <- Sig2 - Sig1
  if (deltaLon == 0 && deltaSig > 0) {
    head <- 0
  }
  else if (deltaLon == 0 && deltaSig < 0) {
    head <- 180
  }
  else if (deltaSig == 0 && deltaLon > 0) {
    head <- 90
  }
  else if (deltaSig == 0 && deltaLon < 0) {
    head <- 270
  }
  else if (deltaSig < 0 && deltaLon < 0) {
    head <- acot(deltaSig/deltaLon) * 180/pi + 180
  }
  else if (deltaSig < 0 && deltaLon > 0) {
    head <- acot(deltaSig/deltaLon) * 180/pi + 180
  }
  else if (deltaSig > 0 && deltaLon > 0) {
    head <- acot(deltaSig/deltaLon) * 180/pi
  }
  else if (deltaSig > 0 && deltaLon < 0) {
    head <- acot(deltaSig/deltaLon) * 180/pi + 360
  }
  else {
    head <-NA}
  return(head)
}

(b <- Sys.time()) 
sp_obj_ls <- lapply(split(mv), function(track){ #for each track,
  print(track@idData$trackID)

  #--STEP 1: thin the data 
  track_th <- track %>%
    thinTrackTime(interval = as.difftime(hr, units = 'mins'),
                  tolerance = as.difftime(tolerance, units = 'mins')) #the un-selected bursts are the large gaps between the selected ones
  
  #--STEP 2: assign burst IDs (each chunk of track with 1 hour intervals is one burst. Longer gaps will divide the bursts) 
  track_th$selected <- c(as.character(track_th@burstId),NA) #assign "selected" as a variable
  track_th$burst_id <-c(1,rep(NA,nrow(track_th)-1)) #define a value for first row
  
  if(nrow(track_th@data) == 1){
    track_th@data$burst_id <- track_th$burst_id
  } else {for(i in 2:nrow(track_th@data)){
    #if(i == nrow(track_th@data)){
    #  track_th@data$burst_id[i] <- NA #why?!
    #} else
    if(track_th@data[i-1,"selected"] == "selected"){
      track_th@data$burst_id[i] <- track_th@data[i-1,"burst_id"]
    } else {
      track_th@data$burst_id[i] <- track_th@data[i-1,"burst_id"] + 1
    }
  }
  }
  
  #convert back to a move object (from move burst)
  track_th <- as(track_th,"Move")
  
  #--STEP 3: calculate step lengths and turning angles 
  #sl_ and ta_ calculations should be done for each burst.
  burst_ls <- split(track_th, track_th$burst_id)
  burst_ls <- Filter(function(x) length(x) >= 3, burst_ls) #remove bursts with fewer than 3 observations
  
  burst_ls <- lapply(burst_ls, function(burst){
    burst$step_length <- c(distance(burst),NA) 
    burst$turning_angle <- c(NA,turnAngleGc(burst),NA)
    burst
  })
  
  #put burst_ls into one dataframe
  bursted_sp <- do.call(rbind, burst_ls)
  
  #reassign values
  if(length(bursted_sp) >= 1){
    bursted_sp$track <- track@idData$trackID
    bursted_sp$individual.id <- track@idData$individual.id
  }
  #bursted_sp$track<-track@idData$seg_id 
  
  return(bursted_sp)
  
}) %>% 
  Filter(function(x) length(x) > 1, .) #remove segments with no observation

if(length(sp_obj_ls) > 0){
  #--STEP 4: estimate step length and turning angle distributions
  #put everything in one df
  burst_df <- sp_obj_ls %>%  
    reduce(rbind) %>% 
    as.data.frame() %>% 
    dplyr::select(-c("coords.x1","coords.x2"))
  
  #estimate von Mises parameters for turning angles
  #calculate the averages (mu).steps: 1) convert to radians. step 2) calc mean of the cosines and sines. step 3) take the arctan. OR use circular::mean.circular
  mu <- mean.circular(rad(burst_df$turning_angle[complete.cases(burst_df$turning_angle)]))
  kappa <- est.kappa(rad(burst_df$turning_angle[complete.cases(burst_df$turning_angle)]))
  
  #estimate gamma distribution for step lengths and CONVERT TO KM!!!
  sl <- burst_df$step_length[complete.cases(burst_df$step_length) & burst_df$step_length > 0]/1000 #remove 0s and NAs
  fit.gamma1 <- fitdist(sl, distr = "gamma", method = "mle")
  
  #plot turning angle and step length distributions
  png(paste0("C:/Users/hbronnvik/Documents/storkSSFs/diagnostic_ssf_plots/", hr, "_", tolerance, ".png"),
      width = 11.5, height = 8, units = "in", res = 2000)
  par(mfrow=c(1,2))
  hist(sl, freq=F, main="", xlab = "Step length (km)")
  plot(function(x) dgamma(x, shape = fit.gamma1$estimate[[1]],
                          rate = fit.gamma1$estimate[[2]]), add = TRUE, from = 0.1, to = 150, col = "blue")
  hist(rad(burst_df$turning_angle[complete.cases(burst_df$turning_angle)]),freq=F,main="",xlab="Turning angles (radians)")
  plot(function(x) dvonmises(x, mu = mu, kappa = kappa), add = TRUE, from = -3.5, to = 3.5, col = "red")
  dev.off()
  
  #diagnostic plots for step length distribution
  png(paste0("C:/Users/hbronnvik/Documents/storkSSFs/diagnostic_ssf_plots/", hr, "_", tolerance, "_diag.png"),
      width = 11.5, height = 8, units = "in", res = 2000)
  plot(fit.gamma1)
  dev.off()
  
  
  #--STEP 5: produce alternative steps
  used_av_track <- lapply(sp_obj_ls, function(track){ #for each trackment
    
    lapply(split(track,track$burst_id),function(burst){ #for each burst,
      
      #assign unique step id
      burst$step_id <- 1:nrow(burst)
      if(kappa != Inf){
      
      
        lapply(c(2:(length(burst)-1)), function(this_point){ #first point has no bearing to calc turning angle, last point has no used endpoint.
          
          current_point<- burst[this_point,]
          previous_point <- burst[this_point-1,] #this is the previous point, for calculating turning angle.
          used_point <- burst[this_point+1,] #this is the next point. the observed end-point of the step starting from the current_point
          
          #calculate bearing of previous point
          #prev_bearing<-bearing(previous_point,current_point) #am I allowing negatives?... no, right? then use NCEP.loxodrome
          prev_bearing <- NCEP.loxodrome.na(previous_point@coords[,2], current_point@coords[,2],
                                            previous_point@coords[,1], current_point@coords[,1])
          
          current_point_m <- spTransform(current_point, meters_proj) #convert to meters proj
          
          #randomly generate n alternative points
          rnd <- data.frame(turning_angle = as.vector(rvonmises(n = n_alt, mu = mu, kappa = kappa)), #randomly generate n step lengths and turning angles
                            step_length = rgamma(n = n_alt, shape = fit.gamma1$estimate[[1]], rate = fit.gamma1$estimate[[2]]) * 1000) %>% 
            #find the gepgraphic location of each alternative point; calculate bearing to the next point: add ta to the bearing of the previous point
            mutate(lon = current_point_m@coords[,1] + step_length*cos(turning_angle),
                   lat = current_point_m@coords[,2] + step_length*sin(turning_angle))
          
          
          #convert back to lat-lon proj
          rnd_sp <- rnd
          coordinates(rnd_sp) <- ~lon+lat
          proj4string(rnd_sp) <- meters_proj
          rnd_sp <- spTransform(rnd_sp, wgs)
          
          #put used and available points together
          df <- used_point@data %>%  
            mutate(x = location.long,
                   y = location.lat) %>% 
            slice(rep(row_number(), n_alt+1)) %>% #paste each row n_alt times for the used and alternative steps
            mutate(location.long = c(head(x,1),rnd_sp@coords[,1]), #the coordinates were called x and y in the previous version
                   location.lat = c(head(y,1),rnd_sp@coords[,2]),
                   turning_angle = c(head(turning_angle,1),deg(rnd_sp$turning_angle)),
                   step_length = c(head(step_length,1),rnd_sp$step_length),
                   used = c(1,rep(0,n_alt)))  %>%
            dplyr::select(-c("x","y")) %>% 
            rowwise() %>% 
            mutate(heading = NCEP.loxodrome.na(lat1 = current_point@coords[,2], lat2 = location.lat, lon1 = current_point@coords[,1], lon2 = location.long)) %>% 
            as.data.frame()
          
          df
        
        }) %>% 
        reduce(rbind)
      }else{
          df <- burst@data %>% 
            mutate(used = NA,
                   heading = NA)
        }
      
    }) %>% 
      reduce(rbind)
    
  }) %>% 
    reduce(rbind)
}
  

Sys.time() - b # Time difference of 1.044653 hours

#create one dataframe with movebank specifications
used_av_df <- used_av_track %>% 
  dplyr::select( c("timestamp", "location.lat", "location.long", "selected", "individual.id",  "burst_id", "step_length", "turning_angle", "track", "step_id", "used", "heading")) %>% 
  # mutate(timestamp = paste(as.character(timestamp),"000",sep = ".")) %>% 
  mutate(stratum = paste(track, burst_id, step_id, sep = "_")) %>% #create unique stratum id
  as.data.frame()

# ---------------------------------------------------------------------------------------------------------

#rename lat and lon columns
colnames(used_av_df)[c(2,3)] <- c("location-lat", "location-long")

# saveRDS(used_av_df, file = "C:/Users/hbronnvik/Documents/storkSSFs/used_av_df_230612.rds")

used_av_df <- used_av_df %>% 
  mutate(splitter = c(rep("A", times = 1000000), rep("B", times = 1000000), rep("C", times = 1000000), rep("D", times = 1000000), 
                      rep("E", times = 1000000), rep("F", times = 1000000), rep("G", times = 1000000), rep("H", times = 1000000), 
                      rep("I", times = 1000000), rep("J", times = 1000000), rep("K", times = 1000000), rep("L", times = 1000000), 
                      rep("M", times = 1000000), rep("N", times = 1000000), rep("O", times = 1000000), rep("P", times = 1000000),
                      rep("Q", times = 1000000), rep("R", times = 1000000), rep("S", times = nrow(used_av_df)-18*1000000)))

# save to files for Movebank annotation
lapply(split(used_av_df, used_av_df$splitter), function(x){
  df <- x %>% 
    dplyr::select(-splitter)
  write.csv(df, file = paste0("C:/Users/hbronnvik/Documents/storkSSFs/burst_data/speed40x80/", unique(x$splitter), "_", Sys.Date(), ".csv"))
})

ang_b <- used_av_df %>% 
  filter(track == "post_2018_504439767") %>% 
  group_by(burst_id) %>% 
  mutate(burst = cur_group_id()) %>% 
  ungroup()

library(maps)
map_data <- map_data('world') %>% 
  filter(long > -20 & long < 10 & lat > 0 & lat < 60)

library(viridis)
ggplot(ang_b, aes(location.long, location.lat)) +
  # geom_raster(data = mask_df, aes(x = x, y = y, fill = tmax1)) +
  geom_polygon(data = map_data,
               aes(x=long, y=lat, group = group),
               color = "gray60", fill = "gray60") +
  # scale_fill_manual(values = c ("gray50", "white")) +
  geom_point(alpha = 0.1, aes(color = burst), data = ang_b %>% filter(used == 0)) +
  geom_point(color = "black", data = ang_b %>% filter(used == 1), cex = 0.5) +
  scale_color_viridis(option = "A") +
  labs(x = "Longitude", y = "Latitude") +
  theme_classic() +
  scale_x_continuous(limits = c(-20, 10), expand = c(0, 0)) +
  # scale_y_continuous(limits = c(12, 60), expand = c(0, 0)) +
  theme(legend.position = "none", 
        axis.text = element_text(color="black", size = 15),
        text = element_text(color="black", size = 15))

ang_b <- full_ml %>% 
  filter(track_id %in% unique(track_id)[17:27])

ggplot(ang_b, aes(location.long, location.lat, color= timestamp)) +
  borders(xlim = c(-10, 10), ylim = c(35,50), fill = "gray50") +
  geom_point() +
  # geom_point(color = "black", data = ang_b %>% filter(used == 1)) +
  scale_color_viridis(option = "A", discrete = F) +
  labs(x = "Longitude", y = "Latitude", title = unique(ang_b$track_id)) +
  theme_classic() +
  theme(legend.position = "none") +
  facet_wrap(~track_id)


check_full <- readRDS("C:/Users/hbronnvik/Documents/storkSSFs/full_data/1176048179_2023-01-28.rds")

ggplot(check_full %>% filter(year(timestamp) == 2021), aes(location.long, location.lat, color= month(timestamp))) +
  geom_point() +
  labs(x = "Longitude", y = "Latitude") +
  theme_classic()


cp <- lapply(split(migration_locations, migration_locations$track_id), function(x){
  x$month <- as.factor(month(x$timestamp))
  png(paste0("C:/Users/hbronnvik/Documents/storkSSFs/m_locs_check/", unique(x$individual.id), ".png"))
  print(ggplot(x, aes(location.long, location.lat, color= month)) +
          borders(xlim = c(-10, 10), ylim = c(round(min(x$location.lat)),50), fill = "gray50") +
          geom_point() +
          labs(x = "Longitude", y = "Latitude", title = unique(x$individual.id)) +
          theme_classic())
  dev.off()
  
})

check <- readRDS("C:/Users/hbronnvik/Documents/storkSSFs/full_data/1176031140_2023-01-28.rds")
check2 <- readRDS("C:/Users/hbronnvik/Documents/storkSSFs/full_data/1578602426_2023-01-28.rds")
ggplot(check2, aes(location.long, location.lat, color = year(timestamp))) +
  borders(xlim = c(-10, 10), ylim = c(round(min(x$location.lat)),50), fill = "gray50") +
  geom_point() +
  labs(x = "Longitude", y = "Latitude") +
  theme_classic()


