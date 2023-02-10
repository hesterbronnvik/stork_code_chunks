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
tolerance <- 15 #minutes; tolerance for sub-samplling
n_alt <- 500 #number of alternative steps.
meters_proj <- CRS("+proj=moll +ellps=WGS84")#Mollweide projection (in meters) for accurate calculation of length
wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")

# the location data during migration
migration_locations <- readRDS("C:/Users/hbronnvik/Documents/storkSSFs/classified_gps/full_migration_locations_40km15min_2023-02-09.rds")
# define birds that took eastern routes as ones that are ever east of 15.7 longitude (East Germany)
eastern_birds <- unique(migration_locations$individual.id[migration_locations$location.long > 15.7])
# remove the eastern birds
migration_locations <- migration_locations %>% 
  filter(!individual.id %in% eastern_birds)
# the total number of migrations attempted and completed by each animal in each season 
meta <- migration_locations %>%
  mutate(season = ifelse(grepl("post", phase), "post-breeding", "pre-breeding")) %>% 
  group_by(individual.id, season) %>% 
  count(track_id) %>% 
  mutate(journey_number = row_number()) %>% 
  ungroup() %>% 
  group_by(individual.id) %>%
  mutate(total_journeys = n()) %>% 
  ungroup()
# add the number of journeys
migration_locations <- migration_locations %>% 
  rowwise() %>% 
  mutate(journey_number = meta$journey_number[which(meta$track_id == track_id)]) %>% 
  ungroup()

# the identities of the complete routes
meta <- readRDS("C:/Users/hbronnvik/Documents/storkSSFs/rs_ids_40km_2023-02-08.rds")

# only locations from complete routes
full_ml <- migration_locations %>% 
  filter(track_id %in% meta$track_id[meta$track_status == "complete"])

# only migratory flight locations
# calculate ground speeds between 15 minute locations
full_ml <- full_ml %>% 
  group_by(track_id) %>% 
  mutate(distance = distVincentyEllipsoid(cbind(location.long, location.lat), cbind(lag(location.long), lag(location.lat))),
         timediff = as.numeric(difftime(timestamp, lag(timestamp), units = "secs")),
         gr_speed = distance/timediff) %>% 
  ungroup()

# filter ground speeds and daily distances (speed is scalar and a few missing data can affect it dramatically)
full_ml <- full_ml %>% 
  filter(daily_dist > 40*1000 & gr_speed > 2)


ggplot(full_ml %>% filter(individual.id == unique(full_ml$individual.id)[1]), aes(location.long, location.lat, color = phase)) +
  borders(xlim = c(-10, 10), ylim = c(30, 48)) +
  geom_point(cex = .25, color = "black") +
  theme_classic()

trk <- full_ml %>% 
  make_track(location.long, location.lat, timestamp, id = individual.id) %>% 
  nest(data = -"id") %>% 
  mutate(steps = map(data, function(x) 
    x %>% track_resample(rate = minutes(hr), tolerance = minutes(tolerance)) %>% steps_by_burst()))

sl_plot <- trk %>% select(id, steps) %>% unnest(cols = steps) %>% 
  ggplot(aes(sl_, fill = factor(id))) +
  geom_density(alpha = 0.4) +
  labs(x = "Step length (degrees)", y = "Density") +
  theme_classic() +
  scale_x_continuous(limits = c(0, 2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 25), expand = c(0, 0)) +
  theme(legend.position = "none")

ta_plot <- trk %>% select(id, steps) %>% unnest(cols = steps) %>% 
  ggplot(aes(ta_, fill = factor(id))) +
  geom_density(alpha = 0.4) +
  labs(x = "Turning angle (radians)", y = "Density") +
  theme_classic() +
  scale_y_continuous(limits = c(0, 15), expand = c(0, 0)) +
  theme(legend.position = "none")

ggarrange(sl_plot, ta_plot)

tracks <- lapply(split(full_ml, full_ml$track_id), function(x){
  trk <- x %>%
    make_track(location.long, location.lat, timestamp, id = individual.id) %>%
    track_resample(rate = minutes(hr), tolerance = minutes(tolerance)) %>% 
    steps_by_burst() %>% 
    random_steps(n_control = 500) %>% 
    mutate(track_id = unique(x$track_id))
})

check <- lapply(split(full_ml, full_ml$track_id), function(x){
  if(nrow(x) < 10){return(nrow(x))}
}) %>% unlist()

move_ls <- lapply(split(full_ml, full_ml$individual.id), function(x){
  mv <- move(x$location.long, x$location.lat, x$timestamp, x, proj = CRS("+proj=longlat +datum=WGS84 +no_defs"), animal = unique(x$individual.id))
})

### Taken from E. Nourani 2021 https://github.com/mahle68/global_seascape_public/blob/main/step_generation.R
#import required functions
full_ml <- full_ml %>% 
  select(-ground.speed, -timediff, -distance)
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
# which(unique(full_ml$individual.id) == 2160738503)
used_av_ls <- lapply(split(full_ml, full_ml$individual.id), function(group){ # for each group (ie. unique species-flyway combination)
  print(unique(group$individual.id))  
  sp_obj_ls <- lapply(split(group, group$track_id), function(track){ #for each track,
    print(unique(track$track_id))
    # make a move object
    mv <- move(track$location.long, track$location.lat, track$timestamp, track, proj = CRS("+proj=longlat +datum=WGS84 +no_defs"), animal = unique(track$track_id))
    mv@data$gr_speed <- track$gr_speed
    mv@data$daily_dist <- track$daily_dist
    mv@data$daily_direction <- track$daily_direction
    print(ncol(mv))
    #--STEP 1: thin the data 
    track_th <- mv %>%
      thinTrackTime(interval = as.difftime(hr, units = 'mins'),
                    tolerance = as.difftime(tolerance, units = 'mins')) #the unselected bursts are the large gaps between the selected ones
    
    #--STEP 2: assign burst IDs (each chunk of track with 1 hour intervals is one burst... longer gaps will divide the brusts) 
    track_th$selected <- c(as.character(track_th@burstId),NA) #assign selected as a variable
    track_th$burst_id <-c(1,rep(NA,nrow(track_th)-1)) #define value for first row
    
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
    burst_ls <- Filter(function(x) length(x) >= 3, burst_ls) #remove bursts with less than 3 observations
    
    burst_ls <- lapply(burst_ls, function(burst){
      burst$step_length <- c(distance(burst),NA) 
      burst$turning_angle <- c(NA,turnAngleGc(burst),NA)
      burst
    })
    
    #put burst_ls into one dataframe
    bursted_sp <- do.call(rbind, burst_ls)
    
    #reassign values
    if(length(bursted_sp) >= 1){
      bursted_sp$track <- mv@idData$track_id
      bursted_sp$individual.id <- mv@idData$individual.id
    }
    #bursted_sp$track<-track@idData$seg_id 
    
    return(bursted_sp)
    
  }) %>% 
    Filter(function(x) length(x) > 1, .) #remove segments with no observation
  
  if(length(sp_obj_ls) > 0){
    #--STEP 4: estimate step length and turning angle distributions
    #put everything in one df
    bursted_df <- sp_obj_ls %>%  
      reduce(rbind) %>% 
      as.data.frame() %>% 
      dplyr::select(-c("coords.x1","coords.x2"))
    
    #estimate von Mises parameters for turning angles
    #calculate the averages (mu).steps: 1) convert to radians. step 2) calc mean of the cosines and sines. step 3) take the arctan. OR use circular::mean.circular
    mu <- mean.circular(rad(bursted_df$turning_angle[complete.cases(bursted_df$turning_angle)]))
    kappa <- est.kappa(rad(bursted_df$turning_angle[complete.cases(bursted_df$turning_angle)]))
    
    #estimate gamma distribution for step lengths and CONVERT TO KM!!! :p
    sl <- bursted_df$step_length[complete.cases(bursted_df$step_length) & bursted_df$step_length > 0]/1000 #remove 0s and NAs
    fit.gamma1 <- fitdist(sl, distr = "gamma", method = "mle")
    
    #plot turning angle and step length distributions
    pdf(paste0("C:/Users/hbronnvik/Documents/storkSSFs/diagnostic_ssf_plots/", hr, "_", tolerance, unique(group$individual.id)[1], ".pdf"))
    par(mfrow=c(1,2))
    hist(sl, freq=F, main="", xlab = "Step length (km)")
    plot(function(x) dgamma(x, shape = fit.gamma1$estimate[[1]],
                            rate = fit.gamma1$estimate[[2]]), add = TRUE, from = 0.1, to = 150, col = "blue")
    hist(rad(bursted_df$turning_angle[complete.cases(bursted_df$turning_angle)]),freq=F,main="",xlab="Turning angles (radians)")
    plot(function(x) dvonmises(x, mu = mu, kappa = kappa), add = TRUE, from = -3.5, to = 3.5, col = "red")
    dev.off()
    
    #diagnostic plots for step length distribution
    pdf(paste0("C:/Users/hbronnvik/Documents/storkSSFs/diagnostic_ssf_plots/", hr, "_", tolerance, unique(group$individual.id)[1], "_diag.pdf"))
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
    used_av_track
  }
  
})
Sys.time() - b # 17.95579 mins

used_av_ls <- used_av_ls %>% 
  Filter(function(x) length(x) > 1, .)

#create one dataframe with movebank specifications
used_av_df <- lapply(c(1:length(used_av_ls)), function(i){
  
  data <- used_av_ls[[i]] %>% 
    dplyr::select( c("timestamp", "location.lat", "location.long", "selected", "individual.id",  "burst_id", "step_length", "turning_angle", "track", "step_id", "used", "heading")) %>% 
    mutate(timestamp = paste(as.character(timestamp),"000",sep = "."), #the movebank
           group = names(used_av_ls)[[i]]) %>% 
    rowwise() %>% 
    mutate(ind = strsplit(track, "_")[[1]][1], #create variable for individual id
           stratum = paste(track, burst_id, step_id, sep = "_")) %>% #create unique stratum id
    as.data.frame()
}) %>% 
  reduce(rbind)

#rename lat and lon columns
colnames(used_av_df)[c(2,3)] <- c("location-lat", "location-long")

used_av_df <- used_av_df %>% 
  mutate(splitter = c(rep("A", times = 1000000), rep("B", times = 1000000), rep("C", times = 1000000), rep("D", times = 1000000), 
                      rep("E", times = 1000000), rep("F", times = 1000000), rep("G", times = 1000000), rep("H", times = 1000000), 
                      rep("I", times = 1000000), rep("J", times = 1000000), rep("K", times = 1000000), rep("L", times = 1000000), 
                      rep("M", times = 1000000), rep("N", times = 1000000), rep("O", times = 1000000), rep("P", times = 1000000),
                      rep("Q", times = 1000000), rep("R", 1000000), rep("S", 543513)))

# save to files for Movebank annotation
lapply(split(used_av_df, used_av_df$splitter), function(x){
  df <- x %>% 
    select(-splitter)
  write.csv(df, file = paste("C:/Users/hbronnvik/Documents/storkSSFs/burst_data/", unique(x$splitter), Sys.Date(), ".csv", sep = "_"))
})

ang_b <- used_av_df %>% 
  filter(track == unique(track)[10])

library(viridis)
ggplot(ang_b, aes(location.long, location.lat)) +
  borders(xlim = c(-10, 10), ylim = c(35,50), fill = "gray50") +
  geom_point(alpha = 0.25, aes(color = stratum), data = ang_b %>% filter(used == 0)) +
  geom_point(color = "black", data = ang_b %>% filter(used == 1)) +
  scale_color_viridis(option = "A", discrete = T) +
  labs(x = "Longitude", y = "Latitude") +
  theme_classic() +
  theme(legend.position = "none")

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



