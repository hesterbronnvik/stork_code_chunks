### Pull in Black Storks and look at their overlap with White Stork selection
### Hester Bronnvik
### 2025-08-06
### hbronnvik@ab.mpg.de

library(move)
library(tidyverse)

theme_set(theme_classic()+
            theme(axis.text = element_text(color = "black", size = 12), 
                  text = element_text(size = 15),
                  strip.background = element_blank()))

# required information
load("/home/hbronnvik/Documents/storkSSFs/loginStored.RData") # Movebank credentials
studies <- c(463766103, 291047293) # Movebank study IDs

# retrieve names and transmission end-dates for available GPS data (496 IDs)
info <- lapply(studies, function(x){
  print(x)
  birds <- getMovebankReferenceTable(study = x, login = loginStored) %>%
    drop_na(animal_id) %>%
    filter(sensor_type_id == 653 & number_of_location_events > 7) %>% 
    dplyr::select(animal_id, study_id, animal_local_identifier, deploy_on_timestamp, timestamp_end)
  return(birds)
}) %>% reduce(rbind)

# download GPS data
# lapply(1:nrow(info), function(x){
#   print(x)
  # df <- getMovebankLocationData(info$study_id[x], sensorID = 653,
  #                               animalName = info$animal_id[x], login = loginStored)
  # saveRDS(df, file = paste0("/home/hbronnvik/Documents/chapter3/black_storks/GPS_data/", info$animal_id[x], "_",gsub("-", "", Sys.Date()),"_full_GPS_data.rds"))
# })

# look at the data
full_df <- lapply(list.files("/home/hbronnvik/Documents/chapter3/black_storks/GPS_data", full.names = T), function(x){
  df <- readRDS(x) %>% 
    dplyr::select(location.lat, location.long, timestamp, individual.id, #eobs.status,
                  # eobs.horizontal.accuracy.estimate, 
                  gps.satellite.count, height.above.ellipsoid,
                  heading)
  # the legbands do not produce information on DOP or activity,
  # some very quick-and-dirty filtering is called for
  if(unique(df$individual.id) == 931941267){
    df <- df %>% 
      mutate(port = location.lat > 36 & location.long < -7) %>% 
      filter(port == F) %>% 
      dplyr::select(-port)
  }
  if(unique(df$individual.id) == 951294665){
    df <- df %>% 
      mutate(port = location.lat > 48 & location.long < 2) %>% 
      filter(port == F) %>% 
      dplyr::select(-port)
  }
  if(unique(df$individual.id) == 930265885){
    df <- df %>% 
      mutate(port = location.lat > 49 & location.long < 5) %>% 
      filter(port == F) %>% 
      dplyr::select(-port)
  }
  if(unique(df$individual.id) == 948389693){
    df <- df %>% 
      mutate(port = date(timestamp)=="2019-09-18") %>% 
      filter(port == F) %>% 
      dplyr::select(-port)
  }
  if(unique(df$individual.id) == 929584364){
    df <- df %>% 
      mutate(port = date(timestamp)=="2019-09-24") %>% 
      filter(port == F) %>% 
      dplyr::select(-port)
  }
  if(unique(df$individual.id) == 502818030){
    df <- df %>% 
      mutate(port = location.long < 5) %>% 
      filter(port == F) %>% 
      dplyr::select(-port)
  }
  return(df)
}) %>% reduce(rbind)

full_sf <- full_df %>% 
  filter(location.lat <= 58) %>% 
  group_by(individual.id) %>% 
  arrange(timestamp) %>% 
  mutate(obs = 1:n(),
         date = date(timestamp),
         td = as.numeric(difftime(timestamp, lag(timestamp), units = "secs"))) %>%
  # down-sample to 15 minutes
  mutate(seq15 = round_date(timestamp, unit = "15 minutes")) %>%
  group_by(individual.id, seq15) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(-seq15, -td) %>% 
  # calculate ground speeds
  mutate(distance = distVincentyEllipsoid(cbind(location.long, location.lat), cbind(lag(location.long), lag(location.lat))),
         timediff = as.numeric(difftime(timestamp, lag(timestamp), units = "secs")),
         ground_speed_15 = distance/timediff) %>% 
  # remove absurd speeds to clean outliers
  filter(ground_speed_15 < 50) %>% 
  st_as_sf(coords = c("location.long", "location.lat"), crs = 4326)
mapview::mapview(full_sf, zcol = "obs")

# segment each individual into migration tracks
full_sf <- full_sf %>% 
  group_by(individual.id) %>% 
  group_split()

d_thresh <- 70 # the number of km per day to count as migration (class_df below)
s_thresh <- 200 # km of space covered to include the point before/after "migration" (edges below)

segmented_tracks <- lapply(1:length(full_sf), function(i){
  bird <- full_sf[[i]]
  print(paste0("Processing ", unique(bird$individual.id), ", bird ", i, " of ", length(full_sf), "."))
  # add on some basic metrics of time and distance
  sub_sf <- bird %>% 
    mutate(date = as.character(date(timestamp))) %>% 
    # get the first location of each day
    group_by(date) %>% 
    slice(1) %>% 
    ungroup() %>% 
    mutate(location.long = st_coordinates(.)[,1],
           location.lat = st_coordinates(.)[,2],
           # see how long it has been since the last observation
           time_lag = as.numeric(difftime(timestamp, lead(timestamp), units = "day")),
           # get the number of km traveled since then
           daily_dist = distVincentyEllipsoid(cbind(location.long, location.lat),
                                              cbind(lead(location.long), lead(location.lat)))/1000,
           # get ground speeds in km/day (to match white storks)
           km_day = daily_dist/abs(time_lag),
           yr = as.character(year(timestamp)),
           # get the bearing between consecutive locations
           direction = bearing(cbind(location.long, location.lat),
                               cbind(lead(location.long), lead(location.lat))),
           # define anything below the mid-line as broadly "southward"
           compass_direction = ifelse(abs(direction) >= 90, "southward", "northward"))
  # now use the distances and directions along with time-of-year to define migration (same as white storks)
  class_df <- sub_sf %>% 
    filter(round(km_day) >= d_thresh) 
  
  if(nrow(class_df) > 0){
    class_df <- class_df %>% 
      mutate(timeLag = as.numeric(difftime(timestamp, lag(timestamp), units = "day")),
             # insert an enormous time lag for the first location to remove the NA
             timeLag = ifelse(is.na(timeLag), 1e5, timeLag),
             new_burst = ifelse(round(timeLag) <= 21, 0, 1),
             burst = cumsum(new_burst)) %>% 
      group_by(burst) %>% 
      mutate(burst_direction = bearing(cbind(location.long[1], location.lat[1]),
                                       cbind(location.long[n()], location.lat[n()])),
             # use the general direction of the burst to classify unclear locations
             burst_direction = ifelse(abs(burst_direction) >= 90, "south", "north"),
             # also get the distance of the burst to ensure rapid movement
             burst_time = as.numeric(difftime(timestamp[n()], timestamp[1], units = "day")),
             # get the number of km traveled since then
             burst_dist = distVincentyEllipsoid(cbind(location.long[1], location.lat[1]),
                                                cbind(location.long[n()], location.lat[n()]))/1000,
             # get ground speeds in km/day (to match white storks)
             burst_km_day = burst_dist/abs(burst_time),
             # Aug-Nov is "late", Jan-Jun is "early" and Jul/Dec are "change"
             obs_class = ifelse(month(timestamp[n()]) %in% c(8:11), "late", 
                                ifelse(month(timestamp[n()]) %in% c(1:6), "early", "change")),
             # base the "change" season on direction (i.e. in Dec. north is spring, but south is fall)
             obs_season = ifelse(obs_class == "early", "spring", 
                                 ifelse(obs_class == "late", "fall",
                                        ifelse(burst_direction == "north" & obs_class == "change", "spring", 
                                               ifelse(burst_direction == "south" & obs_class == "change", "fall", NA)))),
             obs_year = ifelse(month(timestamp[n()]) == 12 & obs_season == "spring", year(years(1) + timestamp[n()]), year(timestamp[n()])),
             trackID = paste(unique(individual.id), obs_season, obs_year, sep = "_")) %>% 
      ungroup() #%>% 
      # filter(round(burst_km_day) >= 21)
    
    # just the start/end of each track
    migrations <- class_df %>% 
      st_drop_geometry()%>% 
      drop_na(obs_season) %>% 
      group_by(trackID) %>% 
      slice(1, n()) %>% 
      dplyr::select(location.long, location.lat, timestamp, trackID) %>% 
      group_by(trackID) %>% 
      group_split()
    
    if(length(migrations) > 0){
      # paste onto the full data
      segged <- lapply(migrations, function(x){
        ind <- bird %>% 
          filter(between(timestamp, x$timestamp[1], x$timestamp[2])) %>% 
          mutate(track_id = unique(x$trackID)) %>% 
          st_drop_geometry()
      }) %>% reduce(rbind)
      
      segged <- bird %>%
        mutate(location.long = st_coordinates(.)[,1],
               location.lat = st_coordinates(.)[,2]) %>% 
        st_drop_geometry() %>% 
        left_join(segged) %>% 
        st_as_sf(coords = c("location.long", "location.lat"), crs = 4326)
      
      # this has reliably pulled out everything we can be confident is migration, 
      # but there are loose ends because the birds will travel 2000km in 70 days, i.e. 30 km/day
      # if the gaps are so big, we need to add in the points before/after
      # this means we have no ability to say *when* migration started/ended, but can guess *where*
      edges <- lapply(migrations, function(x){
        ind <- bird %>% 
          mutate(location.long = st_coordinates(.)[,1],
                 location.lat = st_coordinates(.)[,2])
        index <- ind %>% 
          st_drop_geometry() %>% 
          mutate(index = row_number()) %>% 
          filter(timestamp == x$timestamp[1]|timestamp == x$timestamp[2]) %>% 
          dplyr::select(index) %>% 
          deframe()
        # if there is only one point, duplicate it and treat it as start and end
        if(length(index) == 1){
          index <- rep(index, 2)
        }
        before <- ind %>% 
          slice(index[1]-1, index[1]) %>% 
          mutate(displacement = distVincentyEllipsoid(cbind(location.long, location.lat),
                                                      cbind(lead(location.long), lead(location.lat)))/1000) %>% 
          dplyr::select(displacement, timestamp, location.long, location.lat)
        
        after <- ind %>% 
          slice(index[2], index[2]+1) %>% 
          mutate(displacement = distVincentyEllipsoid(cbind(location.long, location.lat),
                                                      cbind(lag(location.long), lag(location.lat)))/1000) %>% 
          dplyr::select(displacement, timestamp, location.long, location.lat)
        
        # save the rows before and after classified migration if they meet the criterion
        if(index[1] > 1){
          if(round(before$displacement[1]) >= s_thresh){
            x$timestamp[1] <- before$timestamp[1]
            x$location.long[1] <- before$location.long[1]
            x$location.lat[1] <- before$location.lat[1]
          }
        }
        if(index[2] < nrow(ind)){
          if(round(after$displacement[2]) >= s_thresh){
            x$timestamp[2] <- after$timestamp[2]
            x$location.long[2] <- after$location.long[2]
            x$location.lat[2] <- after$location.lat[2]
          }
        }
        
        return(x)
      })
      
      # paste onto the full data
      segged <- lapply(edges, function(x){
        ind <- bird %>% 
          filter(between(timestamp, x$timestamp[1], x$timestamp[2])) %>% 
          mutate(track_id = unique(x$trackID),
                 location.long = st_coordinates(.)[,1],
                 location.lat = st_coordinates(.)[,2]) %>% 
          st_drop_geometry()
      }) %>% reduce(rbind)
      return(segged)
    }
  }
}) %>% reduce(rbind) %>% 
  separate(track_id,  c(NA, "season", NA), remove = F, sep = "_")

# finally, filter out tracks that are not long-distance (either no transmission or death)
track_lengths <- segmented_tracks %>% 
  group_by(track_id) %>% 
  summarize(dist = distVincentyEllipsoid(cbind(location.long[1], location.lat[1]),
                                         cbind(location.long[n()], location.lat[n()]))/1000) %>% 
  filter(round(dist) >= 700)

segmented_tracks <- segmented_tracks %>% 
  filter(track_id %in% track_lengths$track_id)

library(rnaturalearth)
world <- ne_countries(scale = "large", returnclass = "sf")
colfunc <- colorRampPalette(c("#143601", "#1a4301", "#245501", "#538d22", "#73a942", "#aad576"))
# png(filename = "/home/hbronnvik/Documents/chapter3/black_storks/tracks.png",
#     height = 8, width = 11, units = "in", res = 500)
ggplot() +
  geom_sf(data = world, fill = "gray70", color = "gray60") +
  geom_path(data = segmented_tracks, aes(location.long, location.lat, 
                          group = as.factor(track_id),
                          # color = as.numeric(as.factor(track_id))
  ),
  alpha = 0.5) +
  geom_point(data = segmented_tracks, aes(location.long, location.lat, 
                           color = as.numeric(as.factor(individual.id))),
             alpha = 0.5) +
  coord_sf(xlim = c(-20, 50), ylim = c(0, 60)) +
  scale_x_continuous(expand = c(0, 0), labels = function(x) as.character(x)) +
  scale_y_continuous(expand = c(0, 0), labels = function(x) as.character(x)) +
  scale_color_gradientn(colors = colfunc(100)) +
  labs(color = "ID", x = "Longitude", y = "Latitude") +
  facet_wrap(~season)
# dev.off()

segmented_tracks %>% 
  st_as_sf(coords = c("location.long", "location.lat"), crs = 4326) %>% 
  mapview::mapview(zcol = "track_id")

# 62 individuals
length(unique(segmented_tracks$individual.id))
# generating 186 tracks
length(unique(segmented_tracks$track_id))

# for (i in 1:length(unique(segmented_tracks$individual.id))) {
#   x <- segmented_tracks %>% 
#     filter(individual.id == unique(segmented_tracks$individual.id)[i])
#   png(filename = paste0("/home/hbronnvik/Documents/chapter3/black_storks/track_plots/", unique(x$individual.id),"_tracks.png"),
#       height = 11, width = 11, units = "in", res = 500)
#   print(ggplot() +
#           geom_sf(data = world, fill = "gray70", color = "gray60") +
#           geom_path(data = x, aes(location.long, location.lat, 
#                                   group = as.factor(track_id),
#                                   # color = as.numeric(as.factor(track_id))
#           ),
#           alpha = 0.5) +
#           geom_point(data = x, aes(location.long, location.lat, 
#                                    group = as.factor(track_id),
#                                    color = as.factor(track_id)),
#                      alpha = 0.5) +
#           coord_sf(xlim = c(-20, 50), ylim = c(0, 60)) +
#           scale_x_continuous(expand = c(0, 0), labels = function(x) as.character(x)) +
#           scale_y_continuous(expand = c(0, 0), labels = function(x) as.character(x)) +
#           # scale_color_gradientn(colors = colfunc(100)) +
#           labs(color = "Track", x = "Longitude", y = "Latitude"))
#   dev.off()
# };rm(i);rm(x)

# saveRDS(segmented_tracks, file = "/home/hbronnvik/Documents/chapter3/black_storks/segmented_tracks_250807.rds")

# now turn to circuit analysis
# we need sources and grounds for each track
####-----------------------------------------
# the tracks
segmented_tracks <- readRDS("/home/hbronnvik/Documents/chapter3/black_storks/segmented_tracks_250807.rds")

# their info
meta <- segmented_tracks %>%
  group_by(individual.id, season) %>%
  count(track_id) %>% 
  mutate(age = row_number()) %>% 
  ungroup() %>% 
  group_by(individual.id) %>%
  mutate(total_journeys = n()) %>% 
  ungroup()

# the first fall data
fall1s <- meta %>% 
  filter(season == "fall" & age == 1)

grounds <- segmented_tracks %>% 
  filter(track_id %in% fall1s$track_id) %>% 
  arrange(individual.id, timestamp) %>% 
  group_by(track_id) %>% 
  slice(n()) %>% 
  ungroup() %>% 
  # add in an additional filter that removes European "grounds" (probably error/death)
  filter(location.lat < 20) 

sources <- segmented_tracks %>%
  # use only starts corresponding to grounds in sub-Saharan Africa
  filter(track_id %in% grounds$track_id) %>% 
  arrange(individual.id, timestamp) %>% 
  group_by(track_id) %>% 
  slice(1) %>% 
  ungroup() %>% 
  dplyr::select(location.long, location.lat) %>% 
  sf::st_as_sf(coords = c('location.long', 'location.lat'), crs = 4326) %>% 
  sf::st_transform("ESRI:53009") %>%
  mutate(location.long = sf::st_coordinates(.)[,1],
         location.lat = sf::st_coordinates(.)[,2]) %>%
  sf::st_drop_geometry() %>%
  as.data.frame()

# write.table(sources, file = "/home/hbronnvik/Documents/chapter3/black_storks/sources_bs.txt",
#             col.names = F, quote = F, sep = "\t")

grounds <- grounds %>% 
  dplyr::select(location.long, location.lat) %>% 
  sf::st_as_sf(coords = c('location.long', 'location.lat'), crs = 4326) %>% 
  sf::st_transform("ESRI:53009") %>%
  mutate(location.long = sf::st_coordinates(.)[,1],
         location.lat = sf::st_coordinates(.)[,2]) %>%
  sf::st_drop_geometry() %>%
  as.data.frame()
# write.table(grounds, file = "/home/hbronnvik/Documents/chapter3/black_storks/grounds_bs.txt",
#             col.names = F, quote = F, sep = "\t")

# now we need to run Circuitscape in Julia on the HPC

# the current map from the HPC
curmap <- rast("~/Documents/chapter3/black_storks/curconn_res_blh_demBS_curmap.asc")
plot(sqrt(curmap))
names(curmap) <- "current"
crs(curmap) <- "ESRI:53009"

# find all the values of the current greater than 0
v1 <- as.numeric(values(curmap, na.rm = T))
v1 <- v1[v1>0]
ecdf_func <- ecdf(v1)

# annotate the tracks with the current
segmented_tracks <- segmented_tracks  %>% 
  filter(track_id %in% meta$track_id) %>% 
  sf::st_as_sf(coords = c('location.long', 'location.lat'), crs = 4326) %>% 
  st_transform("ESRI:53009") %>% 
  mutate(extract(curmap, .),
         percentile = ecdf_func(current)*100,
         location.long = sf::st_coordinates(.)[,1],
         location.lat = sf::st_coordinates(.)[,2])

segmented_tracks %>% 
  st_drop_geometry() %>% 
  summarize(Mean = mean(percentile),
            Std.Dev. = sd(percentile))

hist(segmented_tracks$percentile)

# a plot of the current
library(rnaturalearth)
world <- ne_coastline(scale = "large", returnclass = "sf")
curp <- ggplot() +
  tidyterra::geom_spatraster(data = sqrt(curmap)) +
  geom_sf(data = world, color = "gray50") +
  coord_sf(xlim = c(-1771942.75609039, 5005879.4100524), 
           ylim = c(-3639482.55008231, 6869064.04589582)) +
  scale_x_continuous(expand = c(0, 0), labels = function(x) as.character(x)) +
  scale_y_continuous(expand = c(0, 0), labels = function(x) as.character(x)) + 
  # tidyterra::scale_fill_grass_c(palette = "grey") +
  tidyterra::scale_fill_grass_c(palette = "haxby") +
  labs(x = "Longitude", y = "Latitude", fill = "sqrt current")
# and a comparison to the tracks
locp <- ggplot() +
  geom_sf(data = segmented_tracks, aes(color = sqrt(current)), color = "white") +
  geom_sf(data = world, color = "gray50") +
  geom_path(data = segmented_tracks, aes(location.long, location.lat, group = track_id), 
            color = "gray60", alpha = 0.5) +
  geom_sf(data = segmented_tracks, aes(color = sqrt(current)), 
          size = 0.5, alpha = 0.5) + 
  coord_sf(xlim = c(-1771942.75609039, 5005879.4100524),
           ylim = c(-3639482.55008231, 6869064.04589582)) +
  scale_x_continuous(expand = c(0, 0), labels = function(x) as.character(x)) +
  scale_y_continuous(expand = c(0, 0), labels = function(x) as.character(x)) + 
  tidyterra::scale_color_grass_c(palette = "haxby") +
  labs(x = "Longitude", y = "Latitude", color = "sqrt current")
# png(filename = "/home/hbronnvik/Documents/chapter3/black_storks/current.png",
#     height = 8, width = 11, units = "in", res = 500)
ggpubr::ggarrange(curp, locp, common.legend = T, legend = "right", labels = "AUTO")
# dev.off()
