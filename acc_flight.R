### ACC of stork flight
### Hester Bronnvik
### 2024-03-21
library(move)
library(moveACC)
library(tidyverse)
theme_set(theme_classic()+theme(axis.text = element_text(color = "black", size = 12), text = element_text(size = 15)))
colfunc <- colorRampPalette(c("#0081A7", "#0098b0", "#00AFB9", "#7fc4b8", "#FED9B7", "#f7a58f", "#f27e71", "#F07167", "#ee5e53"))
# colfunc <- colorRampPalette(c("#8A2846", "#B9375E", "#E05780", "#FF7AA2", "#FF9EBB", "#FFC2D4", "#FFE0E9"))
load("/home/hbronnvik/Documents/storkSSFs/loginStored.RData")
studies <- c(24442409, 212096177, 76367850, 21231406, 1176017658, 173641633, 908232414)

info <- lapply(studies, function(x){
  info <- getMovebankAnimals(x, loginStored) %>% 
    filter(sensor_type_id == 2365683 & grepl("acceleration", sensor_type_ids)) %>% 
    mutate(study = x)
  return(info)
}) %>% reduce(rbind)

records <- readRDS("/home/hbronnvik/Documents/chapter2/migration_dates_2024-01-15.rds") %>% 
  rowwise() %>% 
  mutate(date = as.Date(str_split(id_date, "_")[[1]][2]),
         season = ifelse(grepl("fall", trackID), "Fall", "Spring")) %>% 
  ungroup() %>% 
  group_by(trackID) %>% 
  mutate(ld_day = 1:length(unique(date))) %>% 
  ungroup() %>% 
  dplyr::select(individual.id, trackID, date, ld_day, season)

meta <- readRDS("/home/hbronnvik/Documents/chapter2/migration_dates_2024-01-15.rds") %>% 
  rowwise() %>% 
  mutate(date = as.Date(str_split(id_date, "_")[[1]][2]),
         season = ifelse(grepl("fall", trackID), "Fall", "Spring")) %>% 
  ungroup() %>% 
  group_by(trackID) %>% 
  mutate(ld_day = 1:length(unique(date))) %>% 
  ungroup() %>% 
  dplyr::select(individual.id, trackID, date, ld_day, season) %>% 
  group_by(individual.id, season) %>% 
  count(trackID) %>% 
  mutate(migration = row_number()) %>% 
  ungroup() %>% 
  group_by(individual.id) %>%
  mutate(total_journeys = n()) %>% 
  ungroup()

# GPS flight data
gps_files <- list.files("/home/hbronnvik/Documents/chapter2/classified_data", full.names = T)#[grepl("173673348|173661590|23460438", list.files("/home/hbronnvik/Documents/chapter2/classified_data", full.names = T))]
gps_birds <- lapply(gps_files, function(f){
  gsub("_classified_bursts.rds", "", 
       gsub("/home/hbronnvik/Documents/chapter2/classified_data/", "", f))
}) %>% unlist() %>% paste(collapse = "|")
# all the ACC data
acc_files <- list.files("/home/hbronnvik/Documents/chapter2/acc_data", full.names = T)
acc_birds <- lapply(acc_files, function(f){
  gsub("_acc_2023-10-29.rds", "", 
       gsub("/home/hbronnvik/Documents/chapter2/acc_data/", "", f))
}) %>% unlist() %>% paste(collapse = "|")

gps_files <- gps_files[grepl(acc_birds, gps_files)]

start_time <- Sys.time()
flight_acc <- lapply(1:length(gps_files), function(x){
  print(x)
  # load the GPS data
  gex <- readRDS(gps_files[x]) %>% 
    rename(burstID = burst_id)
  # get the dates of the migration
  rec <- records %>% 
    filter(individual.id == unique(gex$individual.id) & season == "Fall")
  
  if(nrow(rec) > 0){
    # filter the GPS data to migration and add a rounded-off minute
    gex <- lapply(split(rec, rec$trackID), function(track){
      ex <- gex %>% 
        sf::st_drop_geometry() %>% 
        mutate(date = date(timestamp),
               minute = round_date(timestamp, "minute"),
               trackID = unique(track$trackID)) %>% 
        filter(date %in% track$date)
      return(ex)
    }) %>% reduce(rbind)
    
    if(nrow(gex) > 0){
      # load the matching ACC data
      aex <- readRDS(acc_files[grepl(unique(gex$individual.id), acc_files)])
      # s_id <- info$study[info$individual.id == unique(gex$individual.id)]
      # aex <- getMovebankNonLocationData(study = s_id, animalName = unique(gex$individual.id),
      #                                   sensorID = 2365683, login = loginStored)
      # filter ACC data to the migrations
      aex <- lapply(split(rec, rec$trackID), function(track){
        ex <- aex %>% 
          mutate(trackID = unique(track$trackID)) %>% 
          filter(date(timestamp) %in% track$date)
        return(ex)
      }) %>% reduce(rbind)
      
      aex <- aex %>% 
        filter(round_date(timestamp, "minute") %in% unique(round_date(gex$timestamp, "minute")))
      
      # attach whether this is a classified flight burst to the ACC data
      aex <- lapply(1:nrow(aex), function(n){
        # print(n)
        acc <- aex[n,]
        # use the time from the ACC
        ts <- round_date(acc$timestamp, "minute")
        # extract the GPS burst(s) containing that minute
        burst_oi <- unique(gex$burstID[gex$minute == ts])
        # if there is a burst, then the ACC followed GPS classified flight, else it did not
        if(length(burst_oi) == 0){
          acc$flight <- F
          acc$location.long <- NA
          acc$location.lat <- NA
          acc$burstID <- NA
        }else{
          # take the GPS burst with the closest timestamp to the ACC
          burst_oi <- gex %>% 
            filter(burstID %in% burst_oi) %>% 
            dplyr::select(timestamp, burstID) %>% 
            mutate(prox = timestamp - acc$timestamp) %>% 
            filter(abs(prox) == min(abs(prox))) %>% 
            dplyr::select(burstID) %>% 
            deframe() %>% 
            unique()
          acc$flight <- T
          # the last location before the ACC burst started
          acc$location.long <- gex$location.long[gex$burstID == burst_oi][nrow(gex[gex$burstID == burst_oi,])]
          acc$location.lat <- gex$location.lat[gex$burstID == burst_oi][nrow(gex[gex$burstID == burst_oi,])]
          acc$burstID <- burst_oi
        }
        return(acc)
      }) %>% reduce(rbind)
      
      # keep the ACC data that were on migration and in flight
      acc_remaining <- aex %>% 
        filter(flight == T) %>% 
        mutate(date = date(timestamp)) %>% 
        left_join(rec[, c("date", "ld_day")], by = join_by(date)) %>% 
        left_join(meta[, c("trackID", "migration")], by = join_by(trackID)) 
      return(acc_remaining)
    }
  }
})
Sys.time()-start_time # Time difference of 19.37883 mins
flight_acc <- flight_acc[!sapply(flight_acc, function(x) is.null(x))]

# saveRDS(flight_acc, file = "/home/hbronnvik/Documents/chapter2/flight_acc_20240327.rds")

flight_acc <- readRDS(file = "/home/hbronnvik/Documents/chapter2/flight_acc_20240327.rds")


flight_DBA <- lapply(flight_acc, function(acc){
  VeDBA <- lapply(1:nrow(acc), function(x){
    # regardless of the name of the column in the data, use the column including raw data, 
    # the column including the sampling freqency, and the column including the axes (XYZ)
    accRawCol <- grep("accelerationTransformed", names(acc), value=T)
    sampFreqCol <- grep("acceleration_sampling_frequency_per_axis", names(acc), value=T)
    axesCol = grep("acceleration_axes", names(acc), value=T)
    
    # only proceed if the data have all three axes
    
    if(nchar(acc[x, axesCol])<3){stop("The ACC data have fewer than 3 axes.")}
    # row-by-row, turn the raw acc column into a matrix with separate columns for each axis
    accMx <- matrix(as.integer(unlist(strsplit(as.character(acc[x, accRawCol]), " "))), ncol=3, byrow = T)
    n_samples_per_axis <- nrow(accMx)
    acc_burst_duration_s <- n_samples_per_axis/acc[x, sampFreqCol]
    # calculate the vectorial dynamic body acceleration (the vectorial sum of the axes: mathematically correct)
    VeDBA <- sqrt((accMx[,1]-mean(accMx[,1]))^2 + (accMx[,2]-mean(accMx[,2]))^2 + (accMx[,3]-mean(accMx[,3]))^2)
    # calculate the overall dynamic body acceleration (the simple sum of the axes: not mathematically correct)
    # ODBA <- (accMx[,1]-mean(accMx[,1])) + (accMx[,2]-mean(accMx[,2])) + (accMx[,3]-mean(accMx[,3]))
    # return the average VeDBA from that one row
    DBA <- mean(VeDBA)
    return(DBA)
  }) %>% unlist()
  acc$VeDBA <- VeDBA
  return(acc)
}) %>% 
  reduce(rbind) %>% 
  group_by(individual_id) %>% 
  mutate(eastern = max(location.long >= 15)) %>% 
  ungroup() %>% 
  filter(eastern == F) 

flight_DBA <- flight_DBA %>% 
  group_by(individual_id) %>% 
  mutate(dba_res = VeDBA-mean(VeDBA)) %>% 
  ungroup()

flight_DBA %>% 
  group_by(location.long, location.lat) %>% 
  mutate(avg_DBA = mean(dba_res)) %>% 
  slice(1) %>% 
  ungroup() %>% 
  ggplot(aes(location.long, location.lat)) +
  borders(xlim = c(-20,27), ylim = c(10, 50), fill = "white") +
  geom_point(aes(color = avg_DBA), alpha = 0.5) +
  scale_colour_gradientn(colors = colfunc(130)) +
  labs(x = "Longitude", y= "Latitude", color = "VeDBA (g)")

flight_DBA %>% 
  filter(migration < 4) %>% 
  group_by(location.long, location.lat) %>% 
  mutate(avg_DBA = mean(VeDBA)) %>% 
  slice(1) %>% 
  ungroup() %>% 
  # left_join(meta[, c("trackID", "migration")], by = join_by(trackID)) %>% 
  ggplot(aes(as.factor(migration), dba_res, fill = migration)) +
  geom_boxplot() +
  labs(x = "Migration", y= "Mean burst VeDBA - mean tag VeDBA (g)", color = "Migration")

flight_DBA %>% 
  filter(migration < 4) %>%
  group_by(individual_id, date(timestamp)) %>% 
  mutate(daily_mean = mean(dba_res)) %>% 
  ungroup() %>% 
  group_by(trackID) %>% 
  mutate(sum_mean = sum(daily_mean)) %>% 
  ggplot(aes(migration, sum_mean)) +
  geom_point(alpha = 0.5, color = "#B379CD") +
  geom_smooth(method = "lm", color = "black") +
  scale_x_continuous(breaks = c(1, 2, 3)) +
  labs(x = "Migration", y= "Cumulative VeDBA", color = "Migration")

lines <- flight_DBA %>% 
  filter(migration < 4) %>%
  group_by(individual_id, date(timestamp)) %>% 
  mutate(daily_mean = mean(VeDBA)) %>% 
  ungroup() %>% 
  group_by(trackID) %>% 
  mutate(n = length(unique(date(timestamp))),
         sum_mean = sum(daily_mean)/n) %>% 
  ungroup() %>% 
  group_by(migration) %>% 
  summarize(m = median(sum_mean),
            n = length(unique(trackID)))

# png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/poster/data_vis/VeDBA_violins.png",
#     height = 8.5, width = 11, units = "in", res = 500)
flight_DBA %>% 
  filter(migration < 4) %>%
  group_by(individual_id, date(timestamp)) %>% 
  mutate(daily_mean = mean(VeDBA)) %>% 
  ungroup() %>% 
  group_by(trackID) %>% 
  mutate(n = length(unique(date(timestamp))),
         sum_mean = sum(daily_mean)/n) %>% 
  ungroup() %>% 
  ggplot(aes(as.factor(migration), sum_mean)) +
  geom_hline(yintercept = lines$m, color = c("#5F2C77", "#994CBD", "#C9A0DC"), lty = 2) +
  geom_violin(aes(fill = as.factor(migration))) +
  geom_boxplot(width = 0.05, alpha = 0.75) +
  # geom_jitter(width = 0.033, alpha = 0.01) +
  scale_fill_manual(values = c("#A25BC2", "#B379CD", "#C9A0DC")) +
  labs(x = "Migration", y= "Cumulative daily flight VeDBA", fill = "Migration") +
  annotate(geom = "text", x = c(1:3), y = -1, label = paste0("n = ", lines$n), size = 5)
# dev.off()

file <- "/home/hbronnvik/Documents/storkSSFs/ecmwf/single/boundary_layer_height_land_2022.nc"
blh <- try(raster::raster(file))
template <- blh[[1]]

ras_dba <- lapply(split(flight_DBA, flight_DBA$migration), function(m){
  ras_dba <- m
  sp::coordinates(ras_dba) <- ~location.long+location.lat
  sp::proj4string(ras_dba) <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs")
  ras_dba <- raster::rasterize(x = ras_dba, y = template, 
                               field = "dba_res", fun = mean)
  # raster::plot(ras_thermals)
  ras_dba_df <- ras_dba %>% 
    as.data.frame(xy = T) %>% 
    drop_na(layer) %>% 
    mutate(migrations = unique(m$migration))
}) %>% reduce(rbind)

# png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/poster/data_vis/VeDBA_raster.png",
#     height = 11.5, width = 8.5, units = "in", res = 500)
ras_dba %>% 
  ggplot(aes(x, y, fill = layer))+
  borders(xlim = c(-20,15), ylim = c(10, 50), fill = "white") +
  geom_raster() +
  scale_fill_gradientn(colors = colfunc(130),
                       values = scales::rescale(c(-0.4100939, 
                                                  0, 0.359632)),
                       n.breaks = 10) +
  scale_x_continuous(n.breaks = 13) +
  scale_y_continuous(n.breaks = 13) +
  labs(x = "Longitude", y = "Latitude", fill = "VeDBA")
# dev.off()

flight_DBA %>% 
  filter(migration < 4) %>% 
  mutate(lat_bin = ifelse(location.lat > 42.5, "Central Europe",
                          ifelse(between(location.lat, 36, 42.5), "Spain",
                                 ifelse(between(location.lat, 27.4, 36), "Morocco",
                                        ifelse(location.lat < 15, "Sahel", "Sahara")))),
         lat_bin = factor(lat_bin, levels = c("Central Europe", "Spain", "Morocco", "Sahara", "Sahel"))) %>% 
  ggplot(aes(dba_res, fill = as.factor(migration), y = as.factor(migration))) +
  ggridges::geom_density_ridges2(alpha = 0.5) +
  labs(x = "VeDBA", y = "Density", fill = "Migration") +
  scale_fill_manual(values = colfunc(3)) +
  facet_wrap(~lat_bin, nrow = 1)

# relate these DBAs to resting DBA to account for tag and attachment differences
ggplot(statsDF, aes(timestamp, odbaMedian)) +
  geom_segment(aes(xend=timestamp, y = 0, yend=odbaMedian)) + 
  geom_point(size=2, color = "#4F8FB4") + 
  geom_point(data = statsDF %>% filter(odbaMedian == min(odbaMedian)), size=2, color = "#CB8753") 

ggplot(statsDF, aes(timestamp, avgX)) +
  geom_line(color = "#4F8FB4") +
  geom_line(aes(timestamp, avgY), color = "#50A992") +
  geom_line(aes(timestamp, avgZ), color = "#CB8753")

early_DBA <- lapply(split(acc, acc$timestamp), function(f_acc){
  # regardless of the name of the column in the data, use the column including raw data, 
  # the column including the sampling frequency, and the column including the axes (XYZ)
  acc_column <- grep("accelerationTransformed", names(f_acc), value=T)
  sampling_freq <- grep("acceleration_sampling_frequency_per_axis", names(f_acc), value=T)
  axes_column <- grep("acceleration_axes", names(f_acc), value=T)
  
  dba <- lapply(1:nrow(f_acc), function(a){
    ind_dba <- data.frame(individual.id = unique(f_acc$individual_id),
                          timestamp = f_acc$timestamp[a],
                          # location.long = f_acc$location.long[a],
                          # location.lat = f_acc$location.lat[a],
                          # burstID = f_acc$burstID[a],
                          # trackID = f_acc$trackID[a],
                          dba = VeDBA(f_acc[a,], acc_column, sampling_freq, axes_column))
  }) %>% reduce(rbind)
  return(dba)
}) %>% reduce(rbind)

ggplot(early_DBA, aes(timestamp, dba)) +
  geom_segment(aes(xend=timestamp, y = 0, yend=dba)) + 
  geom_point(size=2, color = "#4F8FB4") + 
  geom_point(data = early_DBA %>% filter(dba == min(dba)), size=2, color = "#CB8753") 

flapping <- lapply(1:length(flight_acc), function(a){
  print(paste0("Processing transformed flight ACC data, individual ", a, " of ", length(flight_acc), "."), quote = F)
  acc <- flight_acc[[a]]
  if(nrow(acc) > 5){
    ID <- unique(acc$individual_id)
    waveDF <- ACCwave(acc, transformedData=T)
    # clusterPlot(waveDF, cluster=F)
    # clusterPlot(waveDF, cluster=T, forclustering= c("varWaveWingBeat","eigenValue1"))
    # wingBeatsPlot(dfw=waveDF, forclustering= c("varWaveWingBeat","eigenValue1"))
    # wingBeatsPlot(dfw=waveDF, forclustering= c("varWaveWingBeat","eigenValue1"), interactivePlot=F)
    # wingBeatsHist(dfw=waveDF, forclustering= c("varWaveWingBeat","eigenValue1"), interactivePlot=F)
    wingbeatsDF <- WingBeatsSelection(waveDF, forclustering= c("varWaveWingBeat","eigenValue1"), 
                                      minbeat=2, maxbeat=5)
    wingbeatsDF$individual_id <- ID
    # saveRDS(wingbeatsDF, file = paste0("/home/hbronnvik/Documents/chapter2/flapping_flight/", ID, "_270323.rds"))
    return(wingbeatsDF)
  }
}) %>% reduce(rbind)

build <- lapply(list.files("/home/hbronnvik/Documents/chapter2/flapping_flight", full.names = T), function(file){
  readRDS(file) %>% 
    rename(event_id = event.id)
}) %>% reduce(rbind)
buildDBA <- flight_DBA %>% 
  filter(individual_id %in% unique(build$individual_id)) %>% 
  left_join(build[, c("individual_id", "event_id", "behavior")], by = join_by(individual_id, event_id))


# png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/acc/VeDBA_flapping.png",
#     width = 11.5, height = 8.5, units = "in", res = 500)
ggplot(buildDBA, aes(behavior, dba_res, fill = behavior)) +
  geom_boxplot() +
  # geom_jitter(width = 0.25, alpha = 0.05) +
  scale_fill_manual(values = c("#334371", "#6279B8")) +
  labs(x = "Classification (Individuals = 136, n = 46706)", y = "Difference between VeDBA and its per-tag mean") +
  theme(legend.position = "none")
# dev.off()

# metad <- lapply(studies, function(x){
#   md <- getMovebankReferenceTable(study = x, login = loginStored) %>%
#     drop_na(animal_id) %>%
#     filter(sensor_type_id == 653) %>% 
#     dplyr::select(animal_id, deploy_on_timestamp) %>% 
#     mutate(deploy_on_timestamp = date(deploy_on_timestamp)) %>% 
#     rename(individual.id = animal_id)
#   return(md)
# }) %>% reduce(rbind)
# saveRDS(metad, file = "/home/hbronnvik/Documents/chapter2/deploy_on_times.rds")
metad <- readRDS("/home/hbronnvik/Documents/chapter2/deploy_on_times.rds")

# thermal data with wind vectors attached
classified <- lapply(list.files("/home/hbronnvik/Documents/chapter2/wind_thermal_data/", pattern = "2024-01-16|2024-01-17", full.names = T), readRDS)
classified <- classified[!sapply(classified, function(x) nrow(x)==0)]

classified <- lapply(classified, function(ind){
  ind <- ind %>% 
    mutate(time_diff = as.numeric(difftime(timestamp, lag(timestamp), units = "secs")),
           time_diff = ifelse(is.na(time_diff), 1, time_diff),
           min_split = time_diff > 60,
           thermal_event = paste(ind_burst_id, cumsum(min_split)+1, sep = "_")) %>% 
    dplyr::select(-time_diff, -min_split) %>% 
    group_by(thermal_event) %>% 
    mutate(thermal_duration = n(),
           vspeed_thermal = (height_above_ellipsoid[n()]-height_above_ellipsoid[1])/thermal_duration,
           turn_var_thermal = var(turn_angle, na.rm = T)/thermal_duration) %>% 
    ungroup()
  return(ind)
})

records <- readRDS("/home/hbronnvik/Documents/chapter2/migration_dates_2024-01-15.rds") %>% 
  mutate(season = ifelse(grepl("fall", trackID), "fall", "spring")) %>% 
  group_by(individual.id, season) %>% 
  count(trackID) %>%
  mutate(journey_number = row_number()) %>% 
  ungroup() %>% 
  group_by(individual.id) %>%
  mutate(total_journeys = n()) %>% 
  ungroup() %>% 
  left_join(metad) %>%
  dplyr::select(trackID, season, journey_number, deploy_on_timestamp)

m_days <- readRDS("/home/hbronnvik/Documents/chapter2/migration_dates_2024-01-15.rds") %>% 
  rowwise() %>% 
  mutate(date = as.Date(str_split(id_date, "_")[[1]][2]), 
         id_date = paste(individual.id, date, sep = "_")) %>% 
  ungroup() %>% 
  group_by(trackID) %>% 
  mutate(ld_day = 1:length(unique(date))) %>% 
  ungroup() %>% 
  dplyr::select(id_date, trackID, date, ld_day)

wind_support <- function(u,v,heading) {
  angle <- atan2(u,v) - heading/180*pi
  return(cos(angle) * sqrt(u*u+v*v))
}
cross_wind <- function(u,v,heading) {
  angle <- atan2(u,v) - heading/180*pi
  return(sin(angle) * sqrt(u*u+v*v))
}
wind_speed <- function(u,v) {
  return(sqrt(u*u+v*v))
}
wind_direction <- function(u, v){(90-(atan2(v, u)*(180/pi)))%%360 }

# winds
classified <- data.table::rbindlist(classified, use.names = T, fill = T) %>% 
  # filter(!individual.id %in% eastern_birds$individual.id) %>% 
  mutate(wind_speed = wind_speed(windX, windY),
         # north is 0 clockwise to 360 (the same as the heading from the tags)
         wind_direction = wind_direction(windX, windY),
         cross_wind = cross_wind(windX, windY, heading),
         wind_support = wind_support(windX, windY, heading),
         turn_direction = ifelse(turn.angle > 0, "counter", "clock"),
         directional_change = ifelse(turn_direction == lag(turn_direction), F, T),
         directional_set = cumsum(c(F, directional_change[2:n()]))) %>% 
  group_by(thermal_event) %>% 
  mutate(obs = 1:n(),
         position = as.numeric(scale(obs)),
         n_switches = sum(directional_change),
         exit_height = height_above_ellipsoid[n()]) %>% 
  ungroup() %>% 
  left_join(records, by = join_by(trackID)) %>% 
  left_join(m_days) %>% 
  dplyr::select(individual.id, timestamp, location.long, location.lat, journey_number, season, 
                deploy_on_timestamp, trackID, track_status, track_displacement, ld_day, burst_id,
                thermal_event, thermal_duration, vspeed_thermal, turn_var_thermal, heading, wind_speed,
                wind_direction, cross_wind, wind_support, n_switches, exit_height)
gc()

flapping <- flapping %>% 
  dplyr::select(-burstID) %>% 
  # get the burst id of the flight GPS associated with this ACC info
  left_join(flight_DBA[, c("timestamp", "individual_id", "VeDBA", "dba_res", "burstID")]) %>% 
  rename(burst_id = burstID,
         acc_stamp = timestamp)

flight_wind_acc <- classified %>% 
  group_by(burst_id) %>% 
  mutate(mean_speed = mean(wind_speed, na.rm = T),
         mean_support = mean(wind_support, na.rm = T),
         mean_cross = mean(cross_wind, na.rm = T)) %>% 
  slice(1) %>% 
  ungroup() %>% 
  left_join(flapping)

ggplot(flight_wind_acc %>% drop_na(behavior), aes(sqrt(mean_speed))) +
  geom_histogram(bins = 100, color = "black") +
  # labs() +
  facet_wrap(~journey_number) +
  theme(legend.position = "none")

png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/acc/flapping_wind.png",
    width = 11.5, height = 8.5, units = "in", res = 500)
flight_wind_acc %>% 
  filter(season == "fall") %>% 
  drop_na(behavior) %>% 
  mutate(migration = paste0("Fall migration ", journey_number)) %>% 
  ggplot(aes(behavior, sqrt(mean_speed), fill = behavior)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#334371", "#6279B8")) +
  labs(x = "", y = "sqrt mean wind speed in the GPS burst preceding flapping ACC (m/s)") +
  facet_wrap(~migration) +
  theme(legend.position = "none")
dev.off()

png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/acc/04_DBA_wind.png",
    width = 11.5, height = 8.5, units = "in", res = 500)
flight_wind_acc %>% 
  filter(season == "fall" & journey_number < 4) %>% 
  drop_na(behavior) %>% 
  mutate(migration = paste0("Migration ", journey_number)) %>% 
  ggplot(aes(mean_speed, dba_res)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", aes(group = migration, color = migration)) +
  scale_color_manual(values = c("#334371", "#6279B8", "#8E9FCC")) +
  labs(y = "centered VeDBA", 
       x = " Mean wind speed in the GPS burst preceding flapping ACC (m/s)", 
       color = "Migration")
dev.off()

png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/acc/06_flapping_vspeed.png",
    width = 11.5, height = 8.5, units = "in", res = 500)
flight_wind_acc %>% 
  filter(season == "fall" & journey_number < 4) %>% 
  drop_na(behavior) %>% 
  mutate(migration = paste0("Fall migration ", journey_number)) %>% 
  ggplot(aes(behavior, vspeed_thermal, fill = behavior)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#334371", "#6279B8")) +
  labs(x = "", y = "Vertical speed in thermals (m/s)") +
  facet_wrap(~migration) +
  theme(legend.position = "none")
dev.off()

png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/acc/07_flapping_switching.png",
    width = 11.5, height = 8.5, units = "in", res = 500)
flight_wind_acc %>% 
  filter(season == "fall" & journey_number < 4) %>% 
  drop_na(behavior) %>% 
  mutate(migration = paste0("Fall migration ", journey_number)) %>% 
  ggplot(aes(behavior, log(n_switches), fill = behavior)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#334371", "#6279B8")) +
  labs(x = "", y = "log number of directional switches in thermals") +
  facet_wrap(~migration) +
  theme(legend.position = "none")
dev.off()

png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/acc/08_flapping_cross.png",
    width = 11.5, height = 8.5, units = "in", res = 500)
flight_wind_acc %>% 
  filter(season == "fall" & journey_number < 4) %>% 
  drop_na(behavior) %>% 
  mutate(migration = paste0("Fall migration ", journey_number)) %>% 
  ggplot(aes(behavior, log(abs(mean_cross)), fill = behavior)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#334371", "#6279B8")) +
  labs(x = "", y = "log absolute cross wind speed (m/s)") +
  facet_wrap(~migration) +
  theme(legend.position = "none")
dev.off()

png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/acc/08_flapping_hae.png",
    width = 11.5, height = 8.5, units = "in", res = 500)
flight_wind_acc %>% 
  filter(season == "fall" & journey_number < 4) %>% 
  drop_na(behavior) %>% 
  mutate(migration = paste0("Fall migration ", journey_number)) %>% 
  ggplot(aes(behavior, sqrt(exit_height), fill = behavior)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#334371", "#6279B8")) +
  labs(x = "", y = "sqrt thermal exit height above the ellipsoid (m)") +
  facet_wrap(~migration) +
  theme(legend.position = "none")
dev.off()

png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/acc/06_VeDBA_vspeed.png",
    width = 11.5, height = 8.5, units = "in", res = 500)
flight_wind_acc %>% 
  filter(season == "fall" & journey_number < 4) %>% 
  drop_na(behavior) %>% 
  mutate(migration = paste0("Fall migration ", journey_number)) %>% 
  ggplot(aes(sqrt(vspeed_thermal), dba_res)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", aes(color = migration, group = migration)) +
  scale_color_manual(values = c("#334371", "#6279B8", "#8E9FCC")) +
  labs(y = "Centered VeDBA", x = "Vertical speed in thermals (m/s)", color = "Migration")
dev.off()

# hist(flight_wind_acc$exit_height, breaks = 20)
# library(MASS)
# b <- boxcox(lm(flight_wind_acc$exit_height ~ 1))
# # Exact lambda
# lambda <- b$x[which.max(b$y)]
# lambda
# new_x_exact <- (x ^ lambda - 1) / lambda

thermal_other <- readRDS("/home/hbronnvik/Documents/chapter2/between_thermals_20240325.rds")

thermal_count <- classified %>%
  mutate(date = date(timestamp)) %>% 
  group_by(individual.id, date) %>% 
  mutate(n_thermals = length(unique(thermal_event))) %>% 
  slice(1) %>% 
  ungroup() %>% 
  dplyr::select(date, individual.id, n_thermals)

thermal_other <- thermal_other  %>% 
  group_by(individual.id, date) %>% 
  mutate(n_falls = length(unique(paste0(burst_ind_ID, "_", exit)))) %>% 
  left_join(records, by = join_by(trackID)) %>% 
  left_join(thermal_count, by = join_by(individual.id, date)) %>% 
  mutate(fall_in_ratio = n_falls/n_thermals) %>% 
  filter(journey_number < 4)

hist(thermal_other$fall_in_ratio, breaks= 100, main = "", col = "gray40")
b <- MASS::boxcox(lm(thermal_other$fall_in_ratio ~ 1))
b$x[which.max(b$y)]
lambda <- b$x[which.max(b$y)]
thermal_other$fall_in_bc <- ((thermal_other$fall_in_ratio^lambda-1)/lambda)

thermal_other <- thermal_other %>% 
  filter(between(out_secs, 10, 60)) 

flaps <- flapping %>% 
  drop_na(burst_id) %>% 
  mutate(burst_ind_ID = paste(individual_id, burst_id, sep = "_"))
flaps <- flaps[!duplicated(flaps$burst_ind_ID),]
build <- thermal_other %>% 
  left_join(flaps[, c("burst_ind_ID", "behavior")]) %>% 
  mutate(migration = paste0("Migration ", journey_number))

build %>% 
  drop_na(behavior) %>% 
  ggplot(aes(behavior, )) +
  geom_boxplot()
build %>% 
  drop_na(behavior) %>% 
  ggplot(aes(behavior, exit_vspeed)) +
  geom_boxplot()
build %>% 
  drop_na(behavior) %>% 
  ggplot(aes(behavior, final_height)) +
  geom_boxplot()
build %>% 
  drop_na(behavior) %>% 
  filter(season == "fall") %>% 
  ggplot(aes(behavior, avg_wind_speed_pre)) +
  geom_boxplot() +
  facet_wrap(~migration)

ms <- build %>% 
  group_by(burst_ind_ID) %>% 
  slice(1) %>% 
  ungroup() %>% 
  drop_na(behavior) %>% 
  group_by(migration, behavior) %>% 
  summarize(m = median(fall_in_bc, na.rm = T)) 

png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/acc/flapping_fob.png",
    width = 11.5, height = 8.5, units = "in", res = 500)
build %>% 
  group_by(burst_ind_ID) %>% 
  slice(1) %>% 
  ungroup() %>% 
  drop_na(behavior) %>% 
  ggplot(aes(behavior, fall_in_bc, fill = behavior)) +
  geom_hline(data = ms, aes(yintercept = m, color = behavior), lty = 2, show.legend = T, lwd = 1.5)  + 
  geom_violin() +
  geom_boxplot(width = 0.2, alpha = .75, lwd = 1.5, fill = "white") + 
  scale_fill_manual(values = c("#334371", "#6279B8")) + 
  scale_color_manual(values = c("#334371", "#6279B8")) +
  scale_y_continuous(n.breaks = 7) +
  labs(x = "", y = "Ratio of falls to total thermals per day") +
  facet_wrap(~migration) +
  theme(legend.position = "none")
dev.off()

thermal_other <- readRDS("/home/hbronnvik/Documents/chapter2/between_thermals_20240325.rds")

thermal_other <- thermal_other %>% 
  filter(between(out_secs, 10, 60)) 

thermal_count <- classified %>%
  mutate(date = date(timestamp)) %>% 
  group_by(individual.id, date) %>% 
  mutate(n_thermals = length(unique(thermal_event))) %>% 
  slice(1) %>% 
  ungroup() %>% 
  dplyr::select(date, individual.id, n_thermals)

thermal_other <- thermal_other  %>% 
  group_by(individual.id, date) %>% 
  mutate(n_falls = length(unique(paste0(burst_ind_ID, "_", exit)))) %>% 
  left_join(records, by = join_by(trackID)) %>% 
  left_join(thermal_count, by = join_by(individual.id, date)) %>% 
  mutate(fall_in_ratio = n_falls/n_thermals) %>% 
  filter(journey_number < 4)

flapping_count <- flapping %>%
  drop_na(behavior) %>% 
  mutate(date = date(acc_stamp)) %>% 
  group_by(individual_id, date) %>% 
  mutate(n_flaps = sum(behavior == "Flapping"),
         n_not_flaps = sum(behavior == "NotFlapping"),
         flap_ratio = n_flaps/n()) %>% 
  slice(1) %>% 
  ungroup() %>% 
  dplyr::select(date, individual_id, flap_ratio) %>% 
  rename(individual.id = individual_id)

colfuncBright <- colorRampPalette(c("#FED9B7", "#f7a58f", "#f27e71", "#F07167", "#ee5e53", "#EB3F33", "#CC2014", "#95180F"))
# per-day fall:thermal ratio and per-day flap:behavior ratio
png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/acc/daily_flapping_fob.png",
    width = 11.5, height = 8.5/2, units = "in", res = 500)
thermal_other %>%
  group_by(individual.id, date) %>% 
  slice(1) %>% 
  left_join(flapping_count) %>% 
  filter(flap_ratio != 0) %>% 
  mutate(migration = paste0("Migration ", journey_number)) %>% 
  group_by(individual.id) %>% 
  mutate(mort = ifelse(max(journey_number)==1, "Died on 1",
                       ifelse(max(journey_number) == 2, "Died on 2", "Lived to 3"))) %>% 
  ggplot(aes(log(fall_in_ratio), log(flap_ratio))) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(aes(color = n_thermals), cex = .5) +
  scale_color_gradientn(colors = colfuncBright(200)) +
  geom_smooth(method = "lm", color = "black", aes(lty = migration)) +
  coord_fixed() +
  labs(x = "log falling to thermals ratio (per day)", y = "log flapping to \nbursts ratio (per day)", 
       color = "Sample size", lty = "") +
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 10)) 
dev.off()


model_data <- flight_wind_acc %>% 
  filter(season == "fall" & journey_number < 4) %>% 
  mutate(individual.id = as.factor(individual.id),
         migration = as.factor(journey_number),
         sqrt_vspeed = sqrt(vspeed_thermal),
         log_turn_var_thermal = log(turn_var_thermal),
         behavior = ifelse(behavior == "Flapping", 1, 0)) %>% 
  dplyr::select(timestamp, location.long, location.lat, individual.id, ld_day, 
                migration, journey_number, trackID, thermal_event, track_status, season, 
                sqrt_vspeed, log_turn_var_thermal, mean_speed, dba_res, behavior) %>% 
  # left_join(turn_df, by = join_by(thermal_event)) %>% 
  mutate_at(c("location.lat","ld_day", "sqrt_vspeed", "log_turn_var_thermal", 
              "mean_speed", "journey_number", "dba_res"), list(z = ~(scale(.)))) %>% 
  mutate(thermal_event = as.factor(thermal_event))

library(lme4)

# have skill respond to age et al.
migrations_mod <- lmer(dba_res_z ~ journey_number_z*mean_speed_z*location.lat_z + (1 | individual.id), 
                       data = model_data, REML = T)

summary(migrations_mod)
cis <- confint(migrations_mod) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "pred") %>%
  left_join(fixef(migrations_mod) %>% 
              as.data.frame() %>% 
              rename(fixef = ".") %>% 
              rownames_to_column(var = "pred"), 
            by = join_by(pred))
# png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/acc/10_VeDBA_windAge.png",
#     width = 11.5, height = 8.5, units = "in", res = 500)
cis %>% 
  slice(4:n()) %>% 
  rename(lower = "2.5 %",
         upper = "97.5 %") %>% 
  # filter(pred %in% c("migration1", "migration2", "migration3", "avg_wind_speed_z", "location.lat_z")) %>% 
  mutate(pred = sub("journey_number_z", "Age", 
                    sub("mean_speed_z", "Wind speed",
                        sub("location.lat_z", "Latitude", pred)))) %>%
  ggplot(aes(fixef, pred)) +
  geom_vline(xintercept = 0, lty = 2, color = "gray40", lwd = 2) +
  geom_pointrange(aes(xmin = lower, xmax = upper), color = "#007BA7", fatten = 8, lwd = 2) +
  labs(x = "Effect on VeDBA +/- se", y = "Predictor")
# dev.off()

model_data %>% 
  ggplot(aes(as.factor(journey_number), dba_res)) +
  geom_boxplot()
# 
# model_data <- model_data %>% 
#   filter(track_status == "complete")

migrations_mod <- lmer(dba_res_z ~ journey_number_z*sqrt_vspeed_z*location.lat_z + (1 | individual.id), 
                       data = model_data, REML = T)

summary(migrations_mod)
cis <- confint(migrations_mod) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "pred") %>%
  left_join(fixef(migrations_mod) %>% 
              as.data.frame() %>% 
              rename(fixef = ".") %>% 
              rownames_to_column(var = "pred"), 
            by = join_by(pred))
yr_coefs <- cis %>% 
  slice(4:n()) %>% 
  rename(lower = "2.5 %",
         upper = "97.5 %") %>% 
  # filter(pred %in% c("migration1", "migration2", "migration3", "avg_wind_speed_z", "location.lat_z")) %>% 
  mutate(pred = sub("journey_number_z", "Age", 
                    sub("sqrt_vspeed_z", "Vertical speed",
                        sub("location.lat_z", "Latitude", pred)))) %>%
  ggplot(aes(fixef, pred)) +
  geom_vline(xintercept = 0, lty = 2, color = "gray40", lwd = 2) +
  geom_pointrange(aes(xmin = lower, xmax = upper), color = "#007BA7", fatten = 8, lwd = 2) +
  labs(x = "Effect on VeDBA +/- se", y = "Predictor")

grd_yr <- expand.grid(journey_number = 1:3,
                      sqrt_vspeed = seq(from = quantile(model_data$sqrt_vspeed, 0.025, na.rm = T), 
                                        to = quantile(model_data$sqrt_vspeed, 0.975, na.rm = T), 
                                        length.out = 100)) %>%
  mutate(location.lat = mean(model_data$location.lat, na.rm = T),
         # dba_res = mean(model_data$dba_res, na.rm = T),
         interaction = "age_climb")

set.seed(770)
n <- nrow(grd_yr)

new_data_only <- model_data %>%
  slice_sample(n = n, replace = F) %>% #randomly select n of these strata. Make sure n is less than n_distinct(data$stratum)
  # only keep the columns that I need
  dplyr::select(c("individual.id", "dba_res")) %>% 
  bind_cols(grd_yr) %>% 
  #calculate z-scores. get the mean(center) and sd(scale) from the previous z transformation. for consistency
  mutate(#log_turn_var_thermal_z = (log_turn_var_thermal - mean(model_data$log_turn_var_thermal))/(sd(model_data$log_turn_var_thermal)),
    # boxcox_n.changes_z = (boxcox_n.changes - mean(model_data$boxcox_n.changes))/(sd(model_data$boxcox_n.changes)),
    sqrt_vspeed_z = (sqrt_vspeed - mean(model_data$sqrt_vspeed, na.rm = T))/(sd(model_data$sqrt_vspeed, na.rm = T)),
    journey_number_z = (journey_number - mean(model_data$journey_number))/(sd(model_data$journey_number)),
    location.lat_z = (location.lat - mean(model_data$location.lat, na.rm = T))/(sd(model_data$location.lat, na.rm = T)),
    # avg_wind_speed_z = (avg_wind_speed - mean(model_data$avg_wind_speed, na.rm = T))/(sd(model_data$avg_wind_speed, na.rm = T))
    dba_res_z = (dba_res - mean(model_data$dba_res, na.rm = T))/(sd(model_data$dba_res, na.rm = T)))

new_data <- model_data %>% 
  drop_na(sqrt_vspeed) %>% 
  drop_na(dba_res) %>% 
  mutate(interaction = "OG_data") %>% 
  dplyr::select(names(new_data_only)) %>%  #only keep the columns that are necessary for the model
  bind_rows(new_data_only) 

# now that we have the values to predict, run the model on them
preds <- predict(migrations_mod, newdata = new_data, type = "link")

preds_pr <- new_data %>% 
  mutate(preds = preds) %>% 
  rowwise() %>% 
  mutate(probs = gtools::inv.logit(preds)) #https://rpubs.com/crossxwill/logistic-poisson-prob

inter_preds <- preds_pr %>% 
  filter(interaction != "OG_data") 

yr_tile <- inter_preds %>%
  filter(interaction == "age_climb") %>% 
  ggplot(aes(journey_number, sqrt_vspeed, fill = preds)) +
  geom_raster(interpolate = F) +
  scale_y_continuous(expand = c(0, 0), n.breaks = 8) +
  labs(x = "Migrations", y = "Vertical speed (m/s)", title = " ") +
  scale_fill_gradientn("Prediction", colors = colfunc(200),
                       values = scales::rescale(c(min(inter_preds$preds), 
                                                  0, max(inter_preds$preds))))

png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/acc/11_VeDBA_vertAge.png",
    width = 11.5, height = 8.5, units = "in", res = 500)
ggpubr::ggarrange(yr_coefs, yr_tile + theme(legend.position = "right"), 
                  nrow = 1)
dev.off()

# png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/acc/11_VeDBA_vertAge.png",
#     width = 11.5, height = 8.5, units = "in", res = 500)
# tiles <- ggpubr::ggarrange(yr_tile, yr_tile2+scale_x_reverse(n.breaks = 5), common.legend = T, legend = "right")
# ggpubr::ggarrange(yr_coefs, tiles, 
#                   nrow = 2)
# dev.off()

log_model <- lme4::glmer(behavior~journey_number_z*location.lat_z*mean_speed_z+(1|individual.id), data = model_data, family = "binomial")
summary(log_model)
log_graph <- summary(log_model)[[10]] %>% 
  as.data.frame() %>% 
  cbind(confint(log_model) %>% as.data.frame() %>% slice(2:n()))
colnames(log_graph) <- c("estimate", "std.err", "z_value", "p_value", "lower", "upper")
log_graph

aw_coefs <- log_graph %>% 
  rownames_to_column(var = "predictor") %>%
  filter(predictor != "(Intercept)") %>% 
  mutate(predictor = gsub("journey_number_z", "Migrations",
                          gsub("mean_speed_z", "Wind speed",
                               gsub("location.lat_z", "Latitude", predictor))),
         predictor = factor(predictor, levels = c("Latitude", "Wind speed", "Migrations", 
                                                  "Migrations:Latitude", "Migrations:Wind speed", "Latitude:Wind speed",
                                                  "Migrations:Latitude:Wind speed"))) %>%
  ggplot(aes(estimate, forcats::fct_rev(predictor))) +
  geom_vline(xintercept = 0, lty = 2, linewidth = 1) +
  geom_pointrange(aes(xmin = lower, xmax = upper), color = "#0081A7", linewidth = 1, size = 1) +
  labs(y = "", x = "Effect on flapping") 

grd_3 <- expand.grid(x = 1:3,
                     y = seq(from = quantile(model_data$mean_speed, 0.025, na.rm = T), to = quantile(model_data$mean_speed, 0.975, na.rm = T), length.out = 100),
                     z = seq(from = quantile(model_data$location.lat, 0.025, na.rm = T), to = quantile(model_data$location.lat, 0.975, na.rm = T), length.out = 100)) %>% # quantile(prep_d$w_star, .9, na.rm = T)
  rename(journey_number = x,
         mean_speed = y,
         location.lat = z) %>% 
  mutate(interaction = "map_age_wind")

grd_2 <- expand.grid(x = 1:3,
                     z = seq(from = quantile(model_data$mean_speed, 0.025, na.rm = T), to = quantile(model_data$mean_speed, 0.975, na.rm = T), length.out = 100)) %>% # quantile(prep_d$w_star, .9, na.rm = T)
  rename(journey_number = x,
         mean_speed = z) %>% 
  mutate(location.lat = mean(model_data$location.lat, na.rm = T),
         interaction = "age_wind")

set.seed(770)
n <- nrow(grd_2)

new_data_only <- model_data %>%
  drop_na(dba_res) %>% 
  slice_sample(n = n, replace = F) %>% #randomly select n of these strata. Make sure n is less than n_distinct(data$stratum)
  # only keep the columns that I need
  dplyr::select(c("individual.id", "behavior")) %>% 
  bind_cols(grd_2) %>% 
  #calculate z-scores. get the mean(center) and sd(scale) from the previous z transformation. for consistency
  mutate(#log_turn_var_thermal_z = (log_turn_var_thermal - mean(model_data$log_turn_var_thermal))/(sd(model_data$log_turn_var_thermal)),
    # boxcox_n.changes_z = (boxcox_n.changes - mean(model_data$boxcox_n.changes))/(sd(model_data$boxcox_n.changes)),
    mean_speed_z = (mean_speed - mean(model_data$mean_speed, na.rm = T))/(sd(model_data$mean_speed, na.rm = T)),
    journey_number_z = (journey_number - mean(model_data$journey_number))/(sd(model_data$journey_number)),
    location.lat_z = (location.lat - mean(model_data$location.lat, na.rm = T))/(sd(model_data$location.lat, na.rm = T)),
    # avg_wind_speed_z = (avg_wind_speed - mean(model_data$avg_wind_speed, na.rm = T))/(sd(model_data$avg_wind_speed, na.rm = T))
    # dba_res_z = (dba_res - mean(model_data$dba_res, na.rm = T))/(sd(model_data$dba_res, na.rm = T))
  )

new_data <- model_data %>% 
  # drop_na(sqrt_vspeed) %>% 
  # drop_na(dba_res) %>%
  drop_na(mean_speed) %>%
  drop_na(behavior) %>%
  mutate(interaction = "OG_data") %>% 
  dplyr::select(names(new_data_only)) %>%  #only keep the columns that are necessary for the model
  bind_rows(new_data_only) 

# now that we have the values to predict, run the model on them
preds <- predict(log_model, newdata = new_data, type = "link")

preds_pr <- new_data %>% 
  mutate(preds = preds) %>% 
  rowwise() %>% 
  mutate(probs = gtools::inv.logit(preds)) #https://rpubs.com/crossxwill/logistic-poisson-prob

inter_preds <- preds_pr %>% 
  filter(interaction != "OG_data") 

aw_tile3 <- inter_preds %>%
  filter(interaction == "map_age_wind") %>% 
  mutate(migration = paste0("Migration ", journey_number)) %>% 
  ggplot(aes(location.lat, mean_speed, fill = probs)) +
  geom_raster(interpolate = F) +
  scale_x_reverse(n.breaks = 7) +
  scale_y_continuous(expand = c(0, 0), n.breaks = 8) +
  labs(x = "Latitude", y = "Wind speed (m/s)", title = " ") +
  scale_fill_gradientn("Probability", colors = colfunc(200), 
                       n.breaks = 7,
                       # values = scales::rescale(c(min(inter_preds$probs),
                       #                            mean(inter_preds$probs),
                       #                            max(inter_preds$probs)))
  ) +
  facet_wrap(~migration)

aw_tile2 <- inter_preds %>%
  filter(interaction == "age_wind") %>% 
  mutate(migration = paste0("Migration ", journey_number)) %>% 
  ggplot(aes(journey_number, mean_speed, fill = probs)) +
  geom_raster(interpolate = F) +
  scale_y_continuous(expand = c(0, 0), n.breaks = 8) +
  labs(x = "Migration", y = "Wind speed (m/s)", title = " ") +
  scale_fill_gradientn("Probability", colors = colfunc(200), 
                       n.breaks = 6,
                       # values = scales::rescale(c(min(inter_preds$probs),
                       #                            mean(inter_preds$probs),
                       #                            max(inter_preds$probs)))
  )

png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/acc/12_flapping_windAge.png",
    width = 11.5, height = 8.5/2, units = "in", res = 500)
ggpubr::ggarrange(aw_coefs, aw_tile2 + theme(legend.position = "right"), 
                  nrow = 1, widths = c(2, 1), labels = c("A", "B"))
dev.off()

# model fall out and flapping
fob_mod_dat <- thermal_other %>%
  drop_na(avg_wind_speed_pre) %>% 
  group_by(individual.id, date) %>% 
  slice(1) %>% 
  left_join(flapping_count) %>% 
  filter(flap_ratio != 0) %>% 
  mutate(log_flaps = log(flap_ratio),
         log_falls = log(fall_in_ratio),
         yr_day = yday(date),
         latterdate = as.numeric(paste0(location.lat, as.numeric(date)))) %>% 
  ungroup()

fob_mod_dat$log_flaps_z <- scale((fob_mod_dat$log_flaps))[,1]
fob_mod_dat$log_falls_z <- scale((fob_mod_dat$log_falls))[,1]
fob_mod_dat$location.lat_z <- scale((fob_mod_dat$location.lat))[,1]
fob_mod_dat$journey_number_z <- scale((fob_mod_dat$journey_number))[,1]
fob_mod_dat$yr_day_z <- scale((fob_mod_dat$yr_day))[,1]
fob_mod_dat$avg_wind_speed_pre_z <- scale((fob_mod_dat$avg_wind_speed_pre))[,1]

lapply(split(fob_mod_dat, fob_mod_dat$journey_number), function(falls){
  
  fob_mod <- lmer(log_flaps_z ~ avg_wind_speed_pre_z*log_falls_z*location.lat_z + (1 | individual.id), 
                  data = falls, REML = T)
  # summary(fob_mod)
  cis <- confint(fob_mod) %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "pred") %>%
    left_join(fixef(fob_mod) %>% 
                as.data.frame() %>% 
                rename(fixef = ".") %>% 
                rownames_to_column(var = "pred"), 
              by = join_by(pred))
  # png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/acc/10_VeDBA_windAge.png",
  #     width = 11.5, height = 8.5, units = "in", res = 500)
  coef_plot <- cis %>% 
    slice(4:n()) %>% 
    rename(lower = "2.5 %",
           upper = "97.5 %") %>% 
    # filter(pred %in% c("migration1", "migration2", "migration3", "avg_wind_speed_z", "location.lat_z")) %>% 
    mutate(pred = sub("journey_number_z", "Age", 
                      sub("log_falls_z", "Fall ratio",
                          sub("location.lat_z", "Latitude", 
                              sub("avg_wind_speed_pre_z", "Wind speed", pred))))) %>%
    ggplot(aes(fixef, pred)) +
    geom_vline(xintercept = 0, lty = 2, color = "gray40", lwd = 2) +
    geom_pointrange(aes(xmin = lower, xmax = upper), color = "#007BA7", fatten = 8, lwd = 2) +
    labs(x = "Effect on flapping", y = "Predictor")
  
  grd_3 <- expand.grid(x = seq(from = quantile(falls$avg_wind_speed_pre, 0.025, na.rm = T), to = quantile(falls$avg_wind_speed_pre, 0.975, na.rm = T), length.out = 20),
                       y = seq(from = quantile(falls$log_falls, 0.025, na.rm = T), to = quantile(falls$log_falls, 0.975, na.rm = T), length.out = 20),
                       z = seq(from = quantile(falls$location.lat, 0.025, na.rm = T), to = quantile(falls$location.lat, 0.975, na.rm = T), length.out = 20)) %>% # quantile(prep_d$w_star, .9, na.rm = T)
    rename(avg_wind_speed_pre = x,
           log_falls = y,
           location.lat = z) %>% 
    mutate(interaction = "map_date_fall")
  
  grd_fall_lat <- expand.grid(x = seq(from = quantile(falls$log_falls, 0.025, na.rm = T), 
                                      to = quantile(falls$log_falls, 0.975, na.rm = T), 
                                      length.out = 20),
                              z = seq(from = quantile(falls$location.lat, 0.025, na.rm = T), 
                                      to = quantile(falls$location.lat, 0.975, na.rm = T), 
                                      length.out = 20)) %>% 
    rename(log_falls = x,
           location.lat = z) %>% 
    mutate(avg_wind_speed_pre = mean(falls$avg_wind_speed_pre, na.rm = T),
           interaction = "map_fall")
  
  grd_fall_date <- expand.grid(x = seq(from = quantile(falls$log_falls, 0.025, na.rm = T), 
                                       to = quantile(falls$log_falls, 0.975, na.rm = T), 
                                       length.out = 20),
                               z = seq(from = quantile(falls$avg_wind_speed_pre, 0.025, na.rm = T), 
                                       to = quantile(falls$avg_wind_speed_pre, 0.975, na.rm = T), 
                                       length.out = 20)) %>% 
    rename(log_falls = x,
           avg_wind_speed_pre = z) %>% 
    mutate(location.lat = mean(falls$location.lat, na.rm = T),
           interaction = "date_fall")
  
  grd_map_date <- expand.grid(x = seq(from = quantile(falls$location.lat, 0.025, na.rm = T), 
                                      to = quantile(falls$location.lat, 0.975, na.rm = T), 
                                      length.out = 20),
                              z = seq(from = quantile(falls$avg_wind_speed_pre, 0.025, na.rm = T), 
                                      to = quantile(falls$avg_wind_speed_pre, 0.975, na.rm = T), 
                                      length.out = 20)) %>% 
    rename(location.lat = x,
           avg_wind_speed_pre = z) %>% 
    mutate(log_falls = mean(falls$log_falls, na.rm = T),
           interaction = "date_map")
  
  grd_all <- bind_rows(grd_fall_date, grd_fall_lat, grd_map_date, grd_3) 
  
  set.seed(770)
  n <- nrow(grd_all)
  
  new_data_only <- falls %>%
    # group_by(individual.id) %>% 
    # slice_sample(n = 1) %>% #randomly selects one row (from each stratum)
    # ungroup() %>% 
    slice_sample(n = n, replace = T) %>% #randomly select n of these strata. Make sure n is less than n_distinct(data$stratum)
    # only keep the columns that I need
    dplyr::select(c("individual.id")) %>% 
    bind_cols(grd_all) %>% 
    #calculate z-scores. get the mean(center) and sd(scale) from the previous z transformation. for consistency
    mutate(log_falls_z = (log_falls - mean(fob_mod_dat$log_falls))/(sd(fob_mod_dat$log_falls)),
           avg_wind_speed_pre_z = (avg_wind_speed_pre - mean(fob_mod_dat$avg_wind_speed_pre))/(sd(fob_mod_dat$avg_wind_speed_pre)),
           location.lat_z = (location.lat - mean(fob_mod_dat$location.lat))/(sd(fob_mod_dat$location.lat)))
  
  new_data <- falls %>% 
    mutate(interaction = "OG_data") %>% 
    dplyr::select(names(new_data_only)) %>%  #only keep the columns that are necessary for the model
    bind_rows(new_data_only) 
  
  # now that we have the values to predict, run the model on them
  preds <- predict(fob_mod, newdata = new_data, type = "link")
  
  preds_pr <- new_data %>% 
    mutate(preds = preds) %>% 
    rowwise() %>% 
    mutate(probs = gtools::inv.logit(preds)) #https://rpubs.com/crossxwill/logistic-poisson-prob
  
  inter_preds <- preds_pr %>% 
    filter(interaction != "OG_data") 
  
  aw_tile3 <- inter_preds %>%
    filter(interaction == "map_date_fall") %>% 
    mutate(#avg_wind_speed_pre = sub("1970-", "", as.Date(avg_wind_speed_pre)),
      location.lat = as.factor(round(location.lat))) %>% 
    # mutate(migration = paste0("Migration ", journey_number)) %>% 
    ggplot(aes(avg_wind_speed_pre, log_falls, fill = probs)) +
    geom_tile() +
    geom_hline(lty = 2, yintercept = 0) +
    # scale_x_continuous(n.breaks = 7) +
    # scale_y_continuous(expand = c(0, 0), n.breaks = 8) +
    # labs(x = "Latitude", y = "Wind speed (m/s)", title = " ") +
    scale_fill_gradientn("Probability", colors = colfunc(200), 
                         n.breaks = 5,
                         values = scales::rescale(c(min(inter_preds$probs),
                                                    0.5,
                                                    max(inter_preds$probs)))
    ) + 
    facet_wrap(~fct_rev(location.lat)) +
    theme(axis.text = element_text(size = 12))
  
})





set.seed(770)
n <- nrow(grd_2)

new_data_only <- model_data %>%
  drop_na(dba_res) %>% 
  slice_sample(n = n, replace = F) %>% #randomly select n of these strata. Make sure n is less than n_distinct(data$stratum)
  # only keep the columns that I need
  dplyr::select(c("individual.id", "behavior")) %>% 
  bind_cols(grd_2) %>% 
  #calculate z-scores. get the mean(center) and sd(scale) from the previous z transformation. for consistency
  mutate(#log_turn_var_thermal_z = (log_turn_var_thermal - mean(model_data$log_turn_var_thermal))/(sd(model_data$log_turn_var_thermal)),
    # boxcox_n.changes_z = (boxcox_n.changes - mean(model_data$boxcox_n.changes))/(sd(model_data$boxcox_n.changes)),
    mean_speed_z = (mean_speed - mean(model_data$mean_speed, na.rm = T))/(sd(model_data$mean_speed, na.rm = T)),
    journey_number_z = (journey_number - mean(model_data$journey_number))/(sd(model_data$journey_number)),
    location.lat_z = (location.lat - mean(model_data$location.lat, na.rm = T))/(sd(model_data$location.lat, na.rm = T)),
    # avg_wind_speed_z = (avg_wind_speed - mean(model_data$avg_wind_speed, na.rm = T))/(sd(model_data$avg_wind_speed, na.rm = T))
    # dba_res_z = (dba_res - mean(model_data$dba_res, na.rm = T))/(sd(model_data$dba_res, na.rm = T))
  )



aw_tile2 <- inter_preds %>%
  filter(interaction == "age_wind") %>% 
  mutate(migration = paste0("Migration ", journey_number)) %>% 
  ggplot(aes(journey_number, mean_speed, fill = probs)) +
  geom_raster(interpolate = F) +
  scale_y_continuous(expand = c(0, 0), n.breaks = 8) +
  labs(x = "Migration", y = "Wind speed (m/s)", title = " ") +
  scale_fill_gradientn("Probability", colors = colfunc(200), 
                       n.breaks = 6,
                       # values = scales::rescale(c(min(inter_preds$probs),
                       #                            mean(inter_preds$probs),
                       #                            max(inter_preds$probs)))
  )



thermal_other <- thermal_other %>% 
  filter(between(out_secs, 10, 60)) 

flaps <- flapping %>% 
  drop_na(burst_id) %>% 
  mutate(burst_ind_ID = paste(individual_id, burst_id, sep = "_"))
flaps <- flaps[!duplicated(flaps$burst_ind_ID),]
build <- thermal_other %>% 
  left_join(flaps[, c("burst_ind_ID", "behavior", "acc_stamp")]) %>% 
  mutate(migration = paste0("Migration ", journey_number)) %>% 
  rename(timestamp = acc_stamp)

build %>% 
  drop_na(behavior) %>% 
  ggplot(aes(behavior, )) +
  geom_boxplot()
build %>% 
  drop_na(behavior) %>% 
  ggplot(aes(behavior, exit_vspeed)) +
  geom_boxplot()
build %>% 
  drop_na(behavior) %>% 
  ggplot(aes(behavior, final_height)) +
  geom_boxplot()
build %>% 
  drop_na(behavior) %>% 
  filter(season == "fall") %>% 
  ggplot(aes(behavior, avg_wind_speed_pre)) +
  geom_boxplot() +
  facet_wrap(~migration)

model_data <- build %>% 
  filter(season == "fall") %>% 
  group_by(location.lat) %>% 
  slice(1) %>% 
  ungroup() %>% 
  filter(season == "fall" & journey_number < 4) %>% 
  mutate(individual.id = as.factor(individual.id),
         migration = as.factor(journey_number),
         # sqrt_vspeed = sqrt(vspeed_thermal),
         # log_turn_var_thermal = log(turn_var_thermal),
         sqrt_avg_pre = sqrt(avg_wind_speed_pre),
         behavior = ifelse(behavior == "Flapping", 1, 0)) %>% 
  dplyr::select(timestamp, location.lat, individual.id, ld_day, 
                migration, journey_number, trackID, season, 
                sqrt_avg_pre, behavior) %>% 
  # drop_na(behavior) %>% 
  # drop_na(sqrt_avg_pre) %>% 
  # left_join(turn_df, by = join_by(thermal_event)) %>% 
  mutate_at(c("location.lat","ld_day", "sqrt_avg_pre", "journey_number"), list(z = ~(scale(.))))

library(lme4)

# have skill respond to age et al.
migrations_mod <- glmer(behavior ~ journey_number_z*sqrt_avg_pre_z*location.lat_z + (1|individual.id), 
                        data = model_data, family = "binomial")

summary(migrations_mod)
cis <- confint(migrations_mod, method = "Wald") %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "pred") %>%
  left_join(fixef(migrations_mod) %>% 
              as.data.frame() %>% 
              rename(fixef = ".") %>% 
              rownames_to_column(var = "pred"), 
            by = join_by(pred))
# png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/acc/13_flaps_windAge.png",
#     width = 11.5, height = 8.5, units = "in", res = 500)
flap_coefs <- cis %>% 
  slice(4:n()) %>% 
  rename(lower = "2.5 %",
         upper = "97.5 %") %>% 
  # filter(pred %in% c("migration1", "migration2", "migration3", "avg_wind_speed_z", "location.lat_z")) %>% 
  mutate(pred = sub("journey_number_z", "Age", 
                    sub("sqrt_avg_pre_z", "Wind speed",
                        sub("location.lat_z", "Latitude", pred)))) %>%
  ggplot(aes(fixef, pred)) +
  geom_vline(xintercept = 0, lty = 2, color = "gray40", lwd = 2) +
  geom_pointrange(aes(xmin = lower, xmax = upper), color = "#007BA7", fatten = 8, lwd = 2) +
  labs(x = "Effect on flapping (95% CI)", y = "Predictor", title = "Mean wind speed before a fall")
# dev.off()

grd_age_wind <- expand.grid(x = 1:3,
                            z = seq(from = quantile(model_data$sqrt_avg_pre, 0.025, na.rm = T), 
                                    to = quantile(model_data$sqrt_avg_pre, 0.975, na.rm = T), 
                                    length.out = 100)) %>% 
  rename(journey_number = x,
         sqrt_avg_pre = z) %>% 
  mutate(location.lat = mean(model_data$location.lat, na.rm = T),
         interaction = "age_wind")

# grd_all <- bind_rows(grd_fall_date, grd_fall_lat, grd_map_date, grd_3) 

set.seed(770)
n <- nrow(grd_age_wind)

new_data_only <- model_data %>%
  ungroup() %>% 
  # group_by(individual.id) %>% 
  # slice_sample(n = 1) %>% #randomly selects one row (from each stratum)
  # ungroup() %>% 
  slice_sample(n = n, replace = T) %>% #randomly select n of these strata. Make sure n is less than n_distinct(data$stratum)
  # only keep the columns that I need
  dplyr::select(c("individual.id")) %>% 
  bind_cols(grd_age_wind) %>% 
  #calculate z-scores. get the mean(center) and sd(scale) from the previous z transformation. for consistency
  mutate(journey_number_z = (journey_number - mean(model_data$journey_number, na.rm = T))/(sd(model_data$journey_number, na.rm = T)),
         sqrt_avg_pre_z = (sqrt_avg_pre - mean(model_data$sqrt_avg_pre, na.rm = T))/(sd(model_data$sqrt_avg_pre, na.rm = T)),
         location.lat_z = (location.lat - mean(model_data$location.lat, na.rm = T))/(sd(model_data$location.lat, na.rm = T)))

new_data <- model_data %>% 
  mutate(interaction = "OG_data") %>% 
  dplyr::select(names(new_data_only)) %>%  #only keep the columns that are necessary for the model
  bind_rows(new_data_only) 

# now that we have the values to predict, run the model on them
preds <- predict(migrations_mod, newdata = new_data, type = "link", allow.new.levels = T)

preds_pr <- new_data %>% 
  mutate(preds = preds) %>% 
  rowwise() %>% 
  mutate(probs = gtools::inv.logit(preds)) #https://rpubs.com/crossxwill/logistic-poisson-prob

inter_preds <- preds_pr %>% 
  filter(interaction != "OG_data") 

flap_preds <- inter_preds %>%
  filter(interaction == "age_wind") %>% 
  ggplot(aes(journey_number, sqrt_avg_pre, fill = probs)) +
  geom_tile() +
  labs(x = "Age", y = "Wind speed (m/s)", title = " ") +
  scale_fill_gradientn("Probability", colors = colfunc(200), 
                       n.breaks = 5,
                       values = scales::rescale(c(min(inter_preds$probs),
                                                  mean(inter_preds$probs),
                                                  max(inter_preds$probs)))
  )

png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/acc/13_flaps_fall_windAge.png",
    width = 11.5, height = 8.5, units = "in", res = 500)
ggpubr::ggarrange(flap_coefs +
                    coord_fixed(0.075), 
                  flap_preds +
                    coord_equal(), 
                  nrow = 1)
dev.off()


data3 <- build %>%
  filter(season == "fall") %>% 
  drop_na(behavior) %>% 
  group_by(location.lat) %>% 
  slice(1) %>% 
  ungroup() %>% 
  group_by(journey_number, behavior) %>%
  summarize(obs = n()) %>%
  pivot_wider(names_from = journey_number, values_from = obs) %>%
  dplyr::select(-behavior)

data_percentage <- apply(data3, 2, function(x){x*100/sum(x,na.rm=T)})
# barplot(data_percentage, col = colfunc(2) , border = "white", xlab = "group")

# png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/acc/14_flaps_fall_percent.png",
#     width = 11.5, height = 8.5, units = "in", res = 500)
data_percentage %>% 
  as.data.frame() %>% 
  pivot_longer(cols = c("1", "2", "3"), names_to = "migrations", values_to = "percents") %>% 
  mutate(fob = c(T, T, T, F, F, F)) %>% 
  ggplot(aes(y=percents, x=migrations, fill = migrations)) + 
  geom_bar(position="fill", stat="identity", aes(alpha=fob)) +
  scale_fill_manual(values = colfunc(3)) + 
  guides(fill = "none") +
  scale_alpha_discrete(range = c(0.5, 1)) +
  labs(x = "Migrations", y = "Percentage of falls", alpha = "Followed by \nflapping")
# dev.off()

build %>%
  filter(season == "fall") %>% 
  drop_na(behavior) %>% 
  group_by(location.lat) %>% 
  slice(1) %>% 
  ungroup() %>% 
  ggplot(aes(behavior, exit_height)) +
  geom_violin()+
  geom_boxplot(width = 0.33) +
  facet_wrap(~migration)

build <- build %>% 
  left_join(classified %>% 
              filter(paste0(individual.id, "_", burst_id) %in% build$burst_ind_ID) %>% 
              dplyr::select(trackID, track_status) %>% 
              group_by(trackID) %>% 
              slice(1) %>% 
              ungroup())%>%
  filter(season == "fall") %>% 
  drop_na(behavior) %>% 
  group_by(location.lat) %>% 
  slice(1) %>% 
  ungroup()
table(build$behavior, build$track_status, build$journey_number)
data4 <- table(build$behavior, build$track_status, build$journey_number) %>% 
  as.data.frame() %>% 
  rename(behavior = Var1, status = Var2, migration = Var3) %>% 
  group_by(status, migration) %>% 
  mutate(cent = 100*(Freq/sum(Freq))) %>% 
  ungroup()

# png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/acc/14_flaps_fall_percent.png",
#     width = 11.5, height = 8.5, units = "in", res = 500)
ggplot(data4 %>% mutate(status = ifelse(status == "complete", "Lived", "Died")), aes(migration, cent))+
  geom_bar(aes(fill = migration, alpha = behavior), stat = "identity") +
  geom_hline(yintercept = 90, lty = 2) +
  scale_y_continuous(breaks = c(100, 90, 75, 50, 25, 0)) + 
  guides(fill = "none")  +
  scale_alpha_discrete(range = c(0.5, 1)) +
  scale_fill_manual(values = c("#0081A7", "#FED9B7", "#f27e71")) +
  labs(x = "Fall migration", y = "Percentage", alpha = "") +
  facet_wrap(~status)
# dev.off()


# png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/acc/16_flaps_wind_pre.png",
#     width = 11.5, height = 8.5, units = "in", res = 500)
ggplot(build, aes(behavior, sqrt(avg_wind_speed_pre), fill = behavior)) +
  geom_violin() +
  geom_boxplot(width = 0.33, alpha = 0.5, fill = "white") +
  scale_fill_manual(values = c("#334371", "#6279B8")) +
  labs(x = "", y = "Mean wind speed before falling") +
  facet_wrap(~migration) +
  theme(legend.position = "none")
# dev.off()
# dev.off()
