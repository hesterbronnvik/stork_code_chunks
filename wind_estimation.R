library(parallel)
library(move)
library(moveWindSpeed)
library(geosphere)
# library(data.table)
library(tidyverse)

source("/home/hbronnvik/Documents/chapter2/getWindEstimates_update.R")
source("/home/hbronnvik/Documents/chapter2/thermallingFeaturesFunction.R")
source("/home/hbronnvik/Documents/chapter2/getTrackSegments_updated.R")
source("/home/hbronnvik/Documents/chapter2/getWindEstimate_update.R")
theme_set(theme_classic()+theme(axis.text = element_text(color = "black", size = 12), text = element_text(size = 15)))

wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs") # map projection

# functions:
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
antiwind_direction <- function(u, v){(270-(atan2(v, u)*(180/pi)))%%360 }
# detach("package:plyr", unload = TRUE)

### 1. Identify migrations regardless of success
files <- list.files("/home/hbronnvik/Documents/chapter2/classified_data", full.names = T)

records <- readRDS("/home/hbronnvik/Documents/chapter2/migration_dates_2024-01-15.rds")

f <- files[grepl(2160731880, files)]

# filter for the migration data
migrations_1Hz <- lapply(1:length(files), function(n){
  print(n)
  f <- files[n]
  animal_id <- sub("data/", "", str_split(f, "_")[[1]][2])
  bird <- readRDS(f)
  if(animal_id %in% unique(records$individual.id)){
    iterations <- records %>%
      filter(individual.id == animal_id) %>%
      dplyr::select(trackID) %>% 
      unique() %>% 
      deframe() %>% 
      length()
    record <- records %>% 
      filter(individual.id == animal_id)
    bird <- bird %>% 
      mutate(id_date = paste0(individual.id, "_", date(timestamp))) %>% 
      filter(id_date %in% unique(record$id_date))
    record <- record %>% 
      group_by(trackID) %>% 
      slice(1) %>% 
      ungroup()
    if(nrow(bird) > 0){ # if there are migratory days with 1Hz classified data:
      tracks <- lapply(1:iterations, function(i){
        track <- bird %>% 
          filter(between(timestamp, record$start_migration[i], record$end_migration[i])) %>% 
          mutate(trackID = record$trackID[i],
                 track_status = record$track_status[i],
                 track_displacement = record$track_displacement[i],
                 # the number of days since the onset of this migration
                 m_day = as.numeric(difftime(date(timestamp), date(record$start_migration[i]), units = "days")),
                 # the number of bursts since the onset of the first migration
                 t_burst = burst_id-bird$burst_id[1])
      }) %>% reduce(rbind)
      saveRDS(tracks, file = paste0("/home/hbronnvik/Documents/chapter2/migrations1Hz/", animal_id, "_1Hz.rds"))
      return(tracks)
    }
  }
})

# saveRDS(migrations_1Hz, "/home/hbronnvik/Documents/chapter2/migrations_1Hz_20240112.rds")

# # define the first long-distance displacement for each bird so that we can filter later
# early <- lapply(files, function(file){
#   print(paste0("Processing file ", which(files == file), " of ", length(files), "."), quote = F)
#   data <- readRDS(file)
#   data <- data %>% 
#     as.data.frame() %>% 
#     mutate(date = date(timestamp)) %>% 
#     # reduce the burden for the distance calculation by taking only the first 180 days (migration surely occurs)
#     filter(date < date[1]+days(180)) %>% 
#     group_by(date) %>% 
#     # find the first and last locations of each day and the distance between consecutive locations
#     mutate(start_long = location.long[1],
#            start_lat = location.lat[1],
#            distance = distVincentyEllipsoid(cbind(location.long, location.lat), cbind(lag(location.long), lag(location.lat)))) %>% 
#     ungroup() %>% 
#     # add on the difference between the starting locations for each day
#     mutate(daily_distance = distVincentyEllipsoid(cbind(start_long, start_lat), cbind(lag(start_long), lag(start_lat))))
#   data$daily_distance[1] <- 0 # the first distance always comes out to NA
#   # find the first day of more then 40km displacement
#   first <- data %>% 
#     filter(daily_distance > 40000) %>% 
#     slice(1) %>% 
#     dplyr::select(timestamp) %>%
#     deframe() %>% 
#     date()
#   if(length(first) > 0){
#     record <- data.frame(individual.id = unique(data$individual.id),
#                          individual.local.identifier = unique(data$individual.local.identifier),
#                          start_time = data$timestamp[1],
#                          departure_date = as.numeric(first))
#   }else{
#     record <- data.frame(individual.id = unique(data$individual.id),
#                          individual.local.identifier = unique(data$individual.local.identifier),
#                          start_time = data$timestamp[1],
#                          departure_date = NA)
#   }
#   # saveRDS(data, file = paste0("/home/hbronnvik/Documents/chapter2/clean_data_early/",  unique(data$individual.id),"_21day_migration.rds"))
#   return(record)
# })
# records <- data.table::rbindlist(early) %>% 
#   mutate(departure_date = as.Date(departure_date))
# # saveRDS(records, file = "/home/hbronnvik/Documents/chapter2/start_dates.rds")

# records <- readRDS("/home/hbronnvik/Documents/chapter2/start_dates.rds")
# files <- list.files("/home/hbronnvik/Documents/chapter2/classifiedData/", full.names = T)

# early <- lapply(files, function(file){
#   ind_nodup <- readRDS(file)
#   start <- records %>% 
#     filter(individual.id == ind_nodup@idData$individual.id)
#   ind_nodup <- ind_nodup[timestamps(ind_nodup) < start]
#   # make into a move object
#   mv <- move(x=ind_nodup$location.long, y=ind_nodup$location.lat, 
#              time=ind_nodup$timestamp,
#              proj=crs("+proj=longlat +ellps=WGS84"),
#              animal=ind_nodup$individual.local.identifier,
#              data=ind_nodup)
#   return(mv)
#})


# filter for the migration data
# results <- lapply(files, function(file){
#   print(paste0("Processing file ", which(files == file), " of ", length(files), "."), quote = F)
#   class_df <- readRDS(file)
#   record <- records %>% 
#     filter(individual.id == unique(class_df$individual.id))
#   
#   class_ls <- split(class_df, class_df$trackID)
#   
#   if(nrow(record) > 0){
#     class_df <- class_df %>% 
#       filter(between) 
#     if(nrow(class_df) > 0){
#       class_df <- class_df %>% 
#         mutate(date = date(timestamp)) %>% 
#         group_by(date) %>% 
#         # find the first and last locations of each day and the distance between consecutive locations
#         mutate(start_long = location.long[1],
#                start_lat = location.lat[1]) %>% 
#         ungroup() %>% 
#         # add on the difference between the starting locations for each day
#         mutate(daily_distance = distVincentyEllipsoid(cbind(start_long, start_lat), cbind(lag(start_long), lag(start_lat))),
#                migratory = ifelse(daily_distance > 40000, T, F)) %>% 
#         group_by(date) %>% 
#         mutate(migratory = ifelse(unique(daily_distance)[1] > 40000, T, F)) %>% 
#         ungroup()
#       class_df$daily_distance[1] <- 0 # the first distance always comes out to NA
#     }
#     return(class_df)
#   }
# })
# 155 left (classified and subset)
# early <- early[!sapply(early, function(x) is.null(x))]
# saveRDS(early, file = "/home/hbronnvik/Documents/chapter2/pre_and_21_classified.rds")
# sp <- early[[1]]
# coordinates(sp) <- ~location.long+location.lat
# proj4string(sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
# mapview::mapview(sp, zcol = "migratory")
# detach(package:data.table)
# detach(package:plyr)

# results <- readRDS("/home/hbronnvik/Documents/chapter2/pre_and_21_classified.rds")

results <- readRDS("/home/hbronnvik/Documents/chapter2/migrations_1Hz_20240112.rds")

# the ratios of time spent in each behavior
ratios <- lapply(results, function(r){
  r <- r %>% 
    mutate(date = date(timestamp))
  # the number of classified observations
  n_obs_class <- r %>% 
    group_by(date) %>% 
    summarize(n_obs_class = n())
  r_sum <- r %>% 
    # filter(flightClust_smooth3 != "other") %>% 
    group_by(date, flightClust_smooth3) %>% 
    mutate(obs = n()) %>% 
    ungroup() %>% 
    left_join(n_obs_class, by = join_by(date)) %>% 
    group_by(date, flightClust_smooth3) %>% 
    reframe(flight_ratio = unique(obs)/unique(n_obs_class)) %>% 
    mutate(individual.id = unique(r$individual.id),
           departure_date = records$departure_date[which(records$individual.id == unique(r$individual.id))],
           m_day = as.numeric(difftime(date, departure_date, units = "days")))
  return(r_sum)
}) %>% reduce(rbind)

# find times when there was a glide that lasted at least 30 seconds following a circular soar of 30 seconds
soar_glides <- lapply(results, function(res){
  if(nrow(res)>0){
    record <- records %>% 
      filter(individual.id == unique(res$individual.id))
    # set the first burst of the first day of migration to t_burst == 1
    b1 <- res %>% filter(date(timestamp) == record$departure_date) %>% slice(1) %>% dplyr::select(burstID) %>% deframe()
    # (or just the first classified, migratory burst)
    if(length(b1)==0){
      b1 <- res %>% filter(date(timestamp) > record$departure_date) %>% slice(1) %>% dplyr::select(burstID) %>% deframe()
    }
    if(length(b1) > 0){
      b_no <- data.frame(burstID = unique(res$burstID)) %>% 
        ungroup() %>% 
        mutate(t_burst = which(unique(res$burstID) == burstID)-(which(burstID == b1)-1))
      
      res <- res %>% 
        # set the first day of migration to m_day == 1
        mutate(m_day = as.numeric(difftime(date(timestamp)+days(1), record$departure_date, units = "days"))) %>% 
        group_by(burstID) %>% 
        # produce criteria
        mutate(glide_start = ifelse(flightClust_smooth3 == "gliding" & lag(flightClust_smooth3) == "circular soaring", T, F),
               behav_switch = flightClust_smooth3!=lag(flightClust_smooth3),
               behav_mode = paste0(burstID, "_", cumsum(c(F, behav_switch[2:n()])))) %>% 
        ungroup() %>% 
        group_by(burstID, behav_mode) %>% 
        mutate(obs = n()) %>% 
        ungroup() %>% 
        left_join(b_no, by = join_by(burstID))
      starts <- res %>% 
        filter(lead(glide_start) == T & lead(obs) > 29 & lead(flightClust_smooth3) == "gliding" & flightClust_smooth3 == "circular soaring" & obs > 29) %>% 
        dplyr::select(behav_mode) %>% 
        deframe()
      # filter by those criteria
      res <- res %>% 
        filter(behav_mode %in% starts)
      return(res)
    }
  }
})
soar_glides <- soar_glides[!sapply(soar_glides, function(x) is.null(x))]
soar_glides <- soar_glides[!sapply(soar_glides, function(x) nrow(x)==0)]
# 140 individuals left
soar_heights <- lapply(soar_glides, function(g){
  event <- g %>% 
    group_by(behav_mode) %>% 
    slice(n()) %>% 
    ungroup()
  event
}) %>% reduce(rbind)

p_summary <- soar_heights %>% 
  group_by(m_day) %>% 
  summarize(n = length(unique(behav_mode))) %>% 
  ungroup() %>% 
  mutate(m_day = as.factor(m_day))

png(filename = paste0("/home/hbronnvik/Documents/chapter2/figures/early_thermal_height_day.png"),
    height = 8.3, width = 11.7, units = "in", res = 500)
ggplot(soar_heights %>% filter(m_day %in% c(-6:7)), 
       aes(as.factor(m_day), height.above.ellipsoid, fill = m_day>0)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#0098b0", "#F07167")) +
  labs(x = "Day of the migration 1", y = "Thermal exit height (m)", fill = "Migratory") +
  theme_classic()+
  theme(axis.text = element_text(color = "black", size = 15),
          axis.line = element_line(linewidth = 1.2),
          text = element_text(size = 20),
          legend.text = element_text(size = 15),
          strip.text.x = element_text(size = 15),
          legend.key.size = unit(0.75, 'cm')) +
  annotate(geom="text", x=as.factor(c(-6:7)), y=3200, label=p_summary$n[p_summary$m_day %in% c(-6:7)], color="black", size = 6)
dev.off()
png(filename = paste0("/home/hbronnvik/Documents/chapter2/figures/early_thermal_height_burst.png"),
    height = 8.3, width = 11.7, units = "in", res = 500)
ggplot(soar_heights %>% filter(m_day %in% c(-6:7)), 
       aes(t_burst, height.above.ellipsoid)) +
  geom_point() +
  geom_smooth(method = "lm", aes(color = m_day>0)) +
  scale_color_manual(values = c("#0098b0", "#F07167")) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Recorded flight burst", y = "Thermal exit height (m)", color = "Migratory") +
  theme_classic()+
  theme(axis.text = element_text(color = "black", size = 15),
        axis.line = element_line(linewidth = 1.2),
        text = element_text(size = 20),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(size = 15),
        legend.key.size = unit(0.75, 'cm'))
dev.off()

# info <- info %>%
#   rename(individual.id = animal_id,
#          individual.local.identifier = animal_local_identifier)
# ratios <- ratios %>%
#   left_join(info)
# png(filename = paste0("/home/hbronnvik/Documents/chapter2/figures/flight_time_ratios.png"),
#     height = 8.3, width = 11.7, units = "in", res = 500)
# ggplot(ratios %>% filter(m_day %in% c(-6:7) & flightClust_smooth3 != "other"), 
#        aes(as.factor(m_day), flight_ratio, group = as.factor(m_day), fill = m_day>0)) +
#   geom_boxplot() +
#   scale_fill_manual(values = c("#0098b0", "#F07167")) +
#   labs(x = "Day of the first migration", y = "Ratio of observations to total classified flight", fill = "Migratory") +
#   theme_classic() +
#   facet_wrap(~flightClust_smooth3)+
#   theme(axis.text = element_text(color = "black", size = 15),
#         axis.line = element_line(linewidth = 1.2),
#         text = element_text(size = 20),
#         legend.text = element_text(size = 15),
#         strip.text.x = element_text(size = 15),
#         legend.key.size = unit(0.75, 'cm'))
# dev.off()
# 
# b_data <- ratios %>%
#   filter(m_day %in% c(1:7) & flightClust_smooth3 != "other") %>%
#   group_by(individual.id) %>%
#   mutate(sample_size = n()) %>%
#   ungroup() %>%
#   arrange(sample_size) %>%
#   mutate(individual.id = as.factor(individual.id),
#          individual.id = forcats::fct_reorder(individual.id, sample_size))
# ggplot(b_data, aes(y = as.factor(individual.id), x = flight_ratio, fill = flightClust_smooth3)) +
#   geom_bar(stat = "identity") +
#   theme_classic() +
#   facet_wrap(~m_day)
# 
# colfunc <- colorRampPalette(c("#0081A7", "#0098b0", "#00AFB9", "#7fc4b8", "#FED9B7", "#f7a58f", "#f27e71", "#F07167", "#ee5e53"))
# png(filename = paste0("/home/hbronnvik/Documents/chapter2/figures/ind_var_ratios.png"),
#     height = 8.3, width = 11.7, units = "in", res = 500)
# ggplot(b_data, aes(as.factor(m_day), individual.id, fill = flight_ratio)) +
#   geom_tile(lwd = 0) +
#   scale_fill_gradientn(colors = colfunc(135)) +
#   labs(x = "Day of migration 1", y = "Individual", fill = "Flight time \nspent") +
#   theme_classic() +
#   facet_wrap(~flightClust_smooth3)
# dev.off()
# 
# png(filename = paste0("/home/hbronnvik/Documents/chapter2/figures/mass_ratios.png"),
#     height = 8.3, width = 11.7, units = "in", res = 500)
# ggplot(b_data %>% drop_na(animal_mass) %>% mutate(animal_mass = round(animal_mass, -3)/1000), aes(as.factor(animal_mass), flight_ratio)) +
#   geom_boxplot(aes(fill = as.factor(animal_mass))) +
#   scale_fill_manual(values = c("#0098b0", "#FED9B7", "#F07167")) +
#   labs(x = "Animal mass to the nearest kg", y = "Flight time") +
#   theme_classic() +
#   theme(axis.text = element_text(color = "black", size = 15),
#         axis.line = element_line(linewidth = 1.2),
#         text = element_text(size = 20),
#         legend.text = element_text(size = 15),
#         strip.text.x = element_text(size = 15),
#         legend.position = "none") +
#   facet_wrap(~flightClust_smooth3)
# dev.off()

# png(filename = paste0("/home/hbronnvik/Documents/chapter2/figures/sex_ratios.png"),
#     height = 8.3, width = 11.7, units = "in", res = 500)
# ggplot(b_data %>% drop_na(animal_sex) %>% mutate(animal_sex = ifelse(animal_sex == "m", "Male", "Female")), aes(as.factor(animal_sex), flight_ratio)) +
#   geom_boxplot(aes(fill = as.factor(animal_sex))) +
#   scale_fill_manual(values = c("#0098b0", "#F07167")) +
#   labs(x = "Animal sex", y = "Flight time") +
#   theme_classic() +
#   theme(axis.text = element_text(color = "black", size = 15),
#         axis.line = element_line(linewidth = 1.2),
#         text = element_text(size = 20),
#         legend.text = element_text(size = 15),
#         strip.text.x = element_text(size = 15),
#         legend.position = "none") +
#   facet_wrap(~flightClust_smooth3)
# dev.off()

# results <- readRDS("/home/hbronnvik/Documents/chapter2/migrations_1Hz_20240112.rds")

fls <- list.files("/home/hbronnvik/Documents/chapter2/migrations1Hz", full.names = T)

### 3. Look at metrics of thermal use
thermal_results <- lapply(fls, function(f){
  tr <- readRDS(f)
  tr <- tr %>% 
    filter(flight_clust_sm3 == "circular_soaring") %>% 
    rename(height_above_ellipsoid = height.above.ellipsoid) %>% 
    mutate(ind_burst_id = paste0(individual.id, "_", burst_id))
  tr <- tr %>% 
    group_by(ind_burst_id) %>% 
    mutate(obs = n()) %>% 
    ungroup() %>% 
    filter(obs > 29)
  tr
})
# 154 individuals left
thermal_results <- thermal_results[!sapply(thermal_results, function(x) is.null(x))]

# thermal_results <- data.table::rbindlist(thermal_results)

wind_results <- lapply(101:length(thermal_results), function(n){
  print(n)
  tr <- thermal_results[[n]]
  if(nrow(tr) > 0){
    # start <- records %>% #unique(tr$date[which(tr$migratory==T)])[1]
    #   filter()
    # m_days <- data.frame(date = unique(tr$date), m_day = cumsum(unique(tr$date) >= start))
    # m_days <- m_days %>% 
    #   rownames_to_column(var = "obs") %>% 
    #   mutate(obs = as.numeric(obs),
    #          m_day = ifelse(m_day == 0, obs-nrow(m_days[m_days$m_day == 0,]), m_day)) %>% 
    #   dplyr::select(-obs)
    # t_bursts <- data.frame(burstID = unique(tr$burstID[tr$date >= start]), t_burst = 1:length(unique(tr$burstID[tr$date >= start])))
    # if(length(unique(tr$burstID[tr$date < start]))>0){
    #   t_bursts_pre <- data.frame(burstID = unique(tr$burstID[tr$date < start]), t_burst = 0-(length(unique(tr$burstID[tr$date < start]))-1):0)
    #   t_bursts <- t_bursts_pre %>% rbind(t_bursts)
    # }
    # tr <- tr %>% 
    #   left_join(m_days, by = "date") %>% 
    #   left_join(t_bursts, by = "burstID")
    tr <- tr %>% 
      mutate(time_diff = as.numeric(difftime(timestamp, lag(timestamp), units = "secs")),
             time_diff = ifelse(is.na(time_diff), 1, time_diff),
             min_split = time_diff > 60,
             thermal_event = paste(ind_burst_id, cumsum(min_split)+1, sep = "_")) %>% 
      dplyr::select(-time_diff, -min_split)
    ind <- lapply(split(tr, tr$thermal_event), function(burst){
      burst <- burst %>% as.data.frame()
      mv_burst <- move(x=burst$location.long, y=burst$location.lat, 
                       time=burst$timestamp,
                       proj=wgs,
                       animal=burst$individual.id,
                       data=burst)
      burst_summary <- burst %>%
        mutate(gap = c(1,timeLag(mv_burst)),
               consistent = ifelse(gap == 1, T, F),
               event = cumsum(consistent==0)) %>%
        group_by(event) %>%
        summarize(obs = n()) %>% 
        mutate(sufficient = obs>29)
      if(T %in% unique(burst_summary$sufficient)){
        class_burst <- getWindEstimates(mv_burst)
        class_burst <- as.data.frame(class_burst)
        return(class_burst)
        }
    })
    ind <- ind[!sapply(ind, function(x) is.null(x))]
    ind <- data.table::rbindlist(ind, use.names = T)
    saveRDS(ind, file = paste0("/home/hbronnvik/Documents/chapter2/wind_thermal_data/", unique(ind$individual.id), "_seg_wind_",Sys.Date(),".rds"))
    return(ind)
  }
})
# the wind results using events no more than 60 seconds apart
wind_results <- wind_results[!sapply(wind_results, function(x) is.null(x))]
# saveRDS(wind_results, file = paste0("/home/hbronnvik/Documents/chapter2/wind_thermal_data/short_event_155_seg_wind",Sys.Date(),".rds"))
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
  dplyr::select(trackID, season, journey_number)

thermals <- data.table::rbindlist(classified, use.names = T, fill = T)

# one observation of speed and turning variance per thermal
thermals <- thermals %>% 
  group_by(thermal_event) %>% 
  mutate(wind_speed = wind_speed(windX, windY),
         avg_wind_speed = mean(wind_speed, na.rm = T),
         avg_thermal_str = mean(ThermalStrength, na.rm = T),
         avg_rad = mean(CircRadius, na.rm = T)) %>% 
  slice(1) %>% 
  ungroup() %>% 
  left_join(records, by = join_by(trackID))

library(viridis)
ggplot(thermals, aes(location.lat, avg_wind_speed)) +
  geom_point() +
  scale_x_reverse()
ggplot(thermals %>% filter(location.lat > 46 & m_day < 10), 
       aes(t_burst, log(turn_var_thermal), group = individual.id)) +
  geom_point() +
  geom_smooth(method = "lm", color = "firebrick", se = F)
ggplot(thermals %>% mutate(moday = paste(month(timestamp), day(timestamp), sep = "-"),
                           zone = ifelse(location.lat > 46, "Germany", 
                                         ifelse(location.lat < 17, "Sub-Sahara",
                                                ifelse(location.lat > 17 & location.lat < 36, "Northern Africa", 
                                                       "Southern Europe"))),
                           zone = factor(zone, levels = c("Germany", "Southern Europe", "Northern Africa", "Sub-Sahara"))), 
       aes(t_burst, log(turn_var_thermal), group = zone)) +
  geom_point() +
  geom_smooth(method = "lm", aes(color = zone)) +
  scale_color_viridis(option = "D", discrete = T) +
  scale_x_continuous(n.breaks = 15) +
  labs(x = "Observed burst", y = "log turning variance")
ggplot(thermals %>% filter(season == "fall") %>% mutate(moday = paste(month(timestamp), day(timestamp), sep = "-"),
                                                        zone = ifelse(location.lat > 46, "Germany", 
                                                                      ifelse(location.lat < 17, "Sub-Sahara",
                                                                             ifelse(location.lat > 17 & location.lat < 36, "Northern Africa", 
                                                                                    "Southern Europe"))),
                                                        zone = factor(zone, levels = c("Germany", "Southern Europe", "Northern Africa", "Sub-Sahara"))), 
       aes(zone, vspeed_thermal, fill = zone)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_boxplot() +
  scale_fill_viridis(option = "D", discrete = T) +
  labs(x = "Observed burst", y = "Vertical speed", fill = "Zone") +
  scale_y_continuous(n.breaks = 7) +
  facet_wrap(~journey_number)
ggplot(thermals %>% mutate(moday = paste(month(timestamp), day(timestamp), sep = "-"),
                           zone = ifelse(location.lat > 46, "Germany", 
                                         ifelse(location.lat < 17, "Sub-Sahara",
                                                ifelse(location.lat > 17 & location.lat < 36, "Northern Africa", 
                                                       "Southern Europe"))),
                           zone = factor(zone, levels = c("Germany", "Southern Europe", "Northern Africa", "Sub-Sahara"))), 
       aes(t_burst, avg_thermal_str, group = zone)) +
  geom_point() +
  geom_smooth(method = "lm", aes(color = zone)) +
  scale_color_viridis(option = "D", discrete = T) +
  scale_x_continuous(n.breaks = 15) +
  labs(x = "Observed burst", y = "Thermal strength")

ggplot(thermals %>% filter(season == "fall"), aes(avg_rad, vspeed_thermal)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "gam", color = "firebrick") +
  scale_x_continuous(n.breaks = 13) +
  labs(x = "Average circling radius per thermal (m)", y = "Vertical speed per thermal (m/s)") +
  facet_wrap(~journey_number)
ggplot(thermals %>% filter(season == "fall"), aes(t_burst, avg_rad)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", color = "firebrick") +
  scale_x_continuous(n.breaks = 13) +
  labs(x = "Average circling radius per thermal (m)", y = "Vertical speed per thermal (m/s)") +
  facet_wrap(~journey_number, scales = "free_x")
ggplot(thermals %>% filter(season == "fall"), aes(journey_number, avg_rad, group = journey_number)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_boxplot() +
  labs(y = "Average circling radius per thermal (m)", x = "Fall migration")

my_equation <- paste0("y = ",
                      signif(coef(lm(vspeed_thermal~location.lat, data = thermals))[[2]],3), 
                      "x", " + ",
                      signif(coef(lm(vspeed_thermal~location.lat, data = thermals))[[1]],3))
ggplot(thermals %>% filter(season == "fall"), aes(location.lat, vspeed_thermal)) +
  geom_point() +
  geom_smooth(method = "lm", color = "firebrick") +
  labs(x = "Latitude", y = "Vertical speed per thermal") +
  # scale_x_reverse() + 
  geom_text(x = 20, y = 10, label = my_equation, parse = F)

my_equation <- paste0("y = ",
                      signif(coef(lm(log(turn_var_thermal)~location.lat, data = thermals))[[2]],3), 
                      "x", " + ",
                      signif(coef(lm(log(turn_var_thermal)~location.lat, data = thermals))[[1]],3))
ggplot(thermals %>% filter(season == "fall"), aes(location.lat, log(turn_var_thermal))) +
  geom_point() +
  geom_smooth(method = "lm", color = "firebrick") +
  labs(x = "Latitude", y = "log variance in turning angle per thermal") +
  # scale_x_reverse() + 
  geom_text(x = 20, y = 5, label = my_equation, parse = F)

my_equation <- paste0("y = ",
                     signif(coef(lm(vspeed_thermal~t_burst, data = thermals))[[2]],3), 
                     "x", " + ",
                     signif(coef(lm(vspeed_thermal~t_burst, data = thermals))[[1]],3))

png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/vspeed_thermal.png",
    height = 8.3, width = 11.7, units = "in", res = 500)
ggplot(thermals %>% filter(season == "fall" & t_burst < 9000), aes(t_burst, vspeed_thermal)) +
  geom_point() +
  geom_smooth(method = "lm", color = "firebrick") + 
  # ggpubr::stat_cor()
  geom_text(x = 2000, y = 10, label = my_equation, parse = F) +
  labs(x = "Recorded burst", y = "Vertical speed per thermal")
dev.off()
my_equation <- paste0("y = ",
                      signif(coef(lm(log(turn_var_thermal)~t_burst, data = thermals))[[2]],3), 
                      "x", " + ",
                      signif(coef(lm(log(turn_var_thermal)~t_burst, data = thermals))[[1]],3))

png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/log_turn_var_thermal.png",
    height = 8.3, width = 11.7, units = "in", res = 500)
ggplot(thermals %>% filter(season == "fall"), aes(t_burst, log(turn_var_thermal), group = individual.id)) +
  geom_point() +
  geom_smooth(method = "lm", color = "firebrick") + 
  # ggpubr::stat_cor()
  geom_text(x = 2000, y = 5, label = my_equation, parse = F) +
  labs(x = "Recorded burst", y = "log variance in turning angle per thermal")
dev.off()

png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/ind_log_turn_var_thermal.png",
    height = 8.3, width = 11.7, units = "in", res = 500)
thermals %>% 
  filter(m_day < 8 ) %>% 
  group_by(individual.id) %>% 
  mutate(obs = n()) %>% 
  ungroup() %>% 
  filter(t_burst < 601 & obs > 99) %>% 
  ggplot(aes(t_burst, log(turn_var_thermal), group = individual.id)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "gam", color = "#0081A7", formula = y ~ s(x, bs = "cs", k = 3)) + 
  # ggpubr::stat_cor()
  geom_text(x = 2000, y = 5, label = my_equation, parse = F) +
  labs(x = "Recorded burst", y = "log variance in turning angle per thermal") +
  scale_x_continuous(n.breaks = 13) +
  facet_wrap(~individual.id, scales = "free") +
  theme(axis.text = element_text(color = "black", size = 6), 
        text = element_text(size = 10))
dev.off()

colfunc <- colorRampPalette(c("#0081A7", "#0098b0", "#00AFB9", "#7fc4b8", "#FED9B7", "#f7a58f", "#f27e71", "#F07167", "#ee5e53"))
class_corr <- thermals %>% 
  mutate(week = as.numeric(week(timestamp)),
         turn_var_thermal = log(turn_var_thermal)) %>% 
  drop_na(turn_var_thermal) %>% 
  drop_na(avg_wind_speed) %>% 
  dplyr::select(week, location.lat, t_burst, vspeed_thermal, 
                turn_var_thermal, avg_wind_speed) %>% 
  dplyr::rename(Week = week,
                Latitude = location.lat,
                "Observed soaring event" = t_burst,
                "Vertical speed (m/s)" = vspeed_thermal,
                "Wind speed (m/s)" = avg_wind_speed,
                "Variance in turning angle" = turn_var_thermal) %>% 
  cor()
png(filename = paste0("/home/hbronnvik/Documents/chapter2/figures/look24/correlations.png"),
    height = 8.3, width = 11.7, units = "in", res = 500)
corrplot::corrplot(class_corr, type = "lower", addCoef.col = 'black', col = colfunc(135),tl.col="black")
dev.off()

summary(lm(vspeed_thermal~t_burst+location.lat, data = thermals))
summary(lm(log(turn_var_thermal)~t_burst+location.lat, data = thermals))

thermals %>% 
  group_by(individual.id) %>% 
  slice(1) %>% 
  ungroup() %>% 
  ggplot(aes(yday(timestamp), vspeed_thermal)) +
  geom_point()

thermals %>% 
  filter(trackID == "1176062955_fall_2020" & yday(timestamp) < 241 | trackID == "1176044887_fall_2020" & yday(timestamp) < 241) %>% 
  ggplot(aes(timestamp, vspeed_thermal)) +
  geom_point() +
  geom_smooth(method = "lm", color = "firebrick") +
  facet_wrap(~trackID, scales = "free_x")

# how many observations (thermals) are recorded for each bird
thermals %>% 
  group_by(individual.id) %>% 
  summarize(obs = n()) %>% 
  ggplot(aes(obs)) +
  geom_histogram(bins = 250, color = "black")

thermals %>% 
  group_by(individual.id) %>% 
  summarize(obs = n()) %>% 
  arrange(obs) %>% 
  print(n = 27)

# look at birds with plenty of thermals
thermals_ls <- thermals %>% 
  group_by(individual.id) %>% 
  mutate(obs = n()) %>% 
  ungroup() %>% 
  filter(thermals$journey_number == 1 & thermals$season == "fall" & obs > 44)  %>% 
  group_by(individual.id) %>% 
  group_split()
  
ind_coefs <- lapply(1:length(thermals_ls), function(n){
  i <- thermals_ls[[n]]
  rec <- readRDS("/home/hbronnvik/Documents/chapter2/migration_dates_2024-01-15.rds") %>% 
    filter(trackID == unique(i$trackID))
  
  turn_mod <- lm(log(turn_var_thermal)~t_burst+avg_wind_speed+location.lat, data = i)
  if(!is.na(coef(turn_mod)[2])){
    co_t <- data.frame(estimate = names(coef(summary(turn_mod))[, "Estimate"]),
                       dependent = "log_turn_var_thermal",
                       coef = as.numeric(coef(summary(turn_mod))[, "Estimate"]),
                       std_err = as.numeric(coef(summary(turn_mod))[, "Std. Error"]))
    
    climb_mod <- lm(vspeed_thermal~t_burst+avg_wind_speed+location.lat, data = i)
    co_v <- data.frame(estimate = names(coef(summary(climb_mod))[, "Estimate"]),
                       dependent = "vspeed_thermal",
                       coef = as.numeric(coef(summary(climb_mod))[, "Estimate"]),
                       std_err = as.numeric(coef(summary(climb_mod))[, "Std. Error"]))
    
    ic <- rbind(co_t, co_v)
    ic$individual.id <- unique(i$individual.id)
    ic$sample_size <- length(unique(i$thermal_event))
    ic$track_displacement <- unique(i$track_displacement)
    ic$start_migration <- unique(rec$start_migration)
    ic$track_status <- unique(rec$track_status)
    
    return(ic)
  }
}) %>% 
  reduce(rbind) %>% 
  mutate(significant = sign(coef-std_err) == sign(coef+std_err))

library(viridis)
plot_climb <- ind_coefs %>% 
  filter(estimate == "t_burst" & dependent == "vspeed_thermal") %>% 
  arrange(desc(coef)) %>%
  mutate(id = row_number()) %>% 
  ggplot(aes(id, coef)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = as.numeric(table(ind_coefs$coef[ind_coefs$estimate == "t_burst" & ind_coefs$dependent == "vspeed_thermal"] > 0)[2]),
             lty = 3) +
  geom_pointrange(aes(ymin = coef-std_err, ymax = coef+std_err, color = individual.id), fatten = 0.5) +
  scale_color_viridis(option = "A") +
  labs(y = "Coefficient +/- se", x = "Individual", color = "Individual", title = "Vertical speed by bursts")
plot_turn <- ind_coefs %>% 
  filter(estimate == "t_burst" & dependent == "log_turn_var_thermal") %>% 
  arrange(desc(coef)) %>%
  mutate(id = row_number()) %>% 
  ggplot(aes(id, coef)) +
  geom_hline(yintercept = 0, lty = 2)  +
  geom_vline(xintercept = as.numeric(table(ind_coefs$coef[ind_coefs$estimate == "t_burst" & ind_coefs$dependent == "log_turn_var_thermal"] > 0)[2]),
             lty = 3) +
  geom_pointrange(aes(ymin = coef-std_err, ymax = coef+std_err, color = individual.id), fatten = 0.5) +
  scale_color_viridis(option = "A") +
  labs(y = "Coefficient +/- se", x = "Individual", color = "Individual", title = "Variance in turning by bursts")

# png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/lm_tburst+wind+lat.png",
#     height = 8.3, width = 11.7, units = "in", res = 500)
ggpubr::ggarrange(plot_climb, plot_turn, nrow = 1, common.legend = T, legend = "right")
# dev.off()

plot_climb <- ind_coefs %>% 
  filter(estimate == "t_burst" & dependent == "vspeed_thermal") %>% 
  arrange(sample_size) %>%
  mutate(id = row_number()) %>% 
  ggplot(aes(sample_size, coef)) +
  geom_hline(yintercept = 0, lty = 2) +
  # geom_vline(xintercept = as.numeric(table(ind_coefs$coef[ind_coefs$estimate == "t_burst" & ind_coefs$dependent == "vspeed_thermal"] > 0)[2]),
  #            lty = 3) +
  geom_pointrange(aes(ymin = coef-std_err, ymax = coef+std_err, color = significant), fatten = 0.5) +
  # scale_color_viridis(option = "A") +
  labs(y = "Coefficient +/- se", x = "Sample size", color = "Significant", title = "Vertical speed by bursts")
plot_turn <- ind_coefs %>% 
  filter(estimate == "t_burst" & dependent == "log_turn_var_thermal") %>% 
  arrange(sample_size) %>%
  mutate(id = row_number()) %>% 
  ggplot(aes(sample_size, coef)) +
  geom_hline(yintercept = 0, lty = 2)  +
  # geom_vline(xintercept = as.numeric(table(ind_coefs$coef[ind_coefs$estimate == "t_burst" & ind_coefs$dependent == "log_turn_var_thermal"] > 0)[2]),
  #            lty = 3) +
  geom_pointrange(aes(ymin = coef-std_err, ymax = coef+std_err, color = significant), fatten = 0.5) +
  # scale_color_viridis(option = "A", discrete = T) +
  labs(y = "Coefficient +/- se", x = "Sample size", color = "Significant", title = "Variance in turning by bursts")

# png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/09_coefs_by_n.png",
#     height = 8.3, width = 11.7, units = "in", res = 500)
ggpubr::ggarrange(plot_climb, plot_turn, nrow = 2, common.legend = T, legend = "right")
# dev.off()

ind_coefs %>% 
  filter(estimate == "t_burst" & dependent == "vspeed_thermal") %>% 
  arrange(sample_size) %>%
  mutate(id = row_number()) %>% 
  ggplot(aes(yday(start_migration), coef)) +
  geom_hline(yintercept = 0, lty = 2) +
  # geom_vline(xintercept = as.numeric(table(ind_coefs$coef[ind_coefs$estimate == "t_burst" & ind_coefs$dependent == "vspeed_thermal"] > 0)[2]),
  #            lty = 3) +
  geom_pointrange(aes(ymin = coef-std_err, ymax = coef+std_err, color = sign(coef)), fatten = 0.5) +
  # scale_color_viridis(option = "A") +
  labs(y = "Coefficient +/- se", x = "Sample size", color = "Sign", title = "Vertical speed by bursts")

ind_coefs %>% 
  filter(estimate == "t_burst" & dependent == "vspeed_thermal") %>% 
  arrange(sample_size) %>%
  mutate(id = row_number()) %>% 
  ggplot(aes(log(sample_size), coef, color = sign(coef))) +
  geom_hline(yintercept = 0, lty = 2) +
  # geom_vline(xintercept = as.numeric(table(ind_coefs$coef[ind_coefs$estimate == "t_burst" & ind_coefs$dependent == "vspeed_thermal"] > 0)[2]),
  #            lty = 3) +
  geom_pointrange(aes(ymin = coef-std_err, ymax = coef+std_err), fatten = 0.5) +
  # scale_color_viridis(option = "A") +
  labs(y = "Coefficient +/- se", x = "log sample size", color = "Sign", title = "Vertical speed by bursts") +
  facet_wrap(~track_status)

# png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/08_coefs_by_start.png",
#     height = 8.3, width = 11.7, units = "in", res = 500)
ind_coefs %>% 
  filter(estimate == "t_burst") %>% 
  mutate(dependent = ifelse(dependent == "log_turn_var_thermal", "log turning variance", "Vertical speed (m/s)")) %>% 
  arrange(sample_size) %>%
  mutate(id = row_number()) %>% 
  ggplot(aes(paste0(month(start_migration), "-", day(start_migration)), coef)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_boxplot() +
  labs(x = "Day of the migration onset", y = "lm coefficient") +
  theme(axis.text = element_text(size = 5)) +
  facet_wrap(~dependent)
# dev.off()

aug21 <- ind_coefs %>% 
  filter(estimate == "t_burst" & yday(start_migration) == 233)

check <- thermals %>%
  filter(individual.id %in% unique(aug21$individual.id) & season == "fall" & journey_number == 1) %>% 
  group_by(trackID) %>% 
  sf::st_as_sf(coords = c("location.long", "location.lat"), crs = wgs_sf) #%>%
  # dplyr::summarize(do_union=FALSE) %>%  # do_union=FALSE doesn't work as well
  # sf::st_cast("LINESTRING")
mapview::mapview(check, zcol = "trackID", burst = F)

thermals %>%
  filter(individual.id %in% unique(aug21$individual.id) & season == "fall" & journey_number == 1) %>% 
  group_by(trackID) %>%
  arrange(timestamp) %>% 
  ggplot(aes(location.long, location.lat, color = track_status)) +
  borders(xlim = c(-13, 10), ylim = c(35,50), fill = "gray50") +
  geom_path() +
  facet_wrap(~trackID)

dist_v <- ggplot(ind_coefs %>% filter(estimate == "t_burst" & dependent == "vspeed_thermal"), 
       aes(track_displacement, coef)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_smooth(method = "lm", color = "firebrick") +
  geom_pointrange(aes(ymin = coef-std_err, ymax = coef+std_err), fatten = 0.5) +
  scale_x_continuous(n.breaks = 13) +
  labs(x = "Total displacement (km)", y = "Coefficient +/- se", title = "Vertical speed by bursts")
dist_t <- ggplot(ind_coefs %>% filter(estimate == "t_burst" & dependent == "log_turn_var_thermal"), 
                 aes(track_displacement, coef)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_smooth(method = "lm", color = "firebrick") +
  geom_pointrange(aes(ymin = coef-std_err, ymax = coef+std_err), fatten = 0.5) +
  scale_x_continuous(n.breaks = 13) +
  labs(x = "Total displacement (km)", y = "Coefficient +/- se", title = "Variance in turning by bursts")

# png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/10_coefs_by_dist.png",
#     height = 8.3, width = 11.7, units = "in", res = 500)
ggpubr::ggarrange(dist_v, dist_t)
# dev.off()




classified <- lapply(classified, function(class_df){
  class_df <- class_df %>% 
    filter(dplyr::between(m_day, -6,7))
  class_df
})
classified <- classified[!sapply(classified, function(x) is.null(x))]
classified <- data.table::rbindlist(classified, use.names = T)

# build <- classified %>%
#   mutate(lat = round(location.lat, 3)) %>%
#   ungroup() %>%
#   dplyr::group_by(lat) %>%
#   slice(1) %>%
#   ungroup() %>%
#   dplyr::select(lat, location.long, timestamp) %>%
#   mutate(timestamp = paste0(timestamp, ".000")) %>%
#   rename("location-lat" = lat,
#          "location-long" = location.long)
# write.csv(build, "/home/hbronnvik/Documents/chapter2/sample_locs_wind_mb.csv")
# build <- read.csv("/home/hbronnvik/Documents/chapter2/sample_locs_wind_mb_950-8398241334580478910.csv")
# build <- build %>%
#   mutate(wind_speed_mb = wind_speed(ECMWF.ERA5.PL.U.Wind,ECMWF.ERA5.PL.V.Wind))
# ggplot(build, aes(location.lat, wind_speed)) +
#   geom_point() +
#   geom_smooth(method = "gam") +
#   theme_classic()
# mb <- ggplot(build %>% arrange(wind_speed_mb), aes(location.long, location.lat, color = wind_speed_mb)) +
#   borders(xlim = c(-15,10), ylim = c(0,50), fill = "gray70") +
#   geom_point(alpha = 0.5) +
#   scale_color_gradientn(colors = colfunc(135)) +
#   labs(x = "Longitude", y = "Latitude", color = "ECMWF") +
#   theme_classic()
# mws <- ggplot(build %>% drop_na(wind_speed_mws) %>% arrange(wind_speed_mws), aes(location.long, location.lat, color = wind_speed_mws)) +
#   borders(xlim = c(-15,10), ylim = c(0,50), fill = "gray70") +
#   geom_point(alpha = 0.5) +
#   scale_color_gradientn(colors = colfunc(135)) +
#   labs(x = "Longitude", y = "Latitude", color = "moveWindSpeed") +
#   theme_classic()
# png(filename = paste0("/home/hbronnvik/Documents/chapter2/figures/wind_comparisons.png"),
#     height = 8.3, width = 11.7, units = "in", res = 500)
# ggpubr::ggarrange(mb, mws)
# dev.off()
# sp <- ind
# coordinates(sp) <- ~location.long+location.lat
# proj4string(sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
# m <- mapview::mapview(sp, zcol = "ThermalStrength")
# # saveRDS(m, file = paste0("/home/hbronnvik/Documents/chapter2/figures/", unique(ind$individual.id), "_seg_therm",Sys.Date(),".rds"))
# mapview::mapshot(m, file = "/home/hbronnvik/Documents/chapter2/thermal_strength_map.html")

# just the thermals that lasted over 30 seconds and got wind + thermal estimates
# 149 individuals left
files <- list.files("/home/hbronnvik/Documents/chapter2/wind_thermal_data", full.names = T)
files <- files[grepl("10-10|10-09", files)]

classified <- lapply(files, function(file){
  class_df <- readRDS(file)
  class_df <- class_df %>% 
    filter(dplyr::between(m_day, -6,7))
  class_df
})
classified <- classified[!sapply(classified, function(x) is.null(x))]
classified <- data.table::rbindlist(classified, use.names = T)

classified <- classified %>% 
  mutate(vspeed = vertSpeed_smooth/round(ThermalStrength, digits = 2),
         wind_speed = wind_speed(windX, windY),
         # north is 0 clockwise to 360 (the same as the heading from the tags)
         wind_direction = wind_direction(windX, windY),
         # rotate the wind compass 180 degrees simply to allow ease of comparison to Harel et al. 2016
         antiwind = antiwind_direction(windX, windY),
         cross_wind = cross_wind(windX, windY, heading),
         wind_support = wind_support(windX, windY, heading),
         turn_direction = ifelse(turn.angle > 0, "counter", "clock"),
         directional_change = ifelse(turn_direction == lag(turn_direction), F, T),
         directional_set = cumsum(c(F, directional_change[2:n()]))) %>% 
  group_by(thermal_event) %>% 
  mutate(var_turn = var(turn.angle),
         obs = 1:n(),
         position = as.numeric(scale(obs)),
         n_switches = sum(directional_change)) %>% 
  ungroup()
# %>%
#   left_join(ratios, by = join_by("individual.id", "date", "flightClust_smooth3"))

variables <- data.frame(variable = c("var_turn", "CircRadius", "vspeed", "wind_speed", "vert_speed", "ThermalStrength"),
                        name = c("Variance in turning angle", "Circling radius (m)", "Climb rate", "Wind speed (m/s)", "Vertical speed (m/s)", "Estimated thermal strength (m/s)"))
p_summary <- classified %>% 
  group_by(m_day) %>% 
  summarize(n = length(unique(burst_id))) %>% 
  ungroup() %>% 
  mutate(m_day = as.factor(m_day))

plot_fun <- function(data, x_var, y_var, ylab){
  if((x_var == "m_day")){
    yplace <- ifelse(y_var == "var_turn", 4000, 
                     ifelse(y_var == "CircRadius", 80, 
                            ifelse(y_var == "vspeed", 410, 
                                   ifelse(y_var == "wind_speed", 12, 
                                          ifelse(y_var == "vert_speed", 8, 8.5)))))
    p <- ggplot(data, aes(as.factor(m_day), !!sym(y_var), group = as.factor(m_day), fill = m_day>0)) +
      geom_boxplot() +
      scale_fill_manual(values = c("#0098b0", "#F07167")) +
      labs(x = "Day", y = ylab, fill = "Migrating") +
      theme_classic() + 
      geom_text(data = p_summary, aes(m_day, yplace, label = n)) +
      theme(axis.text = element_text(color = "black"))
  }else{
    p <- ggplot(data %>% drop_na(windX), aes(t_burst, !!sym(y_var))) +
      geom_point(alpha = 0.25) +
      geom_smooth(data = classified %>% filter(m_day<0), color = "#0098b0", method = "lm") +
      geom_smooth(data = classified %>% filter(m_day>0), color = "#F07167", method = "lm") +
      geom_vline(xintercept = 0, lty = 2, color = "gray50") +
      labs(x = "Recorded thermal event", y = ylab) +
      theme_classic() +
      theme(axis.text = element_text(color = "black"))
  }
  return(p)
}
# 
# plot_fun(classified, "m_day", "var_turn", "var")

lapply(1:nrow(variables), function(p){
  png(filename = paste0("/home/hbronnvik/Documents/chapter2/figures/early_", variables$variable[p],"_day.png"),
      height = 8.3, width = 11.7, units = "in", res = 500)
  print(plot_fun(data = classified, x_var = "m_day", y_var = variables$variable[p], ylab = variables$name[p]))
  dev.off()
  png(filename = paste0("/home/hbronnvik/Documents/chapter2/figures/early_", variables$variable[p],"_burst.png"),
      height = 8.3, width = 11.7, units = "in", res = 500)
  print(plot_fun(data = classified, x_var = "t_burst", y_var = variables$variable[p], ylab = variables$name[p]))
  dev.off()
})

png(filename = paste0("/home/hbronnvik/Documents/chapter2/figures/early_strength_shift_burst.png"),
    height = 8.3, width = 11.7, units = "in", res = 500)
ggplot(classified %>% drop_na(ThermalStrength) %>% mutate(label = ifelse(m_day>0, "Migration", "Pre-migration")), aes(position, ThermalStrength)) + 
  geom_point(alpha = 0.05) +
  geom_smooth(method = "gam", color = "firebrick") +
  labs(x = "Relative observation within the soaring event (zero is the middle of the burst)", y = "Estimated thermal strength (m/s)") +
  theme_classic() +
  facet_wrap(~label)
dev.off()


# the ratio of moveWind thermals to MS thermals
nrow(classified[!is.na(classified$ThermalStrength),])/nrow(classified[classified$flightClust_smooth3 == "circular soaring",])

# the proportion of time spent
# wind speed relative to vertical speed
ggplot(classified %>% drop_na(ThermalStrength) %>% mutate(climb_wind = vert_speed/wind_speed), aes(t_burst, climb_wind)) + 
  geom_point(alpha = 0.05) +
  geom_smooth(method = "lm", color = "firebrick") +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Recorded thermal event", y = "Vertical speed given wind speed") +
  theme_classic()
# wind speed relative to thermal speed
ggplot(classified %>% drop_na(ThermalStrength) %>% mutate(thermal_strength_wind = ThermalStrength/wind_speed), aes(t_burst, thermal_strength_wind)) + 
  geom_point(alpha = 0.05) +
  geom_smooth(method = "lm", color = "firebrick") +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Recorded thermal event", y = "Thermal strength given wind speed") +
  theme_classic()
# wind speed relative to circling radii
ggplot(classified %>% drop_na(ThermalStrength) %>% mutate(cirle_rad_wind = CircRadius/wind_speed), aes(t_burst, cirle_rad_wind)) +
  geom_point(alpha = 0.05) +
  geom_smooth(method = "lm", color = "firebrick") +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Recorded thermal event", y = "Circling radius given wind speed") +
  theme_classic()

colfunc <- colorRampPalette(c("#0081A7", "#0098b0", "#00AFB9", "#7fc4b8", "#FED9B7", "#f7a58f", "#f27e71", "#F07167", "#ee5e53"))
# the direction of the wind relative to vertical speed
ggplot(classified %>% drop_na(wind_support) %>% filter(m_day %in% c(-6:7)) %>% mutate(wind_diff = ((antiwind-heading)+180)%%360-180), 
       aes(wind_diff, vert_speed, group = m_day, color = m_day)) +
  geom_smooth() +
  scale_color_gradientn(colors = colfunc(42)) +
  theme_classic()

# (wind_direction+180)%%360-180
sp <- classified %>% filter(thermal_event == unique(thermal_event)[100]) %>% mutate(wind_diff = (((wind_direction-heading)+180)%%360-180))
coordinates(sp) <- ~location.long+location.lat
proj4string(sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
# mapview::mapview(sp, zcol = "flightNum_smooth3")
mapview::mapview(sp, zcol = "wind_diff")

info <- info %>% 
  dplyr::rename(individual.id = animal_id, 
                individual.local.identifier = animal_local_identifier)

classified <- classified %>% 
  left_join(info)
classified %>% 
  mutate(week = week(timestamp)) %>% 
  drop_na(ThermalStrength) %>% 
  dplyr::select(week, location.lat, vert_speed, ThermalStrength, CircRadius, wind_speed, 
                vspeed, var_turn, wind_support, cross_wind, n_switches) %>% 
  corrr::correlate()
class_corr <- classified %>% 
  drop_na(ThermalStrength) %>%
  mutate(week = as.numeric(week(timestamp)),
         n_switches = ifelse(is.na(n_switches), 0, n_switches)) %>% 
  dplyr::select(week, location.lat, t_burst, vert.speed, ThermalStrength, CircRadius, wind_speed, 
                var_turn, wind_support, n_switches) %>% 
  dplyr::rename(Week = week,
         Latitude = location.lat,
         "Observed soaring event" = t_burst,
         "Vertical speed (m/s)" = vert.speed,
         "Thermal strength (m/s)" = ThermalStrength,
         "Circling radius (m)" = CircRadius,
         "Wind speed (m/s)" = wind_speed,
         "Variance in turning angle" = var_turn, 
         "Wind support (m/s)" = wind_support,
         # "Cross wind (m/s)" = cross_wind,
         "Directional switches" = n_switches) %>% 
  cor()
png(filename = paste0("/home/hbronnvik/Documents/chapter2/figures/correlations.png"),
    height = 8.3, width = 11.7, units = "in", res = 500)
corrplot::corrplot(class_corr, type = "lower", addCoef.col = 'black', col = colfunc(135),tl.col="black")
dev.off()

png(filename = paste0("/home/hbronnvik/Documents/chapter2/figures/wind_lat.png"),
    height = 8.3, width = 11.7, units = "in", res = 500)
ggplot(winds %>% drop_na(wind_speed) %>% filter(lat > 35), aes(as.factor(lat), wind_speed, fill = "#0098b0")) +
  geom_violin(trim = F) +
  geom_boxplot(width=0.1, fill="white") +
  # scale_fill_manual(values = colfunc(42)) +
  # geom_smooth(method = "gam", color = "firebrick") +
  labs(x = "Latitude", y = "Wind speed (m/s)") +
  theme_classic() +
  theme(legend.position = "none")
dev.off()



# png(filename = paste0("/home/hbronnvik/Documents/chapter2/figures/sex_speed.png"),
#     height = 8.3, width = 11.7, units = "in", res = 500)
# ggplot(classified %>% 
#          drop_na(vert_speed) %>% 
#          filter(animal_sex %in% c("m", "f")) %>% 
#          mutate(animal_sex = ifelse(animal_sex == "m", "Male", "Female")), 
#        aes(as.factor(animal_sex), vert_speed, fill = animal_sex)) + 
#   geom_boxplot() + 
#   scale_fill_manual(values = c("#0098b0", "#F07167")) +
#   theme_classic() + 
#   labs(x = "Sex", y = "Vertical speed (m/s)") +
#   theme(legend.position = "none")
# dev.off()

build <- classified %>% 
  group_by(thermal_event, directional_set) %>% 
  filter(n()>30) %>% 
  ungroup() %>% 
  group_by(thermal_event) %>% 
  mutate(turn_pos = ifelse(directional_change == F & lead(directional_change) == T, "before",
                           ifelse(directional_change == F & lag(directional_change) == T, "after", 
                                  ifelse(directional_change == T, "during", "unchanged")))) %>% 
  slice(3:(n()-1)) %>% 
  ungroup() %>% 
  drop_na(ThermalStrength)
y <- split(build, build$thermal_event)
strengths <- lapply(1:length(y), function(n)try({
  # print(n)
  x <- y[[n]]
  change <- x %>% 
    filter(directional_change == T) %>% 
    dplyr::select(timestamp) %>% 
    deframe()
  if(length(change) > 0){
    pre <- x %>% 
      filter(timestamp < change)
    l <- ifelse(nrow(pre)>= 15, 15, nrow(pre))
    pre<- pre %>% 
      slice(l) %>% 
      dplyr::select(ThermalStrength) %>% 
      deframe() %>% 
      mean(na.rm = T)
    post <- x %>% 
      filter(timestamp > change) %>% 
      slice(1:15) %>% 
      dplyr::select(ThermalStrength) %>% 
      deframe() %>% 
      mean(na.rm = T)
    df <- data.frame(thermal_event = unique(x$thermal_event),
                     m_day = unique(x$m_day),
                     t_burst = unique(x$t_burst),
                     identity = c("Before", "After"),
                     strength = c(pre, post))
    return(df)
  }
}))
strengths <- strengths[!sapply(strengths, function(x) is.null(x))]
build <- data.table::rbindlist(strengths)

build <- build %>% 
  group_by(thermal_event) %>% 
  mutate(gain_loss = diff(strength)) %>% 
  ungroup()

# library(ggridges)
png(filename = paste0("/home/hbronnvik/Documents/chapter2/figures/strength_a_switch.png"),
    height = 8.3, width = 11.7, units = "in", res = 500)
ggplot(build, aes(x = strength, y = as.factor(m_day), fill = m_day)) + 
  geom_density_ridges(scale = 5) +
  geom_vline(xintercept = 0, lty = 2) +
  scale_fill_gradientn(colors = colfunc(135)) +
  labs(x = "Thermal strength (m/s)", y = "Day of migration 1", fill = "Day") +
  theme_classic() +
  facet_wrap(~identity)
dev.off()

png(filename = paste0("/home/hbronnvik/Documents/chapter2/figures/strength_switch.png"),
    height = 8.3, width = 11.7, units = "in", res = 500)
ggplot(build, aes(x = as.factor(m_day), y = strength, fill = identity)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, lty = 2) +
  scale_fill_manual(values = c("#0098b0", "#F07167")) +
  labs(y = "Thermal strength (m/s)", x = "Day of migration 1", fill = "Position relative \nto directional switch") +
  theme_classic()
dev.off()

ggplot(build, aes(x = as.factor(m_day), y = gain_loss, fill = m_day>0)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, lty = 2) +
  scale_fill_manual(values = c("#0098b0", "#F07167")) +
  labs(y = "Gain or loss of thermal strength (m/s)", x = "Day of migration 1", fill = "Migratory") +
  theme_classic()

png(filename = paste0("/home/hbronnvik/Documents/chapter2/figures/strength_gain_loss.png"),
    height = 8.3, width = 11.7, units = "in", res = 500)
ggplot(build, aes(x = t_burst, y = gain_loss)) +
  geom_point() +
  geom_smooth(method = "lm", aes(color = m_day>0)) +
  scale_fill_manual(values = c("#0098b0", "#F07167")) +
  labs(y = "Gain or loss of thermal strength (m/s)", x = "Recorded thermal event", color = "Migratory") +
  theme_classic()
dev.off()

burst <- build[which(build$thermal_event == build$thermal_event[10000]),]
sp <- burst
coordinates(sp) <- ~location.long+location.lat
proj4string(sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
mapview::mapview(sp, zcol = "turn_pos")

png(filename = paste0("/home/hbronnvik/Documents/chapter2/figures/vspeed_length_stay.png"),
    height = 8.3, width = 11.7, units = "in", res = 500)
ggplot(classified %>% filter(m_day == 1), aes(time_since_tagging, vert.speed)) +
  geom_point() +
  geom_smooth(method = "lm", color = "firebrick") +
  labs(y = "Vertical speed (m/s)", x = "Days since tagging") +
  theme_classic()
dev.off()
ggplot(classified %>% filter(m_day == 1), aes(time_since_tagging, n_switches)) +
  geom_point() +
  geom_smooth(method = "lm", color = "firebrick") +
  labs(y = "Wind support", x = "Days since tagging") +
  theme_classic()

### -----------------------------------------------------------------------------------------
burst <- classified %>% 
  filter(burst_id %in% unique(burst_id)[11:21]) %>% 
  mutate(turn.direction = turn.angle<0, 
         turn_changed = ifelse(turn.direction == lag(turn.direction), F, T))
burst$turn_changed[1] <- F
burst <- burst %>% 
  group_by(thermal_event) %>% 
  mutate(turn_set = cumsum(turn_changed))
# sp <- burst
# coordinates(sp) <- ~location.long+location.lat
# proj4string(sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
# # mapview::mapview(sp, zcol = "flightNum_smooth3")
# mapview::mapview(sp, zcol = "turn.direction")

build <- burst %>% 
  group_by(burst_id,flightNum_smooth3) %>% 
  arrange(timestamp) %>% 
  # slice(1, n()) %>% 
  ungroup() %>% 
  mutate(shift = as.numeric(difftime(timestamp, lag(timestamp), units = "secs")),
         shift = c(0, shift[2:n()]))  %>% 
  group_by(flightNum_smooth3) %>% 
  # slice(1) %>% 
  ungroup() %>% 
  dplyr::select(flightNum_smooth3, timestamp, shift,burstID, burst_id) %>% 
  mutate(soar_changed = ifelse(shift < 180, F, T),
         event = paste(burst_id, as.factor(cumsum(soar_changed)), sep = "_")) 
# %>% 
#   group_by(burst_id,thermal) %>% 
#   mutate(obs = 1:n()) %>% 
#   filter(obs %in% c(1, n())) %>% 
#   mutate(change = as.numeric(difftime(timestamp[n()], timestamp[1], units = "secs"))) %>% 
#   filter(obs == obs[1])
burst <- burst %>%
  left_join(build)

build <- lapply(split(burst, burst$event), function(ev){
  center <- centroid(cbind(ev$location.long, ev$location.lat))
  ev$cent_long <- center[1]
  ev$cent_lat <- center[2]
  ev
}) %>% reduce(rbind)

build <- build %>% 
  mutate(slide = distHaversine(cbind(cent_long, cent_lat), cbind(lag(cent_long), lag(cent_lat)))) %>% 
  group_by(event) %>% 
  mutate(slide = slide[1]/1000,
         slide = ifelse(is.na(slide), 0, slide),
         shift = shift[1]/60,
         new_location = ifelse(round(slide,1) <= 2 & shift <= 180, F, T)) %>% 
  ungroup()

check <- build %>% 
  group_by(event) %>% 
  slice(1) %>% 
  ungroup() %>%
  mutate(new = cumsum(new_location),
         thermal = paste0(burst_id, "_", new)) %>% 
  dplyr::select(cent_long, cent_lat, thermal)


sp <- build %>% left_join(check)
coordinates(sp) <- ~location.long+location.lat
proj4string(sp) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
# mapview::mapview(sp, zcol = "flightNum_smooth3")
mapview::mapview(sp, zcol = "event")
# mapview::mapshot(mapview::mapview(sp, zcol = "event"), url = "/home/hbronnvik/Documents/chapter2/thermal_id_map.html")
