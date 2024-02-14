### Figures required for BLS8
### Hester Bronnvik, M.Sc.
### hbronnvik@ab.mpg.de 
### 2024-02-08


####################################################################################################
### readying the environment:
####################################################################################################
# packages:
library(parallel)
library(move)
library(moveWindSpeed)
library(geosphere)
library(raster)
library(rnaturalearth)
library(rnaturalearthdata)
library(tidyverse)

# functions:
source("/home/hbronnvik/Documents/chapter2/getWindEstimates_update.R")
source("/home/hbronnvik/Documents/chapter2/thermallingFeaturesFunction.R")
source("/home/hbronnvik/Documents/chapter2/getTrackSegments_updated.R")
source("/home/hbronnvik/Documents/chapter2/getWindEstimate_update.R")

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

wgs <- sf::st_crs("+proj=longlat +datum=WGS84 +no_defs")

# metadata:
records <- readRDS("/home/hbronnvik/Documents/chapter2/migration_dates_2024-01-15.rds")
load("/home/hbronnvik/Documents/storkSSFs/loginStored.RData")

file <- "/home/hbronnvik/Documents/storkSSFs/ecmwf/single/boundary_layer_height_land_2022.nc"
blh <- try(raster::raster(file))
template <- blh[[1]]

# plotting:
theme_set(theme_classic()+theme(axis.text = element_text(color = "black", size = 15), text = element_text(size = 17)))
colfunc <- colorRampPalette(c("#0081A7", "#0098b0", "#00AFB9", "#7fc4b8", "#FED9B7", "#f7a58f", "#f27e71", "#F07167", "#ee5e53"))
colfunc_circ <- colorRampPalette(c("#0081A7", "#0098b0", "#00AFB9", "#7fc4b8", "#FED9B7", "#f7a58f", "#f27e71", "#F07167", "#EB3F33",
                                   "#F07167", "#f27e71", "#f7a58f", "#FED9B7", "#7fc4b8", "#00AFB9", "#0098b0", "#0081A7"))
my_labels <- c("Northward", "Westward", "Southward", "Eastward", "Northward")
world <- ne_countries(scale = "medium", returnclass = "sf")

studies <- c(24442409, 212096177, 76367850, 21231406, 1176017658, 173641633)



####################################################################################################
### shaping the data:
####################################################################################################
metad <- lapply(studies, function(x){
  md <- getMovebankReferenceTable(study = x, login = loginStored) %>%
    drop_na(animal_id) %>%
    filter(sensor_type_id == 653) %>% 
    dplyr::select(animal_id, deploy_on_timestamp) %>% 
    mutate(deploy_on_timestamp = date(deploy_on_timestamp)) %>% 
    rename(individual.id = animal_id)
  return(md)
}) %>% reduce(rbind)

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
  mutate(date = as.Date(str_split(id_date, "_")[[1]][2])) %>% 
  ungroup() %>% 
  group_by(trackID) %>% 
  mutate(ld_day = 1:length(unique(date))) %>% 
  ungroup() %>% 
  dplyr::select(trackID, date, ld_day)

thermals <- data.table::rbindlist(classified, use.names = T, fill = T)

# one observation of speed, turning variance, and radius per thermal
thermals <- thermals %>% 
  group_by(thermal_event) %>% 
  mutate(obs = n(),
         date = date(timestamp),
         wind_speed = wind_speed(windX, windY),
         avg_wind_speed = mean(wind_speed, na.rm = T),
         avg_thermal_str = mean(ThermalStrength, na.rm = T),
         avg_rad = mean(CircRadius, na.rm = T),
         avg_windX = mean(windX, na.rm = T),
         avg_windY = mean(windY, na.rm = T),
         avg_wind_direction = wind_direction(avg_windX, avg_windY)) %>% 
  slice(1) %>% 
  ungroup() %>% 
  filter(obs > 29) %>% 
  left_join(records, by = join_by(trackID)) %>% 
  left_join(m_days, by = join_by(trackID, date)) %>% 
  mutate(dst = as.numeric(difftime(timestamp, deploy_on_timestamp, units = "days")),
         migration = as.factor(journey_number))

eastern_birds <- thermals %>% 
  filter(location.long > 16) %>% 
  group_by(individual.id) %>% 
  summarize(eastern = T)

thermals <- thermals %>% 
  filter(!individual.id %in% eastern_birds$individual.id & journey_number != 4) %>% 
  mutate(season = ifelse(season == "fall", "Fall", "Spring"))

# winds
classified <- data.table::rbindlist(classified, use.names = T, fill = T) %>% 
  filter(!individual.id %in% eastern_birds$individual.id) %>% 
  mutate(wind_speed = wind_speed(windX, windY),
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
  ungroup() %>% 
  left_join(records, by = join_by(trackID))

thermals_ls <- thermals %>% 
  group_by(season, journey_number) %>% 
  group_split()

ras_thermals_df <- lapply(thermals_ls, function(r){
  ras_thermals <- r
  sp::coordinates(ras_thermals) <- ~location.long+location.lat
  sp::proj4string(ras_thermals) <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs")
  ras_thermals <- raster::rasterize(x = ras_thermals, y = template, 
                                    field = "avg_wind_direction", fun = mean)
  # raster::plot(ras_thermals)
  ras_thermals_df <- ras_thermals %>% 
    as.data.frame(xy = T) %>% 
    drop_na(layer) %>% 
    mutate(season = unique(r$season),
           migrations = unique(r$journey_number))
  return(ras_thermals_df)
}) %>% reduce(rbind)

# the first and last location of soaring bouts lasting longer than 30 seconds
# annotated using Movebank's Env-DATA service
start_heights <- read.csv("/home/hbronnvik/Documents/chapter2/start_heights-4024091921683662660/start_heights-4277582604051158152.csv") %>% 
  mutate(timestamp = as.POSIXct(timestamp, tz = "UTC"),
         date = date(timestamp),
         altitude = height_above_ellipsoid-ASTER.Elevation) %>% 
  filter(journey_number != 4) %>% 
  left_join(m_days, by = join_by(trackID, date)) %>% 
  mutate(migration = as.factor(journey_number))



####################################################################################################
### formulating plots:
####################################################################################################
thermals %>% 
  filter(season == "Fall" & journey_number == 1) %>% 
  group_by(trackID) %>% 
  arrange(timestamp) %>% 
  slice(1) %>% 
  ungroup() %>% 
  ggplot(aes(dst)) +
  geom_histogram(bins = 15, color = "black") +
  scale_x_continuous(n.breaks = 13) +
  labs(x = "Days since tagging", y = "Count")

hist(thermals$avg_rad, breaks= 100, main = "", col = "gray40")

png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/poster/circling.png",
    height = 8.2, width = 11.7, units = "in", res = 500)
ggplot(thermals, aes(ld_day, avg_rad)) +
  geom_point(alpha = 0.25) +
  geom_path(data = thermals[1,], aes(color = migration), lwd = 2) +
  geom_smooth(method = "lm", aes(color = migration), show.legend = F) +
  scale_color_manual(values = colfunc(3)) +
  scale_x_continuous(breaks = seq(0, 50, by = 10)) +
  labs(x = "Day of migration", y = "Circling radius (m)", color = "Migration") +
  facet_wrap(~season)
dev.off()

hist(thermals$turn_var_thermal, breaks= 100, main = "", col = "gray40")
hist(log(thermals$turn_var_thermal), breaks= 100, main = "", col = "gray40")

png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/poster/turning.png",
    height = 8.2, width = 11.7, units = "in", res = 500)
ggplot(thermals, aes(ld_day, log(turn_var_thermal))) +
  geom_point(alpha = 0.25) +
  geom_path(data = thermals[1,], aes(color = migration), lwd = 2) +
  geom_smooth(method = "lm", aes(color = migration), show.legend = F) +
  scale_color_manual(values = colfunc(3)) +
  scale_x_continuous(breaks = seq(0, 50, by = 10)) +
  labs(x = "Day of migration", y = "log variance in turning angle (deg/s)^2", color = "Migration") +
  facet_wrap(~season)
dev.off()

hist(thermals$vspeed_thermal, breaks= 100, main = "", col = "gray40")
b <- MASS::boxcox(lm(thermals$vspeed_thermal[thermals$vspeed_thermal > 0] ~ 1))
b$x[which.max(b$y)]
hist(sqrt(thermals$vspeed_thermal), breaks= 100, main = "", col = "gray40")

png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/poster/vspeed.png",
    height = 8.2, width = 11.7, units = "in", res = 500)
ggplot(thermals, aes(ld_day, sqrt(vspeed_thermal))) +
  geom_point(alpha = 0.25) +
  geom_path(data = thermals[1,], aes(color = migration), lwd = 2) +
  geom_smooth(method = "lm", aes(color = migration), show.legend = F) +
  scale_color_manual(values = colfunc(3)) +
  scale_x_continuous(breaks = seq(0, 50, by = 10)) +
  labs(x = "Day of migration", y = "sqrt vertical speed (m/s)", color = "Migration") +
  facet_wrap(~season)
dev.off()

hist(thermals$avg_wind_speed, breaks= 100, main = "", col = "gray40")
b <- MASS::boxcox(lm(thermals$avg_wind_speed ~ 1))
b$x[which.max(b$y)]
hist(sqrt(thermals$avg_wind_speed), breaks= 100, main = "", col = "gray40")

png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/poster/wind_speed.png",
    height = 8.2, width = 11.7, units = "in", res = 500)
ggplot(thermals, aes(ld_day, (avg_wind_speed))) +
  geom_point(alpha = 0.25) +
  geom_path(data = thermals[1,], aes(color = migration), lwd = 2) +
  geom_smooth(method = "lm", aes(color = migration), show.legend = F) +
  scale_color_manual(values = colfunc(3)) +
  scale_x_continuous(breaks = seq(0, 50, by = 10)) +
  labs(x = "Day of migration", y = "Mean wind speed per thermal (m/s)", color = "Migration") +
  facet_wrap(~season)
dev.off()

hist(start_heights$altitude, breaks= 100, main = "", col = "gray40")
b <- MASS::boxcox(lm(start_heights$altitude[start_heights$altitude > 0] ~ 1))
lambda <- b$x[which.max(b$y)]
hist((start_heights$altitude^lambda - 1)/lambda, breaks= 100, main = "", col = "gray40")

start_heights <- start_heights %>% 
  mutate(boxcox_alt = (altitude^lambda - 1)/lambda,
         loco_ID = ifelse(loco_ID == "base", "Entering a thermal", "Exiting a thermal"))

png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/poster/heights_lat.png",
    height = 8.2, width = 11.7, units = "in", res = 500)
ggplot(start_heights, aes(location.lat, boxcox_alt)) +
  geom_point(alpha = 0.25) +
  geom_path(data = start_heights[1,], aes(color = migration), lwd = 2) +
  geom_smooth(method = "gam", aes(color = migration), show.legend = F) +
  scale_color_manual(values = colfunc(3)) +
  scale_x_reverse(breaks = c(50, 42, 36, 30, 25, 20, 13),
                  labels = c("50\nFrankfurt", "42\nZaragoza", "36\nTarifa", "30\nMarrakesh",
                             "25\n ", "20\nTimbuktu", "13\nGambia")) +
  labs(x = "Day of migration", y = "Height above ground (m)^0.263", color = "Migration") +
  facet_wrap(~loco_ID) +
  theme(axis.text = element_text(size = 12))
dev.off()

harel_df <- classified %>% 
  filter(journey_number < 4 & season == "fall") %>% 
  drop_na(wind_support) %>% 
  mutate(wind_diff = ((antiwind-heading)+180)%%360-180,
         wind_bin = ifelse(wind_speed < 2, "< 2 m/s",
                           ifelse(wind_speed > 10, "> 10 m/s", 
                                  ifelse(between(wind_speed, 2, 6), "2 to 6 m/s", "6 to 10 m/s"))),
         wind_bin = factor(wind_bin, levels = c("< 2 m/s", "2 to 6 m/s", "6 to 10 m/s", "> 10 m/s")),
         migrations = as.factor(journey_number),
         migrations = factor(migrations, levels = c(1, 2, 3)))

png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/poster/harel.png",
    height = 8.2, width = 11.7, units = "in", res = 500)
ggplot(harel_df, aes(wind_diff, vert_speed, color = migrations)) +
  geom_path(data = harel_df[1,], aes(color = migrations), lwd = 2) +
  geom_smooth(method = "gam", aes(color = migrations), show.legend = F) +
  scale_y_continuous(n.breaks = 10) +
  scale_x_continuous(breaks = c(-180, -90, 0, 90, 180),
                     labels = c("-180\nDirection:", "-90\nLee", "0\nHead", "90\nWind", 180)) +
  scale_color_manual(values = colfunc(3)) +
  labs(x = "Flight vs wind direction (deg)", y = "Vertical speed (m/s)", color = "Migration") +
  facet_wrap(~wind_bin, nrow = 1) +
  theme(axis.text = element_text(size = 10))
dev.off()


