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
         avg_wind_direction = wind_direction(avg_windX, avg_windY),
         avg_bank = mean(BankAngle, na.rm = T)) %>% 
  slice(1) %>% 
  ungroup() %>% 
  filter(obs > 29) %>% 
  left_join(records, by = join_by(trackID)) %>% 
  left_join(m_days, by = join_by(trackID, date)) %>% 
  drop_na(ld_day) %>% 
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

png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/poster/turning_colors.png",
    height = 8.2, width = 11.7, units = "in", res = 500)
ggplot(thermals %>% filter(migration == 1), aes(ld_day, log(turn_var_thermal), color = migration)) +
  geom_point(alpha = 0.5) +
  geom_point(data = thermals %>% filter(migration == 2), alpha = 0.5) +
  geom_point(data = thermals %>% filter(migration == 3), alpha = 0.5) +
  geom_path(data = thermals[1,], aes(color = migration), lwd = 2) +
  geom_smooth(method = "lm", aes(color = migration), show.legend = F) +
  geom_smooth(data = thermals %>% filter(migration == 2), method = "lm", aes(color = migration), show.legend = F)  +
  geom_smooth(data = thermals %>% filter(migration == 3), method = "lm", aes(color = migration), show.legend = F)  +
  scale_color_manual(values = colfunc(3)) +
  scale_x_continuous(breaks = seq(0, 50, by = 10)) +
  labs(x = "Day of migration", y = "log variance in turning angle (deg/s)^2", color = "Migration") +
  facet_wrap(~season)
dev.off()

png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/poster/vspeed_colors.png",
    height = 8.2, width = 11.7, units = "in", res = 500)
ggplot(thermals %>% filter(migration == 1), aes(ld_day, sqrt(vspeed_thermal), color = migration)) +
  geom_point(alpha = 0.5) +
  geom_point(data = thermals %>% filter(migration == 2), alpha = 0.5) +
  geom_point(data = thermals %>% filter(migration == 3), alpha = 0.5) +
  geom_path(data = thermals[1,], aes(color = migration), lwd = 2) +
  geom_smooth(method = "lm", aes(color = migration), show.legend = F) +
  geom_smooth(data = thermals %>% filter(migration == 2), method = "lm", aes(color = migration), show.legend = F)  +
  geom_smooth(data = thermals %>% filter(migration == 3), method = "lm", aes(color = migration), show.legend = F)  +
  scale_color_manual(values = colfunc(3)) +
  scale_x_continuous(breaks = seq(0, 50, by = 10)) +
  labs(x = "Day of migration", y = "sqrt vertical speed (m/s)", color = "Migration") +
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
  labs(x = "Latitude", y = "Height above ground (m)^0.263", color = "Migration") +
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

turn_df <- classified %>% 
  group_by(thermal_event) %>% 
  filter(between(vert.speed, -10, 10) & journey_number < 4) %>% 
  mutate(turn_direction = ifelse(turn.angle > 0, "counter", "clock"),
         directional_change = ifelse(turn_direction == lag(turn_direction), F, T),
         directional_set = cumsum(c(F, directional_change[2:n()])),
         n.changes = length(unique(directional_set))) %>% 
  ungroup() %>% 
  group_by(thermal_event, directional_set) %>% 
  mutate(last_vspeed = vert.speed[n()],
         first_vspeed = vert.speed[1]) %>% 
  slice(1) %>% 
  ungroup() %>% 
  dplyr::select(timestamp, location.long, location.lat, individual.id, trackID, ind_burst_id, thermal_event,
                turn_direction, directional_change, directional_set, n.changes, #before_vspeed, after_vspeed,
                last_vspeed, first_vspeed, season, journey_number)

turn_df$change_vspeed <- NA

for(d in 2:nrow(turn_df)){
  turn_df$change_vspeed[d] <- turn_df$first_vspeed[d]-turn_df$last_vspeed[d-1]
}

hist(turn_df$change_vspeed, breaks= 100, main = "", col = "gray40")
b <- MASS::boxcox(lm(turn_df$change_vspeed[turn_df$change_vspeed > 0] ~ 1))
lambda <- b$x[which.max(b$y)]
hist((turn_df$change_vspeed^lambda - 1)/lambda, breaks= 100, main = "", col = "gray40")

turn_df$boxcox_change <- (turn_df$change_vspeed^lambda - 1)/lambda

png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/poster/direction.png",
    height = 8.2, width = 11.7, units = "in", res = 500)
ggplot(turn_df %>% drop_na(directional_change), 
       aes(x = change_vspeed, y = as.factor(journey_number), fill = as.factor(journey_number))) +
  ggridges::geom_density_ridges2(alpha = 0.5, bandwidth = 0.1) +
  geom_vline(xintercept = 0, lty = 2, color = "gray50") +
  scale_fill_manual(values = colfunc(3)) +
  labs(fill = "Migrations", x = "Change in vertical speed after a turn (m/s)", y = "Density")
dev.off()


turn_df <- classified %>% 
  filter(between(vert.speed, -10, 10) & journey_number < 4) %>% 
  group_by(thermal_event) %>% 
  mutate(turn_direction = ifelse(turn.angle > 0, "counter", "clock"),
         directional_change = ifelse(turn_direction == lag(turn_direction), F, T),
         directional_set = cumsum(c(F, directional_change[2:n()])),
         n.changes = length(unique(directional_set))) %>% 
  ungroup() %>% 
  group_by(thermal_event, directional_set) %>% 
  mutate(obs = n()) %>% 
  filter(obs > 5) %>% 
  slice(1:3, (n()-2):n()) %>%
  mutate(position = c(rep("after", 3), rep("before", 3))) %>% 
  ungroup() %>% 
  group_by(thermal_event, directional_set, position) %>% 
  mutate(avg_pos_speed = mean(vert.speed, na.rm = T)) %>% 
  slice(1) %>% 
  ungroup() %>% 
  group_by(thermal_event) %>% 
  mutate(change_vspeed = avg_pos_speed-lag(avg_pos_speed)) %>% 
  dplyr::select(timestamp, location.long, location.lat, individual.id, trackID, ind_burst_id, thermal_event,
                turn_direction, directional_change, directional_set, n.changes, position,
                avg_pos_speed, n.changes, change_vspeed, season, journey_number)

turn_df <- turn_df %>% 
  filter(thermal_event == lag(thermal_event) & directional_set-lag(directional_set) == 1 & position == "after")

ggplot(turn_df, aes(x = change_vspeed, fill = as.factor(journey_number))) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, lty = 2, color = "gray50") +
  scale_fill_manual(values = colfunc(3)) +
  labs(fill = "Migrations", x = "Change in vertical speed after a turn (m/s)", y = "Density") +
  facet_wrap(~season)

turn_df <- classified %>% 
  group_by(thermal_event) %>% 
  mutate(turn_direction = ifelse(turn.angle > 0, "counter", "clock"),
         directional_change = ifelse(turn_direction == lag(turn_direction), F, T),
         directional_set = cumsum(c(F, directional_change[2:n()])),
         n.changes = length(unique(directional_set)),
         date = date(timestamp)) %>% 
  slice(1) %>% 
  ungroup() %>% 
  left_join(m_days, by = join_by(trackID, date)) %>% 
  dplyr::select(timestamp, location.long, location.lat, individual.id, trackID, ind_burst_id, 
                thermal_event, turn_direction, directional_change, directional_set, n.changes, 
                thermal_duration, season, journey_number, m_day, ld_day)

hist(turn_df$n.changes, breaks= 100, main = "", col = "gray40")
b <- MASS::boxcox(lm(turn_df$n.changes ~ 1))
lambda <- b$x[which.max(b$y)]
hist((turn_df$n.changes^lambda - 1)/lambda, breaks= 100, main = "", col = "gray40")

turn_df <- turn_df %>% 
  mutate(boxcox_n = (n.changes^lambda - 1)/lambda,
         season = ifelse(season == "fall", "Fall", "Spring")) %>% 
  filter(journey_number < 4)


png(filename = "/home/hbronnvik/Documents/chapter2/figures/look24/poster/change_direction.png",
    height = 8.2, width = 11.7, units = "in", res = 700)
ggplot(turn_df, aes(ld_day, boxcox_n)) +
  geom_point(alpha = 0.25) +
  geom_path(data = turn_df[1,], aes(color = as.factor(journey_number)), lwd = 2) +
  geom_smooth(method = "gam", aes(color = as.factor(journey_number)), show.legend = F) +
  scale_color_manual(values = colfunc(3)) +
  labs(x = "Day of migration", y = "Number of directional switches (n)^-0.101", color = "Migration") +
  facet_wrap(~season)
dev.off()

ggplot(turn_df, aes(location.lat, boxcox_n)) +
  geom_point(alpha = 0.25) +
  geom_path(data = turn_df[1,], aes(color = as.factor(journey_number)), lwd = 2) +
  geom_smooth(method = "gam", aes(color = as.factor(journey_number)), show.legend = F) +
  scale_color_manual(values = colfunc(3)) +
  scale_x_reverse() +
  labs(x = "Latitude", y = "Number of directional switches (n)^-0.101", color = "Migration") +
  facet_wrap(~season)




####################################################################################################
### basic statistics:
####################################################################################################

model_data <- thermals %>% 
  mutate(individual.id = as.factor(individual.id),
         migration = as.factor(journey_number),
         sqrt_vspeed = sqrt(vspeed_thermal)) %>% 
  dplyr::select(timestamp, location.long, location.lat, individual.id, ld_day, m_day, t_burst, 
                migration, trackID, ind_burst_id, thermal_event, track_status, season, 
                sqrt_vspeed, turn_var_thermal, avg_wind_speed, avg_rad, avg_bank) %>% 
  mutate_at(c("location.lat","ld_day", "sqrt_vspeed", "turn_var_thermal", 
              "avg_wind_speed", "avg_rad", "avg_bank"), list(z = ~(scale(.)))) %>% 
  mutate(thermal_event = as.factor(thermal_event),
         ind_burst_id = as.factor(ind_burst_id))

library(lme4)

# have skill respond to age et al.
migrations_mod <- lmer(sqrt_vspeed_z ~ migration + avg_wind_speed_z + location.lat_z + (1 | individual.id), 
                       data = model_data, REML = F)
# AIC      BIC   logLik deviance df.resid 
# 124248.0 124309.6 -62117.0 124234.0    49373 
# vertical speed is higher for younger birds (negative effects of age)
# wind speeds increase vertical speeds
 
summary(migrations_mod)
cis <- confint(migrations_mod) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "pred") %>%
  left_join(fixef(migrations_mod) %>% 
              as.data.frame() %>% 
              rename(fixef = ".") %>% 
              rownames_to_column(var = "pred"), 
            by = join_by(pred))
cis %>% 
  mutate(pred = ifelse(pred == "(Intercept)", "migration1", pred)) %>% 
  rename(lower = "2.5 %",
         upper = "97.5 %") %>% 
  filter(pred %in% c("migration1", "migration2", "migration3", "avg_wind_speed_z", "location.lat_z")) %>% 
  mutate(pred = ifelse(pred == "location.lat_z", "Latitude", 
                       ifelse(pred == "avg_wind_speed_z", "Wind speed", pred))) %>% 
  ggplot(aes(fixef, pred)) +
  geom_vline(xintercept = 0, lty = 2, color = "gray40") +
  geom_pointrange(aes(xmin = lower, xmax = upper)) +
  labs(x = "Effect +/- se", y = "Predictor")


long_dist_mod <- lmer(sqrt_vspeed_z ~ ld_day_z + avg_wind_speed_z + location.lat_z + avg_rad + turn_var_thermal +
                         (1 | individual.id) + (1 | migration) + (1 | trackID), 
                       data = model_data, REML = F)
# AIC      BIC   logLik deviance df.resid 
# 124159.4 124221.1 -62072.7 124145.4    49340 
summary(long_dist_mod)
cis2 <- confint(long_dist_mod) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "pred") %>%
  left_join(fixef(long_dist_mod) %>% 
              as.data.frame() %>% 
              rename(fixef = ".") %>% 
              rownames_to_column(var = "pred"), 
            by = join_by(pred))
cis2 %>% 
  mutate(pred = ifelse(pred == "(Intercept)", "intercept", pred)) %>% 
  rename(lower = "2.5 %",
         upper = "97.5 %") %>% 
  filter(pred %in% c("ld_day", "avg_wind_speed_z", "location.lat_z", "avg_rad", "turn_var_thermal")) %>% 
  mutate(pred = ifelse(pred == "location.lat_z", "Latitude", 
                       ifelse(pred == "avg_wind_speed_z", "Wind speed", 
                              ifelse(pred == "ld_day", "Day of \nmigration", pred)))) %>% 
  ggplot(aes(fixef, pred)) +
  geom_vline(xintercept = 0, lty = 2, color = "gray40") +
  geom_pointrange(aes(xmin = lower, xmax = upper)) +
  labs(x = "Effect on vertical speed +/- se", y = "Predictor")

randoms <- ranef(long_dist_mod, postVar = TRUE)
qq <- attr(ranef(long_dist_mod, postVar = TRUE)[[1]], "postVar")

rand.interc<-randoms$individual.id

df<-data.frame(Intercepts=randoms$individual.id[,1],
               sd.interc=2*sqrt(qq[,,1:length(qq)]),
               lev.names=rownames(rand.interc)) %>% 
  arrange(Intercepts) %>% 
  mutate(order = 1:n())

ggplot(df, aes(Intercepts, order)) +
  geom_vline(xintercept = 0, lty = 2, color = "gray40") +
  geom_pointrange(aes(xmin = Intercepts-sd.interc, xmax = Intercepts+sd.interc, color = as.numeric(lev.names))) +
  labs(color = "ID", x = "Effect +/- sd", y = "Individual")


randoms2 <- ranef(long_dist_mod, postVar = TRUE)
qq2 <- attr(ranef(long_dist_mod, postVar = TRUE)[[2]], "postVar")

rand.interc2 <- randoms$migration

df2 <- data.frame(Intercepts=randoms2$migration[,1],
               sd.interc=2*sqrt(qq2[,,1:length(qq2)]),
               lev.names=rownames(rand.interc2)) %>% 
  arrange(Intercepts) %>% 
  mutate(order = 1:n())

ggplot(df2, aes(Intercepts, lev.names)) +
  geom_vline(xintercept = 0, lty = 2, color = "gray40") +
  geom_pointrange(aes(xmin = Intercepts-sd.interc, xmax = Intercepts+sd.interc, color = lev.names)) +
  scale_color_manual(values = colfunc(3)) +
  labs(color = "Migration", x = "Effect +/- sd", y = "Migration")

# install.packages("brms")
# library(brms)
# fit1 <- brm(sqrt_vspeed_z ~ ld_day_z + avg_wind_speed_z + location.lat_z + avg_rad_z + turn_var_thermal_z + 
#               (1|individual.id) + (1|migration),
#             data = model_data, family = gaussian())
# # Warning messages:
# #   1: Rows containing NAs were excluded from the model. 
# # 2: There were 181 divergent transitions after warmup. See
# # https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
# # to find out why this is a problem and how to eliminate them. 
# # 3: Examine the pairs() plot to diagnose sampling problems
# # 
# # 4: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
# # Running the chains for more iterations may help. See
# # https://mc-stan.org/misc/warnings.html#bulk-ess 
# summary(fit1)
# 
# plot(fit1)
# plot(conditional_effects(fit1))

# interaction

model_data <- thermals %>% 
  filter(season == "Fall") %>% 
  mutate(individual.id = as.factor(individual.id),
         migration = as.factor(journey_number),
         sqrt_vspeed = sqrt(vspeed_thermal),
         log_turn_var_thermal = log(turn_var_thermal)) %>% 
  dplyr::select(timestamp, location.long, location.lat, individual.id, ld_day, m_day, t_burst, 
                migration, journey_number, trackID, ind_burst_id, thermal_event, track_status, season, 
                sqrt_vspeed, log_turn_var_thermal, avg_wind_speed, avg_rad, avg_bank) %>% 
  mutate_at(c("location.lat","ld_day", "sqrt_vspeed", "log_turn_var_thermal", 
              "avg_wind_speed", "avg_rad", "avg_bank", "journey_number"), list(z = ~(scale(.)))) %>% 
  mutate(thermal_event = as.factor(thermal_event),
         ind_burst_id = as.factor(ind_burst_id))

migrations_int <- lmer(sqrt_vspeed_z ~ journey_number_z*location.lat_z + (1 | individual.id), 
                       data = model_data, REML = F)
# AIC      BIC   logLik deviance df.resid 
# 105078.5 105129.9 -52533.2 105066.5    38650 
summary(migrations_int)
cis3 <- confint(migrations_int) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "pred") %>%
  left_join(fixef(migrations_int) %>% 
              as.data.frame() %>% 
              rename(fixef = ".") %>% 
              rownames_to_column(var = "pred"), 
            by = join_by(pred))
speed <- cis3 %>% 
  mutate(pred = ifelse(pred == "(Intercept)", "intercept", pred)) %>% 
  rename(lower = "2.5 %",
         upper = "97.5 %") %>% 
  filter(pred %in% c("journey_number_z", "location.lat_z", "journey_number_z:location.lat_z")) %>% 
  mutate(pred = ifelse(pred == "location.lat_z", "Latitude", 
                       ifelse(pred == "journey_number_z", "Migration", "Migration:Latitude"))) %>% 
  ggplot(aes(fixef, pred)) +
  geom_vline(xintercept = 0, lty = 2, color = "gray40") +
  geom_pointrange(aes(xmin = lower, xmax = upper)) +
  labs(x = "Effect on vertical speed +/- se", y = "Predictor")


migrations_int2 <- lmer(turn_var_thermal_z ~ journey_number_z*location.lat_z + (1 | individual.id), 
                       data = model_data, REML = F)
# AIC      BIC   logLik deviance df.resid 
# 111251.4 111303.2 -55619.7 111239.4    41469 
turning <- confint(migrations_int2) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "pred") %>%
  left_join(fixef(migrations_int2) %>% 
              as.data.frame() %>% 
              rename(fixef = ".") %>% 
              rownames_to_column(var = "pred"), 
            by = join_by(pred)) %>% 
  mutate(pred = ifelse(pred == "(Intercept)", "intercept", pred)) %>% 
  rename(lower = "2.5 %",
         upper = "97.5 %") %>% 
  filter(pred %in% c("journey_number_z", "location.lat_z", "journey_number_z:location.lat_z")) %>% 
  mutate(pred = ifelse(pred == "location.lat_z", "Latitude", 
                       ifelse(pred == "journey_number_z", "Migration", "Migration:Latitude"))) %>% 
  ggplot(aes(fixef, pred)) +
  geom_vline(xintercept = 0, lty = 2, color = "gray40") +
  geom_pointrange(aes(xmin = lower, xmax = upper)) +
  labs(x = "Effect on variance in turning angle +/- se", y = "Predictor")

migrations_int3 <- lmer(avg_wind_speed_z ~ journey_number_z*location.lat_z + (1 | individual.id), 
                        data = model_data, REML = F)
# AIC      BIC   logLik deviance df.resid 
# 101308.7 101359.7 -50648.4 101296.7    36267 
winds <- confint(migrations_int3) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "pred") %>%
  left_join(fixef(migrations_int3) %>% 
              as.data.frame() %>% 
              rename(fixef = ".") %>% 
              rownames_to_column(var = "pred"), 
            by = join_by(pred)) %>% 
  mutate(pred = ifelse(pred == "(Intercept)", "intercept", pred)) %>% 
  rename(lower = "2.5 %",
         upper = "97.5 %") %>% 
  filter(pred %in% c("journey_number_z", "location.lat_z", "journey_number_z:location.lat_z")) %>% 
  mutate(pred = ifelse(pred == "location.lat_z", "Latitude", 
                       ifelse(pred == "journey_number_z", "Migration", "Migration:Latitude"))) %>% 
  ggplot(aes(fixef, pred)) +
  geom_vline(xintercept = 0, lty = 2, color = "gray40") +
  geom_pointrange(aes(xmin = lower, xmax = upper)) +
  labs(x = "Effect on flown wind speeds +/- se", y = "Predictor")

ggpubr::ggarrange(speed, turning, winds, nrow = 3)

### limiting

limit_mod <- lme4::lmer(log_turn_var_thermal_z ~ journey_number_z*location.lat_z* + (1 | individual.id), 
           data = model_data[model_data$ld_day < 22,], REML = F)
# AIC      BIC   logLik deviance df.resid 
# 106163.3 106214.7 -53075.7 106151.3    38866 
conf_limit <- confint(limit_mod) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "pred") %>%
  left_join(lme4::fixef(limit_mod) %>% 
              as.data.frame() %>% 
              rename(fixef = ".") %>% 
              rownames_to_column(var = "pred"), 
            by = join_by(pred)) %>% 
  mutate(pred = ifelse(pred == "(Intercept)", "intercept", pred)) %>% 
  rename(lower = "2.5 %",
         upper = "97.5 %") %>% 
  filter(pred %in% c("journey_number_z", "location.lat_z", "journey_number_z:location.lat_z")) %>% 
  mutate(pred = ifelse(pred == "location.lat_z", "Latitude", 
                       ifelse(pred == "journey_number_z", "Migration", "Migration:Latitude")))
ggplot(conf_limit, aes(fixef, pred)) +
  geom_vline(xintercept = 0, lty = 2, color = "gray40") +
  geom_pointrange(aes(xmin = lower, xmax = upper)) +
  labs(x = "Effect on turning variance +/- se", y = "Predictor")

colfunc <- colorRampPalette(c("#0098b0", "#00AFB9", "#7fc4b8", "#FED9B7", "#f7a58f", "#f27e71", "#F07167", "#ee5e53"))
plot_data <- model_data %>% 
  filter(ld_day < 26)
plot3D::scatter3D(plot_data$journey_number_z, plot_data$location.lat, plot_data$log_turn_var_thermal, 
          colvar = plot_data$journey_number_z, col = colfunc(135), alpha = 0.75, pch=16,bty="b2", nticks=5, 
          ticktype="detailed", theta=45, phi=20, xlab="Migration", ylab="Latitude", zlab="Turning")

library(effects)
e <- allEffects(limit_mod)
print(e)
plot(e,multiline=TRUE,confint=TRUE,ci.style="bars")

e1 <- e[[1]]
e.df <- e1 %>% 
  as.data.frame() %>% 
  mutate(Latitude = (location.lat_z*(sd(model_data$location.lat)))+mean(model_data$location.lat))
g <- ggplot(e.df,aes(x=journey_number_z,y=fit,color=Latitude,ymin=lower,ymax=upper)) + 
  geom_pointrange(position=position_dodge(width=.1)) + 
  # viridis::scale_color_viridis(option = "D", limits = c(10, 50)) +
  scale_color_gradientn(colors = rev(viridis::magma(6)[1:5])) +
  labs(x = "Scaled migrations", y = "Scaled turning variance", color = "Latitude")
plot(g)

grd_lat <- expand.grid(journey_number = 1:max(model_data$journey_number),
                       location.lat = seq(from = quantile(model_data$location.lat, 0.025, na.rm = T),
                                          to = quantile(model_data$location.lat, 0.975, na.rm = T),
                                          length.out = 500)) %>% # n = 135
  mutate(log_turn_var_thermal = mean(model_data$log_turn_var_thermal, na.rm = T), #set other variables to their mean
         interaction = "lat_migration")

set.seed(770)
n <- nrow(grd_lat)

new_data_only <- model_data %>%
  slice_sample(n = n, replace = F) %>% #randomly select n of these strata. Make sure n is less than n_distinct(data$stratum)
  # only keep the columns that I need
  dplyr::select(c("individual.id")) %>% 
  bind_cols(grd_lat) %>% 
  #calculate z-scores. get the mean(center) and sd(scale) from the previous z transformation. for consistency
  mutate(log_turn_var_thermal_z = (log_turn_var_thermal - mean(model_data$log_turn_var_thermal))/(sd(model_data$log_turn_var_thermal)),
         journey_number_z = (journey_number - mean(model_data$journey_number))/(sd(model_data$journey_number)),
         location.lat_z = (location.lat - mean(model_data$location.lat, na.rm = T))/(sd(model_data$location.lat, na.rm = T)))

new_data <- model_data %>% 
  mutate(interaction = "OG_data") %>% 
  dplyr::select(names(new_data_only)) %>%  #only keep the columns that are necessary for the model
  bind_rows(new_data_only) 

# now that we have the values to predict, run the model on them
preds <- predict(limit_mod, newdata = new_data, type = "link")

preds_pr <- new_data %>% 
  mutate(preds = preds) %>% 
  rowwise() %>% 
  mutate(probs = gtools::inv.logit(preds)) #https://rpubs.com/crossxwill/logistic-poisson-prob

inter_preds <- preds_pr %>% 
  filter(interaction != "OG_data") 

colfunc <- colorRampPalette(c("#00AFB9", "#7fc4b8", "#FED9B7", "#f7a58f", "#f27e71", "#F07167", "#ee5e53", "#EB3F33"))

tiled <- ggplot(inter_preds, aes(journey_number, location.lat, fill = probs)) +
  geom_tile(lwd = 0) +
  scale_x_continuous(expand = c(0, 0), breaks = c(1, 2, 3)) +
  scale_y_continuous(expand = c(0, 0), n.breaks = 10) +
  labs(x = "Migration", y = "Latitude", title = "") +
  scale_fill_gradientn("Prediction", colors = colfunc(135))

ggpubr::ggarrange(g + ggtitle("Effect on turning variance"), tiled)


limit_mod3 <- lme4::lmer(log_turn_var_thermal_z ~ journey_number_z*location.lat_z*avg_wind_speed_z + 
                          (1 | individual.id), 
                        data = model_data[model_data$ld_day < 22,], REML = F)
# AIC      BIC   logLik deviance df.resid 
# 74886.9  74971.2 -37433.5  74866.9    33877 
conf_limit3 <- confint(limit_mod3) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "pred") %>%
  left_join(lme4::fixef(limit_mod3) %>% 
              as.data.frame() %>% 
              rename(fixef = ".") %>% 
              rownames_to_column(var = "pred"), 
            by = join_by(pred)) %>% 
  mutate(pred = ifelse(pred == "(Intercept)", "intercept", pred)) %>% 
  rename(lower = "2.5 %",
         upper = "97.5 %") %>% 
  slice(4:n())  %>% 
  filter(pred != "journey_number_z") %>%  # this has no real meaning and can't be interpreted
  mutate(pred = gsub("avg_wind_speed", "Wind speed", 
                                       gsub("location.lat", "Latitude", 
                                            gsub("journey_number", "Migrations", 
                                                 gsub("_z", "", pred)))),
         pred = factor(pred, levels = c("Latitude", "Wind speed", "Latitude:Wind speed",
                                        "Migrations:Latitude", "Migrations:Wind speed", 
                                        "Migrations:Latitude:Wind speed")),)

coefs3 <- ggplot(conf_limit3, aes(fixef, forcats::fct_rev(pred))) +
  geom_vline(xintercept = 0, lty = 2, color = "gray40") +
  geom_pointrange(aes(xmin = lower, xmax = upper), color = "#482A65") +
  labs(x = "Effect on turning variance +/- se", y = "Predictor")

grd_3 <- expand.grid(journey_number = (1:max(model_data$journey_number)),
                     avg_wind_speed = seq(from = quantile(model_data$avg_wind_speed, 0.025, na.rm = T), to = quantile(model_data$avg_wind_speed, 0.975, na.rm = T), length.out = 100),
                     location.lat = seq(from = quantile(model_data$location.lat, 0.025, na.rm = T), to = quantile(model_data$location.lat, 0.975, na.rm = T), length.out = 100)) %>% # quantile(prep_d$w_star, .9, na.rm = T)
  mutate(log_turn_var_thermal = mean(model_data$log_turn_var_thermal, na.rm = T),
         interaction = "wind_age_map")

set.seed(770)
n <- nrow(grd_3)

new_data_only <- model_data %>%
  slice_sample(n = n, replace = F) %>% #randomly select n of these strata. Make sure n is less than n_distinct(data$stratum)
  # only keep the columns that I need
  dplyr::select(c("individual.id")) %>% 
  bind_cols(grd_3) %>% 
  #calculate z-scores. get the mean(center) and sd(scale) from the previous z transformation. for consistency
  mutate(log_turn_var_thermal_z = (log_turn_var_thermal - mean(model_data$log_turn_var_thermal))/(sd(model_data$log_turn_var_thermal)),
         journey_number_z = (journey_number - mean(model_data$journey_number))/(sd(model_data$journey_number)),
         location.lat_z = (location.lat - mean(model_data$location.lat, na.rm = T))/(sd(model_data$location.lat, na.rm = T)),
         avg_wind_speed_z = (avg_wind_speed - mean(model_data$avg_wind_speed, na.rm = T))/(sd(model_data$avg_wind_speed, na.rm = T)))

new_data <- model_data %>% 
  mutate(interaction = "OG_data") %>% 
  dplyr::select(names(new_data_only)) %>%  #only keep the columns that are necessary for the model
  bind_rows(new_data_only) 

# now that we have the values to predict, run the model on them
preds <- predict(limit_mod3, newdata = new_data, type = "link")

preds_pr <- new_data %>% 
  mutate(preds = preds) %>% 
  rowwise() %>% 
  mutate(probs = gtools::inv.logit(preds)) #https://rpubs.com/crossxwill/logistic-poisson-prob

inter_preds <- preds_pr %>% 
  filter(interaction == "wind_age_map") 

colfunc <- colorRampPalette(c("#0081A7", "#0098b0", "#00AFB9", "#7fc4b8", "#FED9B7", "#f7a58f", "#f27e71", "#F07167", "#ee5e53"))

tiled <- inter_preds %>% 
  mutate(journey_number = ifelse(journey_number == 1, "Migration 1", 
                                  ifelse(journey_number == 2, "Migration 2", "Migration 3"))) %>% 
  ggplot(aes(location.lat, avg_wind_speed, fill = probs)) +
  geom_raster(interpolate = F) +
  scale_y_continuous(expand = c(0, 0), breaks = 1:10) +
  scale_x_reverse(expand = c(0, 0), breaks = c(45, 35, 25)) +
  labs(y = "Wind speed (m/s)", x = "Latitude", title = "") +
  scale_fill_gradientn("Trend", colors = colfunc(10)) +
  facet_wrap(~journey_number)

ggpubr::ggarrange(coefs3 + ggtitle("Effect on turning variance"), tiled)










ggplot(thermals, aes(location.lat, log(turn_var_thermal))) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm", aes(color = as.factor(journey_number))) +
  scale_x_reverse() +
  scale_color_manual(values = colfunc(3), name = "Migration")

ggplot(thermals %>% filter(season == "Fall"), aes(as.factor(ld_day), sqrt(vspeed_thermal), fill = as.factor(journey_number))) +
  geom_boxplot() +
  # geom_smooth(method = "gam", aes(ld_day, sqrt(vspeed_thermal), lty = as.factor(journey_number)),
  #             color = "black", show.legend = T) +
  scale_fill_manual(values = colfunc(3)) +
  # scale_color_manual(values = colfunc(3)) +
  labs(x = "Day of migration", y = "sqrt vertical speed", fill = "Migrations", lty = "Migrations")

ggplot(thermals %>% filter(season == "Fall"), aes(as.factor(ld_day), log(turn_var_thermal), fill = as.factor(journey_number))) +
  geom_boxplot() +
  geom_smooth(method = "gam", aes(ld_day, log(turn_var_thermal), lty = as.factor(journey_number)),
              color = "black", show.legend = T) +
  scale_fill_manual(values = colfunc(3)) +
  # scale_color_manual(values = colfunc(3)) +
  labs(x = "Day of migration", y = "log turning variance", fill = "Migrations", lty = "Migrations")

# library(brms)
bayes_migrations <- brm(sqrt_vspeed_z ~ journey_number_z*location.lat_z + (1|individual.id),
                        data = model_data[!is.na(model_data$sqrt_vspeed_z),], family = gaussian())
plot(bayes_migrations)
summary(bayes_migrations)
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)     0.19      0.02     0.16     0.22 1.01      634     1234

# library(posterior)

migration_sum <-
  as_draws_df(bayes_migrations) %>% 
  summarise_draws()

migration_sum %>% 
  head(n = 10)

fac_labs <- c("Bulk ESS", "Tail ESS")
names(fac_labs) <- c("ess_bulk", "ess_tail")

migration_sum %>% 
  pivot_longer(starts_with("ess")) %>% 
  ggplot(aes(x = value)) +
  geom_vline(xintercept = 4000, lty = 2) +
  geom_histogram(bins = 50, color = "black") +
  xlim(0, NA) +
  theme(panel.grid = element_blank()) +
  facet_wrap(~ name, labeller = labeller(name = fac_labs))

migration_sum %>% 
  ggplot(aes(x = rhat)) +
  geom_vline(xintercept = 1, color = "white") +
  geom_histogram(bins = 50, color = "black") +
  theme(panel.grid = element_blank())

# saveRDS(model_data, file = "/home/hbronnvik/Documents/chapter2/model_data.rds")
model_data <- readRDS("/home/hbronnvik/Documents/chapter2/model_data.rds")
# AR1_migrations <- brms::brm(sqrt_vspeed_z ~ journey_number_z*location.lat_z + (1|individual.id) + ar(p = 1, cov = TRUE),
#                         data = model_data, family = gaussian())
AR1_migrations <- brms::brm(sqrt_vspeed_z ~ journey_number_z*location.lat_z + (1|individual.id),
                            data = model_data[!is.na(model_data$sqrt_vspeed_z),], family = gaussian())
