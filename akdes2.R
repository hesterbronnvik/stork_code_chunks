### CTMM AKDE building
### Hester Bronnvik 
### hbronnvik@ab.mpg.de
### 2023-01-23

library(move)
library(ctmm)
library(lubridate)
library(raster)
library(terra)
library(tidyverse)

wind_support <- function(u,v,heading) {
  angle <- atan2(u,v) - heading/180*pi
  return(cos(angle) * sqrt(u*u+v*v))
  14
}
cross_wind <- function(u,v,heading) {
  angle <- atan2(u,v) - heading/180*pi
  return(sin(angle) * sqrt(u*u+v*v))
}
# CubeRoot<-function(x){
#   sign(x)*abs(x)^(1/3)
# }
# 
# deardorff <- function(blh, T2m, s_flux){
#   g <- 9.80665
#   wT <- (-1*s_flux)/1013/1.2 # upward is negative, reverse the sign
#   w_star <- CubeRoot(g*blh*wT/T2m)
#   return(w_star)
# }
# 
# w_star <- function(g = 9.81, blh, T2m, s_flux, m_flux) {
#   
#   z <- blh
#   T_k_2m <- T2m
#   T_c_2m <- T_k_2m - 273.15
#   Thetav_k_z <- (T_k_2m) + 0.006 * z
#   wT <- (s_flux * -1) / 1013 / 1.2 #reverse the sign. ECMWF upward fluxes are negative
#   wq <- (m_flux * -1) *1000 /1.2 #reverse the sign. ECMWF upward fluxes are negative
#   
#   wthetav <- wT + 0.61 * T_c_2m * wq
#   
#   w_star <- CubeRoot(g*z*(wthetav/Thetav_k_z))
#   
#   return(w_star)
#   
# }

## Get the data
# read in the saved data
locs <- readRDS("C:/Users/hbronnvik/Documents/storkSSFs/soc_info_migration_locations_2023_05_02.rds")

# the total number of migrations attempted and completed by each animal in each season 
meta <- locs %>%
  mutate(season = ifelse(grepl("fall", trackID), "post", "pre")) %>% 
  group_by(individual.id, season) %>% 
  count(trackID) %>% 
  mutate(journey_number = row_number()) %>% 
  ungroup() %>% 
  group_by(individual.id) %>%
  mutate(total_journeys = n()) %>% 
  ungroup() #%>% 
#   rowwise() %>% 
#   mutate(track_status = unique(locs$track_status[locs$trackID == trackID])) %>% 
#   ungroup()
# table(meta$track_status)
# add the number of journeys
locs <- locs %>% 
  rowwise() %>% 
  mutate(journey_number = meta$journey_number[which(meta$trackID == trackID)]) %>% 
  ungroup()

# align the locations in time
locs$datestamp <- locs$timestamp
year(locs$datestamp) <- 2024
second(locs$datestamp) <- 00
locs$datestamp <- round_date(locs$datestamp, "hour")

# get the outlines of land to remove the possibility of information over water
# tmax <- raster::getData('worldclim', var = "tmax", res = 10)
# mask <- crop(raster(tmax, 1), raster::extent(min(first_locs$location.long), max(first_locs$location.long), min(first_locs$location.lat), max(first_locs$location.lat)))
# raster::plot(mask)
# # the Straits of Gibraltar are at most 60km across
# mm <- buffer(mask, 30000)
# raster::plot(mm)
# writeRaster(mm, "C:/Users/hbronnvik/Documents/storkSSFs/buffered_water_ras.grd")
# mm <- raster::raster("C:/Users/hbronnvik/Documents/storkSSFs/buffered_water_ras.grd")
# outlines <- raster::rasterToPolygons(mm, dissolve=TRUE)
# 
# # make a list of the unique hours
# dates <- split(locs, locs$datestamp)
# # keep only list elements with at least 10 data points so that the model can fit
# dates <- dates[lapply(dates, nrow) > 5]
# # choose an extent over which the kdes are calculated
# ee1 <- extent(-1721608, 2004455, 0, 6876759)
# # build AKDEs
# start_time <- Sys.time()
# full_akdes <- lapply(dates, function(x){
#   tryCatch({
#     print(unique(x$datestamp))
#     # manipulate duplicated timestamps because these data are not from a single bird
#     # thus, there are by design duplicated times that make fitting the model impossible
#     # here, we add 5 minutes to the second of a duplicate so that two birds that transmitted simultaneously are offset
#     if(sum(duplicated(round_date(x$timestamp, "minute"))) > 0){
#       x$timestamp[which(x$timestamp %in% x$timestamp[duplicated(round_date(x$timestamp, "minute"))])] <- x$timestamp[which(x$timestamp %in% x$timestamp[duplicated(round_date(x$timestamp, "minute"))])] + minutes(5)
#     }
#     # create a telemetry object to use the ctmm functions
#     ind <- x %>%
#       # erase identities because ctmm automatically detects them
#       mutate(individual.id == 1,
#              individual.local.identifier = 1) %>% 
#       # align the locations in space
#       as.telemetry(projection = "ESRI:54009")
#     # fit a model to the data
#     guess <- ctmm.guess(ind, interactive=FALSE)
#     fit <- ctmm.select(ind, guess)
#     # calculate the KDE using the land outlines and the model
#     UD <- akde(ind, fit, SP = outlines, grid = list(dr = 30000, extent = ee1))
#     saveRDS(UD, file = paste0("C:/Users/hbronnvik/Documents/storkSSFs/KDES/40x70x30day/grid_", gsub("-|:| ", "", unique(x$datestamp)), ".rds"))
#     return(UD)
#     # gc()
#   }, error = function(e){print(geterrmessage())})
# })
# Sys.time() - start_time

### annotate the locations from the burst tracks (step_generation.R) with akde values
# the tracks:
# burst_data <- list.files("C:/Users/hbronnvik/Documents/storkSSFs/burst_data/speed40x80", full.names = T)
# burst_data <- lapply(burst_data, read.csv) %>% 
#   reduce(rbind) 
# 
# burst_data <- burst_data %>% 
#   mutate(timestamp = as.POSIXct(timestamp, tz = "UTC"),
#          alignment = round_date(timestamp, unit = "hour"))
# year(burst_data$alignment) <- 2024

# # the KDEs: 
# kde_files <- list.files("C:/Users/hbronnvik/Documents/storkSSFs/KDES/gridded/list", pattern = "grid", full.names = T)
# kde_f <- str_sub(kde_files, 49, 52)
# kde_files <- kde_files[order(str_sub(kde_files, 48, 52))]
# kdes <- readRDS(kde_files[1])
# for (i in 2:length(kde_files)) {
#   kde <- readRDS(kde_files[i])
#   kdes <- c(kdes, kde)
#   rm(kde)
# }
# 

# 
# library(gifski)
# gif_file <- tempfile(fileext = ".gif")
# gifski(imgs, gif_file)
# unlink(imgs)
# utils::browseURL(gif_file)
# save_gif(gifski(imgs, gif_file), gif_file, 1000, 1500, res = 1000)
# 
# kdes_names <- lapply(kdes, function(x){
#   n <- x@info$identity
# }) %>% unlist()
# 
# # find missing timestamps:
# ts <- burst_data %>% 
#   group_by(alignment) %>% 
#   slice(1) %>% 
#   ungroup() %>% 
#   mutate(alignment = as.character(alignment)) %>% 
#   filter(!alignment %in% kdes_names) %>% 
#   dplyr::select(alignment) %>% 
#   as.vector()

# retrieve the weather data annotations 
files_sl <- list.files("C:/Users/hbronnvik/Documents/storkSSFs/movebank/used/", pattern = "SL-", full.names = T)
con <- file(files_sl[1], "r")
df <- read.csv(con, nrows=10)
close(con)
coln <- colnames(df)

a_data_SL <- lapply(files_sl, function(x){
  sl <- read.csv(x, row.names = NULL)
  sl <- sl %>% 
    dplyr::select(coln)
  return(sl)
}) %>% 
  reduce(rbind)

colnames(a_data_SL)[15:18] <- c("blh", "s_flux", "tempK", "m_flux")
a_data_SL$X <- NULL

files_pl <- list.files("C:/Users/hbronnvik/Documents/storkSSFs/movebank/used/", pattern = "PL-", full.names = T)
con <- file(files_pl[10], "r")
df <- read.csv(con, nrows=10)
close(con)
coln <- colnames(df)

a_data_PL <- lapply(files_pl, function(x){
  df <- read.csv(x, row.names = NULL)
  df <- df %>% 
    dplyr::select(coln)
  # print(x)
  # print(table(is.na(df[, 17])))
  return(df)
}) %>% 
  reduce(rbind)

colnames(a_data_PL)[15:17] <- c("u_wind", "v_wind", "pvv")
a_data_PL$X <- NULL

# combine the data at single and at pressure levels, calculate the wind support on each step, and add the AKDEs
a_data <- a_data_SL %>%
  full_join(a_data_PL) %>% 
  mutate(s_wind = wind_support(u_wind, v_wind, heading),
         c_wind = cross_wind(u_wind, v_wind, heading),
         timestamp = as.POSIXct(timestamp, tz = "UTC"),
         alignment = round_date(timestamp, unit = "hour"))
# align the burst data locations in time
year(a_data$alignment) <- 2024
second(a_data$alignment) <- 00
a_data$alignment <- round_date(a_data$alignment, "hour")

# make the home ranges for all but a focal animal 
HR <- lapply(unique(a_data$individual.id), function(x){
  ind <- locs %>% 
    mutate(individual.id = as.character(datestamp),
           individual.local.identifier = as.character(datestamp)) %>% 
    group_by(individual.id) %>% 
    mutate(count = n()) %>% 
    ungroup() %>% 
    filter(count >= 5)%>% 
    filter(datestamp %in% unique(a_data$alignment))
  year(ind$timestamp) <- 2024
  
  # take all of the migratory locations within the hours of migrations by birds that survived
  ind <- ind %>% 
    filter(datestamp %in% unique(a_data$alignment))
  # make a list of telemetry objects
  ind_tele <- as.telemetry(ind, projection = "ESRI:54009")
  # saveRDS(ind_tele, file = "C:/Users/hbronnvik/Documents/storkSSFs/ind_tele.rds")
  # make a list of model fits
  start_time <- Sys.time()
  it <- lapply(ind_tele, function(x){
    tryCatch({
      guess <- ctmm.guess(x, interactive=FALSE)
      # fit the models
      fit <- ctmm.select(x, guess)
      print(paste0("Fitted a model for ", x$timestamp[1], "."), quote = F)
      pair <- list(x, fit)
      # saveRDS(pair, file = paste0("C:/Users/hbronnvik/Documents/storkSSFs/ctmm_fits/", gsub("-|:| ", "",x$timestamp[1]), ".rds"))
      return(pair)
    }, error = function(e){print(geterrmessage(), quote = F)})
  })
  Sys.time() - start_time #Time difference of 35.47424 mins
  # split these apart for speed
  fits <- lapply(it, function(fit){
    f <- fit[[2]]
  })
  tels <- lapply(it, function(data){
    d <- data[[1]]
  })
  # call in a buffered raster of land outlines + 30km
  mm <- raster::raster("C:/Users/hbronnvik/Documents/storkSSFs/buffered_water_ras.grd")
  outlines <- raster::rasterToPolygons(mm, dissolve=TRUE)
  # estimate utilization distributions of all the storks in each hour cropped to the buffered map
  start_time <- Sys.time()
  kdes <- akde(tels, fits, SP = outlines, grid = list(dr = 30000, extent = extent(-2468787, 3006683, 0, 6876759)))
  Sys.time() - start_time
})
# saveRDS(kdes, file = "C:/Users/hbronnvik/Documents/storkSSFs/40x70x30day_kdes_ls_20230502.rds")
kdes <- readRDS("C:/Users/hbronnvik/Documents/storkSSFs/40x70x30day_kdes_ls_20230502.rds")
## Plot the AKDEs to see that they are as expected
library(maptools)
library(broom)
library(magick)
data("wrld_simpl")
ws <- crop(wrld_simpl, extent(-20,20,0,60))
spdf_fortified <- tidy(ws, region = "ISO3")

# kdes_sub <- kdes[1:1000]
lapply(kdes, function(x){
  # make plotting data:
  ras <- raster::raster(x, DF = "PDF")
  ras <- terra::rast(ras)
  ras <- terra::project(ras, "+proj=longlat +datum=WGS84 +no_defs")
  ras <- raster::raster(ras)
  ras <- crop(ras, extent(-20,20,0,60))
  rp <- rasterToPoints(ras)
  df <- as.data.frame(rp)
  colnames(df)[3] <- "Probability"
  # make the plot:
  png(paste0("C:/Users/hbronnvik/Documents/storkSSFs/may_ras_plot/", gsub("-| |:", "", x@info$identity), ".png"),
      width = 8, height = 11, units = "in", res = 100)
  print(ggplot() +
          geom_raster(data = df, aes(x = x, y = y, fill = Probability)) +
          scale_fill_viridis(option = "B") +
          geom_polygon(data = spdf_fortified, aes(x = long, y = lat, group = group), fill= NA, color="white")  +
          coord_cartesian(xlim = c(-20, 20), ylim = c(0, 60)) +
          labs(title = sub("2024-", "     Time:  ", x@info$identity)) +
          # geom_polygon(data = spdf_fortified, aes(x = long, y = lat, group = group), fill= NA, color="white")  +
          theme_void() +
          theme(panel.background = element_rect(fill = "black", colour = "black"),
                plot.background = element_rect(fill = "black", colour = "black"),
                plot.title = element_text(colour = "white"),
                legend.text=element_text(color="white"),
                legend.title = element_text(color = "white", size = 10),
                plot.margin = unit(c(1,1,1,1), "cm"),
                legend.position = "none"))
  dev.off()
})
df <- as.data.frame(rp)
colnames(df)[3] <- "Probability"
df$Probability[which(df$Probability < 1.5e-32)] <- NA
ggplot() +
  geom_raster(data = df, aes(x = x, y = y, fill = Probability)) +
  scale_fill_continuous_sequential(palette = "terrain", na.value = "white") +
  geom_polygon(data = spdf_fortified, aes(x = long, y = lat, group = group), fill= NA, color="black")  +
  coord_cartesian(xlim = c(-20, 20), ylim = c(0, 60)) +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", colour = "white"))

imgs <- list.files("C:/Users/hbronnvik/Documents/storkSSFs/may_ras_plot/", full.names = TRUE)
imgs <- c(imgs[grep("20240505", imgs)[1]:length(imgs)], imgs[1:grep("20240505", imgs)[1]-1])
av::av_encode_video(imgs, 'C:/Users/hbronnvik/Documents/storkSSFs/may_ras_plot.mp4', framerate = 100)
utils::browseURL('C:/Users/hbronnvik/Documents/storkSSFs/may_ras_plot.mp4')
# go through the KDEs to annotate the tracks
akde_data <- lapply(kdes, function(x){
  # the hour covered by this KDE
  kde_hour <- x@info$identity
  print(kde_hour)
  # the locations in that hour
  tracks <- a_data %>% 
    filter(alignment == kde_hour)
  if(nrow(tracks) > 0){
    # convert to a spatial object
    coordinates(tracks) <- ~ location.long + location.lat
    proj4string(tracks) <- CRS("EPSG:4326")
    tracks <- spTransform(tracks, "ESRI:54009")
    # the KDE as a raster ("PDF" gives the average probability density per cell)
    ud <- raster(x, DF = "PDF")
    # extract the UD value
    vals <- raster::extract(ud, tracks)
    # append
    df <- a_data %>% 
      filter(alignment == kde_hour) %>% 
      mutate(UD_PDF = vals)
    return(df)
  }
}) %>% 
  discard(is.null)
library(data.table)
akde_data <- rbindlist(akde_data)
# akde_data$X <- NULL

a_data <- a_data %>% 
  full_join(akde_data)

# calculate w*, the proxy for thermal uplift velocity
# a_data$w_star <- deardorff(a_data$blh, a_data$tempK, a_data$s_flux)
a_data$w_star <- w_star(blh = a_data$blh, T2m = a_data$tempK, s_flux = a_data$s_flux, m_flux = a_data$m_flux)

ggplot(a_data, aes(x = w_star, group = ind, fill = ind)) +
  geom_density(adjust = 1.5, alpha = 0.4) +
  scale_fill_manual(values = c("#7F9AE5", "#B8DE29FF")) +
  labs(x = "w*", y = "Density", fill = "Season") +
  theme_classic()

# add on the ages of the birds and the ratio of used to average available per stratum
rs_ids <- meta %>% 
  dplyr::filter(trackID %in% a_data$track) %>% 
  rename(track = trackID)

a_data <- a_data %>% 
  full_join(rs_ids[, c("journey_number", "track")])%>% 
  mutate(timestamp = as.POSIXct(timestamp, tz = "UTC"),
         alignment = round_date(timestamp, unit = "hour"))
year(a_data$alignment) <- 2024 

# the absolute values of w*
# ggplot(a_data %>% 
#          mutate(ind = ifelse(ind == "post", "Post-breeding", "Pre-breeding")), aes(as.character(journey_number), w_star, fill = as.factor(used))) +
#   geom_boxplot() +
#   scale_fill_manual(values = c("#C2474E", "#EE8434")) +
#   labs(x = "Migration", y = "w*", fill = "Use case") +
#   theme_classic() +
#   facet_wrap(~ind)
# # per stratum, is used w* higher or lower than the average available value?
# pd <- lapply(split(a_data, a_data$stratum), function(x){
#   u <- x[x$used == 1,]
#   a <- x[x$used == 0,]
#   if(u$w_star/mean(a$w_star) > 10){
#     print(u$w_star/mean(a$w_star))
#     print(unique(x$stratum))
#   }
#   x$rat <- u$w_star/median(a$w_star)
#   return(x)
# })
# library(data.table)
# pd <- rbindlist(pd)
# pd <- pd %>% 
#   group_by(stratum) %>% 
#   slice(1) %>% 
#   ungroup() %>% 
#   mutate(ind = ifelse(ind == "post", "Post-breeding", "Pre-breeding"))
# ggplot(pd , aes(as.character(journey_number), rat, fill = ind)) +
#   geom_hline(yintercept = 1, lty = 2) +
#   geom_boxplot() +
#   scale_fill_manual(values = c("#B8DE29FF", "#7F9AE5")) +
#   ylim(c(0, 3)) +
#   labs(x = "Migration", y = "Used w*/ median available", fill = "Use case") +
#   theme_classic()

# data missing the akdes
check <- a_data %>% 
  filter(is.na(UD_PDF)) %>% 
  group_by(alignment) %>% 
  slice(1)

# todo <- locs %>% 
#   filter(datestamp %in% check$alignment)
# 
# fit_names <- lapply(fits, function(x){
#   n <- x@info$identity
# }) %>% unlist()
# 
# fits <- fits[which(fit_names %in% as.character(todo$datestamp))]
# tels <- tels[which(fit_names %in% as.character(todo$datestamp))]
# 
# leftover_kdes <- lapply(which(fit_names %in% as.character(todo$datestamp))[1:30], function(x){
#   print(x)
#   kde <- akde(tels[[x]], fits[[x]], SP = outlines,  grid = list(dr = 30000, extent = ee1))
#   saveRDS(kde, file = paste0("C:/Users/hbronnvik/Documents/storkSSFs/KDES/gridded/grid_", x, ".rds"))
#   return(kde)
# })

# set the missing UD values to zero
a_data$UD_PDF[is.na(a_data$UD_PDF)] <- 0

mod_data <- a_data %>% 
  mutate_at(c("s_wind", "c_wind", "blh", "pvv", "UD_PDF", "step_length", "turning_angle"),
            list(z = ~(scale(.))))

i_mod <- function(df) {
  require(survival)
  clogit(used ~ UD_PDF_z + blh_z + s_wind_z + c_wind_z + step_length_z + turning_angle_z + strata(step_id), data = df)
}

# 2023 tracks came back from annotation empty, remove them
mod_data <- mod_data %>% 
  filter(!grepl("2023", track))

res <- lapply(split(mod_data, mod_data$track), function(x){
  print(unique(x$track))
  model <- as.data.frame(i_mod(x)$coefficients) %>% 
    rownames_to_column()
  colnames(model) <- c("variable", "coefficient")
  model$season <- ifelse(grepl("fall", unique(x$track)), "Fall", "Spring")
  model$journey_number <- as.character(unique(x$journey_number))
  model$track <- unique(x$track)
  return(model)
}) %>% reduce(rbind)

samples <- res %>%
  group_by(season, journey_number) %>% 
  summarize(n = length(unique(track)))

pd <- res %>% 
  filter(variable != "step_length_z" & variable != "turning_angle_z") %>% 
  mutate(variable = ifelse(variable == "UD_PDF_z", "Social information", 
                           ifelse(variable == "blh_z", "Uplift", "Supporting wind")))

png("C:/Users/hbronnvik/Documents/storkSSFs/clogit_UD.BLH.WS.png", width = 33.87, height = 19.05, units = "cm", res = 500)
print(ggplot(pd, aes(journey_number, coefficient, fill = variable)) +
        geom_hline(yintercept = 0, lty = 2) +
        geom_boxplot() +
        scale_fill_manual(values = c("#6b4596ff",  "#de7065ff", "#f7cb44ff", "#ef789999")) +
        # scale_fill_viridis(option = "B", discrete = T) +
        ylim(-3, 3) +
        labs(y = "Coefficient", x = "Migration") +
        theme_classic() +
        theme(legend.position = "none",
              text = element_text(size = 20),
              axis.text = element_text(color = "black")) +
        geom_text(data = samples, inherit.aes = F, aes(journey_number, Inf, label = n), vjust = 1.5) +
        facet_wrap(~season+variable))
dev.off()

# check the daily variance
v <- mod_data %>% 
  mutate(date = yday(alignment)) %>% 
  group_by(stratum) %>% 
  mutate(varU = var(na.omit(blh))) %>% 
  ungroup() %>% 
  group_by(date) %>% 
  mutate(n = length(unique(stratum))) %>% 
  ungroup() %>% 
  mutate(sem = varU/sqrt(n))
png("C:/Users/hbronnvik/Documents/storkSSFs/stratum_sdsqrtn_blhz_per_date.png", width = 33.87, height = 19.05, units = "cm", res = 500)
print(ggplot(v, aes(journey_number, sem, group = journey_number)) +
        geom_violin() +
        # geom_bar(stat = "identity") +
        labs(x = "Migration", y = "Stratum variance") +
        theme_classic() +
        theme(text = element_text(size = 20, color = "black"),
              axis.text = element_text(size = 15, color = "black")))
dev.off()

#Welch’s Test for Unequal Variances
t.test(v$varU[v$journey_number == 1], v$varU[v$journey_number == 5])
#Bartlett's test 
bartlett.test(varU ~ journey_number, data = v)
#perform Welch's ANOVA
model <- oneway.test(blh ~ factor(journey_number), data = mod_data, var.equal = FALSE)
#perform a normal ANOVA
model <- aov(varU ~ factor(journey_number), data = v)
#load DescTools package
library(DescTools)
#perform Scheffe's test
ScheffeTest(model)
# compute Kruskal-Wallis test
kruskal.test(varU ~ factor(journey_number), data = v)
pairwise.wilcox.test(v$varU, factor(v$journey_number),
                     p.adjust.method = "BH")

library(ggridges)
v %>% 
  filter(varU < 1e+05) %>%
  ggplot(aes(varU, factor(journey_number), fill = journey_number)) +
  geom_density_ridges2(alpha = 0.8,quantile_lines=TRUE,
                       quantile_fun=function(x,...)mean(x)) +
  labs(x = "Per stratum variance in BLH", y = "Migration", fill = "Migration", title = "")+
  scale_fill_viridis(option = "B")+
  coord_cartesian(expand = F)+
  theme_minimal() +
  theme(legend.position = "none", text = element_text(size=15),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 20),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"))

# add variables for ID and year, and scale the predictors
a_data <- a_data %>% 
  filter(!grepl("2023", track)) %>% 
  mutate(season = ifelse(grepl("fall", track), "fall", "spring")) %>% 
  group_by(season) %>%
  mutate(migration = journey_number,
         id1 = factor(individual.id),
         id2 = factor(individual.id),
         id3 = factor(individual.id),
         id4 = factor(individual.id)#,
         # yr1 = factor(year(timestamp)),
         # yr2 = factor(year(timestamp))
         ) %>% 
  ungroup() %>% 
  mutate_at(c("s_wind", "blh", "pvv", "UD_PDF", "migration"),
            list(z = ~(scale(.))))

meta_sub <- meta %>% 
  filter(total_journeys > 3)

subset <- a_data %>% 
  filter(!grepl("2023", track))  %>% 
  filter(track %in% meta$trackID) %>% 
  mutate(season = ifelse(grepl("fall", track), "fall", "spring")) %>% 
  filter(season != "spring")

# export the annotated and prepared data for modeling on the MPCDF Raven HPC cluster
saveRDS(subset, file = paste0("C:/Users/hbronnvik/Documents/storkSSFs/annotated_data_subset_repeat_", Sys.Date(), ".rds"))

# the following is run by the MPCDF: ---------------------------------------------------------------

library(INLA)
a_data <- readRDS("C:/Users/hbronnvik/Documents/storkSSFs/annotated_data_2023-05-03.rds")
a_data <- a_data %>% 
  mutate_at(c("step_length", "turning_angle", "s_wind", "blh", "UD_PDF"),
            list(z = ~(scale(.))))
a_data1 <- a_data[which(a_data$season == "fall"),]
a_data1 <- a_data1 %>% 
  filter(!grepl("2023", track)) %>% 
  mutate(migration = journey_number,
         id1 = factor(individual.id),
         id2 = factor(individual.id),
         id3 = factor(individual.id)) %>% 
  ungroup()

formula_w <- used ~ -1 + s_wind_z + w_star_z + UD_PDF_z + migration_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) +
  f(id1, w_star_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(id2, s_wind_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(id3, UD_PDF_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(id4, UD_PDF_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))


mean.beta <- 0
prec.beta <- 1e-4 

formula_w <- used ~ -1 + UD_PDF_z*blh_z +
  f(stratum, model = "iid",
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) +
  f(id1, UD_PDF_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))),
    group = migration, control.group = list(model = "ar1", cyclic = F, scale.model = TRUE)) +
  f(id2, blh_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))),
    group = migration, control.group = list(model = "ar1", cyclic = F, scale.model = TRUE))

#model without missing values
M_post <- inla(formula_w, family = "Poisson",
               control.fixed = list(
                 mean = mean.beta,
                 prec = list(default = prec.beta)),
               data = a_data1,
               num.threads = 10,
               control.predictor = list(compute = TRUE), #this means that NA values will be predicted.
               control.compute = list(openmp.strategy = "huge", config = TRUE, cpo = T))

saveRDS(M_post, file = "M_post.rds")

#extract model validation results
eval <- data.frame(CPO = mean(M_post$cpo$cpo, na.rm = T),
                   Mlik = as.numeric(M_post$mlik[,1][2]))

graph <- as.data.frame(summary(M_post)$fixed)
colnames(graph)[which(colnames(graph)%in%c("0.025quant","0.975quant"))]<-c("Lower","Upper")
colnames(graph)[which(colnames(graph)%in%c("0.05quant","0.95quant"))]<-c("Lower","Upper")
colnames(graph)[which(colnames(graph)%in%c("mean"))]<-c("Estimate")
graph$Factor <- rownames(graph)

ggplot(graph, aes(Estimate, Factor)) +
  geom_pointrange(aes(xmin = Lower, xmax = Upper), color = "#6B0F68") +
  geom_point(color = "#6B0F68") +
  geom_vline(xintercept = 0, lty = 2) +
  theme_classic()

summ_rndm  <- M_post$summary.random
tab <-  M_post$summary.random$id1
tab$year <- c(tab$ID[1], unique(M_post$.args$data$journey_number))

names(tab)[c(4, 6)] <- c("q025", "q975")

ggplot(tab[2:4,], aes(y = year, x = mean)) + 
  geom_pointrange(aes(xmin = q025, xmax = q975), color = "#6B0F68") +
  geom_point(color = "#6B0F68") +
  geom_vline(xintercept = 0, lty = 2) +
  theme_classic()

M_pre <- inla(formula_w, family = "Poisson",
              control.fixed = list(
                mean = mean.beta,
                prec = list(default = prec.beta)),
              data = a_data2,
              num.threads = 10,
              control.predictor = list(compute = TRUE), #this means that NA values will be predicted.
              control.compute = list(openmp.strategy = "huge", config = TRUE, cpo = T))
saveRDS(M_pre, file = "M_pre.rds")

# Set up the model, but do not yet fit it
TMBStruc <- glmmTMB(used ~ s_wind_z + w_star_z + UD_PDF_z + migration_z +
                     (1|stratum) + (0 + s_wind_z | id1) + (0 + w_star_z | id1) + (0 + UD_PDF_z | id1) + (0 + migration_z | id1),
                   family=poisson,
                   data=a_data1,
                   doFit=FALSE)
# TMBStruc <- glmmTMB(used ~ blh_z*UD_PDF_z*migration_z +
#                       (1|stratum) + (0 + blh_z | id1) + (0 + UD_PDF_z | id1) + (0 + migration_z | id1),
#                     family=poisson,
#                     data=subset,
#                     doFit=FALSE)
# # Fix the standard deviation of the first random term, which is the (1|stratum) component
# in the above model equation
TMBStruc$parameters$theta[1] <- log(1e3)
TMBStruc$mapArg <- list(theta=factor(c(NA,1:4)))
# Fit the model
m_post <- glmmTMB::fitTMB(TMBStruc)
summary(m_post)

library(TwoStepCLogit)
log_post <- Ts.estim(used ~ blh_z*UD_PDF_z + strata(stratum) + cluster(as.factor(journey_number)),
         random = ~ blh_z*UD_PDF_z,
         data = subset,
         D="UN(1)")
f <- data.frame(beta = log_post$beta, se = log_post$se)
f$variable <- rownames(f) 
f <- f %>% 
  mutate(variable = ifelse(grepl("\\:", variable), "Storks X Uplift",
                           ifelse(grepl("blh", variable), "Uplift", "Stork density")))
ggplot(f, aes(beta, variable, color = variable)) +
  geom_errorbar(aes(xmin = beta-se, xmax = beta+se), width = 0, lwd = 5) +
  geom_point(size = 15) +
  scale_color_manual(values = c("#A63A50", "#f7cb44ff", "#f9a242ff")) +
  labs(x = "β", y = "", title = "Two step conditional logistic regression", color = "Variable") +
  scale_x_continuous(breaks = c(-1, 0, 1, 2, 3, 4)) +
  geom_vline(xintercept = 0, lty = 2, lwd = 2) +
  theme_classic() +
  theme(text = element_text(color = "black", size = 60),
        axis.text = element_text(color = "black", size = 50),
        axis.line = element_line(colour = 'black', linewidth = 3))

r <- as.data.frame(log_post$r.effect) %>% 
  rownames_to_column(var = "migration") %>% 
  rename(Uplift = blh_z,
         "Stork density"  = UD_PDF_z,
         "Storks X Uplift" = "blh_z:UD_PDF_z") %>% 
  pivot_longer(cols = Uplift:"Storks X Uplift", names_to = "variable", values_to = "beta")
ggplot(r, aes(migration, beta, color = variable, group = variable)) +
  # geom_pointrange(aes(xmin = beta-se, xmax = beta+se), color = "#A63A50") +
  geom_hline(yintercept = 0, lty = 2, lwd = 2) +
  geom_point(size = 15) +
  geom_line(lty = 2, lwd = 5) +
  scale_color_manual(values = c("#A63A50", "#f7cb44ff", "#f9a242ff")) +
  labs(y = "β", x = "Migration", title = "Two step conditional logistic regression", color = "Variable") +
  # scale_x_continuous(breaks = c(-1, 0, 1, 2, 3, 4)) +
  geom_vline(xintercept = 0, lty = 2, lwd = 2) +
  theme_classic() +
  theme(text = element_text(color = "black", size = 60),
        axis.text = element_text(color = "black", size = 50),
        axis.line = element_line(colour = 'black', linewidth = 3))

w_mod <- function(df) {
  require(survival)
  clogit(used ~ s_wind_z + w_star_z + UD_PDF_z + step_length_z + turning_angle_z + strata(step_id), data = df)
}
res <- lapply(split(a_data, a_data$track), function(x){
  model <- as.data.frame(w_mod(x)$coefficients) %>% 
    rownames_to_column()
  colnames(model) <- c("variable", "coefficient")
  model$season <- ifelse(unique(x$ind) == "post", "Post-breeding", "Pre-breeding")
  model$journey_number <- as.character(unique(x$journey_number))
  model$track <- unique(x$track)
  return(model)
}) %>% reduce(rbind)
pd <- res %>% 
  filter(variable != "step_length_z" & variable != "turning_angle_z") %>% 
  mutate(variable = ifelse(variable == "s_wind_z", "Wind support",
                           ifelse(variable == "w_star_z", "Thermal uplift", 
                                  "Social density")))
ggplot(pd, aes(journey_number, coefficient, fill = variable)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_boxplot() +
  scale_fill_manual(values = c("#6b4596ff", "#de7065ff", "#f7cb44ff")) +
  ylim(-3, 3) +
  labs(y = "Model coefficient", x = "Migration") +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size = 20),
        axis.text = element_text(color = "black")) +
  facet_wrap(~season+variable)

# ---------------------------------------------------------------------------------------------------

start_time <- Sys.time()
kdes <- lapply(fitz, function(fit){
  # fit <- readRDS(fitz[x])
  UD <- akde(fit[[1]], fit[[2]], SP = outlines, grid = list(dr = 30000))
})
Sys.time() - start_time
# plot(kdes[[10]])
check <- lapply(fitz, readRDS)
checker <- lapply(1:length(check), function(x){
  df <- check[[x]][[2]]
  df <- data.frame(len = length(df), no = x)
  return(df)
}) %>% reduce(rbind)

# create a list of hours so that an AKDE can be built for each hour
# create a list of sampling periods so that an AKDE can be built for each one
# first_ls <- first_locs %>% 
#   # temporarily reduce the data to test the function:
#   slice(1:1000) %>% 
#   mutate(day_time = paste(month(timestamp), day(timestamp), sampling_period, sep = "_"))
first_ls <- split(locs[which(locs$journey_number == 1),], locs$datestamp[which(locs$journey_number == 1)])
first_ls <- first_ls[lapply(first_ls, nrow) > 10]
# test 
late_ls <- late_ls[which(names(late_ls) %in% names(first_ls))]

start_time <- Sys.time()
first_akdes <- lapply(first_ls, function(x){
  # create a telemetry object to use the ctmm functions
  ind <- x %>%
    # erase identities because ctmm automatically detects them
    mutate(individual.id == 1,
           individual.local.identifier = 1) %>% 
    # align the locations in space
    as.telemetry(projection = "ESRI:54009")
  # fit a model to the data
  fit <- ctmm.fit(ind)
  # calculate the KDE using the land outlines and the model
  UD <- akde(ind, fit, SP = outlines)
  return(UD)
})

# split the data according to stage
# late_locs <- locs %>% 
#   filter(stage == "adult") %>% 
#   # align the locations in hour
#   mutate(datestamp = round_date(datestamp, unit = "hour"),
#          hour = hour(timestamp),
#          sampling_period = ifelse(hour >= 0 & hour < 6, 1, 
#                                   ifelse(hour >= 6 & hour < 12, 2, 
#                                          ifelse(hour >= 12 & hour < 18, 3, 4))))


# create a list of hours so that an AKDE can be built for each hour
late_ls <- split(late_locs, late_locs$datestamp)
# create a list of sampling periods so that an AKDE can be built for each one
# late_ls <- late_locs %>% 
#   # temporarily reduce the data to test the function:
#   slice(1:1000) %>% 
#   mutate(day_time = paste(month(timestamp), day(timestamp), sampling_period, sep = "_"))
late_ls <- split(locs[which(locs$journey_number != 1),], locs$datestamp[which(locs$journey_number != 1)])
late_ls <- late_ls[lapply(late_ls, nrow) > 10]
# test
late_ls <- late_ls[which(names(late_ls) %in% names(first_ls))]

late_akdes <- lapply(late_ls, function(x){
  # create a telemetry object to use the ctmm functions
  ind <- x %>%
    # erase identities because ctmm automatically detects them
    mutate(individual.id == 1,
           individual.local.identifier = 1) %>% 
    # align the locations in space
    as.telemetry(projection = "ESRI:54009")
  # fit a model to the data
  fit <- ctmm.fit(ind)
  # calculate the KDE using the land outlines and the model
  UD <- akde(ind, fit, SP = outlines)
  return(UD)
})
Sys.time() - start_time

first_akdes <- readRDS("C:/Users/hbronnvik/Documents/storkSSFs/first_akdes_02.02.rds")
late_akdes <- readRDS("C:/Users/hbronnvik/Documents/storkSSFs/late_akdes_02.03.rds")

# combine the adult and juvenile akdes in time
akdes <- lapply(1:length(names(late_akdes)), function(x){
  # the juvenile akde that matches an adult akde
  id <- which(names(first_akdes) == names(late_akdes)[x])
  if(length(id) == 1){# conversion to a raster of the PDF
    f <- raster(first_akdes[[id]], DF = "PDF")
    # using the terra package for reprojection
    f <- terra::rast(f)
    f <- terra::project(f, "EPSG:4326")
    l <- raster(late_akdes[[x]], DF = "PDF")
    l <- terra::rast(l)
    l <- terra::project(l, "EPSG:4326")
    l <- terra::project(l, f)
    # combining the two akdes
    fl <- terra::merge(f, l)
    # adding the date as the name
    names(fl) <- names(first_akdes)[[id]]
    print(paste0("Combined ", names(late_akdes)[x], "."), quotes = F)
    return(fl)}
})
# the hours that have information from both adults and juveniles:
akdes <- akdes[lapply(akdes, is.null) == F]

only_first <- which(!names(first_akdes) %in% names(late_akdes))
only_late <- which(!names(late_akdes) %in% names(first_akdes))

akdes_first <- lapply(1:length(first_akdes[only_first]), function(x){
  # do the same transformations to the juvenile only akdes as done to the combination above
  f <- raster::raster(first_akdes[only_first][[x]], DF = "PDF")
  # using the terra package for reprojection
  f <- terra::rast(f)
  f <- terra::project(f, "EPSG:4326")
  # adding the date as the name
  names(f) <- names(first_akdes[only_first])[x]
  return(f)
})

akdes_late <- lapply(1:length(late_akdes[only_late]), function(x){
  # do the same transformations to the adult only akdes as done to the combination above
  f <- raster::raster(late_akdes[only_late][[x]], DF = "PDF")
  # using the terra package for reprojection
  f <- terra::rast(f)
  f <- terra::project(f, "EPSG:4326")
  # adding the date as the name
  names(f) <- names(late_akdes[only_late])[x]
  return(f)
})

akdes <- c(akdes, akdes_first)#, akdes_late)
# saveRDS(akdes, file = "C:/Users/hbronnvik/Documents/storkSSFs/akdes_combined_2023-02-03.rds")

# plot all of the akdes to see how they move through time  
# p <- vect(outlines)
# png(file="%02d.png", width=480, height=480)
# 
# for (i in 1:length(akdes)){
#   par(mar = c(0, 0, 0, 0))
#   ud <- akdes[[i]]
#   ud <- subst(ud, 0, NA)
#   terra::plot(p, main = names(akdes[[i]]))
#   terra::plot(ud, add = T, buffer = T)
#   terra::plot(p, add = T)
# }
# dev.off()
# 
# library(av)
# library(gtools)
# imgs <- list.files("C:/Users/hbronnvik/Documents/storkSSFs/jan_ras_plot", full.names = TRUE)
# av_encode_video(mixedsort(sort(imgs)), framerate = 115, output = "test2.mp4")

library(terra)
terra::plot(late_akdes[[1]])
projection(late_akdes[[1]])
check <- raster(late_akdes[[1]], DF = "PDF")
check <- rast(check)
check <- terra::project(check, "EPSG:4326")
p <- vect(outlines)
check <- akdes[[6000]]
check <- akdes[[1904]]
check <- subst(check, 0, NA)
terra::plot(p, main = names(akdes[[1904]]))
terra::plot(check, add = T, buffer = T)
terra::plot(p, add = T)

f <- raster(first_akdes[[1]], DF = "PDF")
f <- terra::rast(f)
f <- terra::project(f, "EPSG:4326")
l <- raster(late_akdes[[1]], DF = "PDF")
l <- terra::rast(l)
l <- terra::project(l, "EPSG:4326")
l <- terra::project(l, f)
fl <- terra::merge(f, l)


locs <- locs %>% 
  arrange(track_id)
start_time <- Sys.time()
interpolated_tracks <- lapply(split(locs, locs$track_id), function(x){
  
  # print(paste0("Interpolating the track: ", x, ", ", which(unique(fall_tracks$track_id) == x),
  #              " of ", length(unique(fall_tracks$track_id)), "."))
  print(paste0("Interpolating the track: ", unique(x$track_id), ", ", which(unique(locs$track_id) == unique(x$track_id)),
               " of ", length(unique(locs$track_id)), "."))
  
  # ind <- fall_tracks %>%
  #   filter(track_id == x)
  
  ind <- as.telemetry(x)#, keep = T)
  
  GUESS1 <- ctmm.guess(ind, interactive = FALSE)
  
  print("fitting model")
  
  FIT1_pHREML <- ctmm.select(ind, GUESS1, method = "pHREML", verbose = TRUE)
  
  print("filling the gaps")
  
  filled <- predict(ind, CTMM = FIT1_pHREML[[1]], res = .25)#quarter hour
  
  cn <- colnames(filled)
  
  filled <- as.data.frame(filled@.Data)
  
  colnames(filled) <- cn
  
  filled$track_id <- unique(x$track_id)
  
  return(filled)
}) #%>% reduce(rbind)
Sys.time()-start_time






# pull out autumn data
fall_tracks <- tracks %>% 
  filter(phase == "autumn_migration")  %>% 
  # remove erroneous first fix of the day (only a half hour earlier than the next location)
  group_by(individual.id, year(timestamp), date(timestamp)) %>% 
  slice(-1) %>% 
  ungroup() %>% 
  # get the tracks in order of their ids
  arrange(track_id) %>% 
  group_by(track_id) %>%
  mutate(locs = n()) %>%
  ungroup() %>%
  filter(locs > 99)

# reduce the data to check results
fall_tracks <- fall_tracks %>%
  filter(!track_id %in% unique(fall_tracks$track_id)[101:length(unique(fall_tracks$track_id))])

## interpolate the missing hours in the tracks of each individual

# only 20 birds on their autumn migrations in August and September
locs <- locations_thin %>% 
  filter(month(datestamp) %in% c(8,9))


# fit the akde for each individual
# predict the missing locations using the total time from start to end in hourly intervals
# fit akde for each hour using one location for each individual
start_time <- Sys.time()
interpolated_tracks <- lapply(split(fall_tracks, f = factor(fall_tracks$track_id)), function(x){
  
  # print(paste0("Interpolating the track: ", x, ", ", which(unique(fall_tracks$track_id) == x),
  #              " of ", length(unique(fall_tracks$track_id)), "."))
  print(paste0("Interpolating the track: ", unique(x$track_id), ", ", which(unique(fall_tracks$track_id) == unique(x$track_id)),
               " of ", length(unique(fall_tracks$track_id)), "."))
  
  # ind <- fall_tracks %>%
  #   filter(track_id == x)
  
  ind <- as.telemetry(x)#, keep = T)
  
  GUESS1 <- ctmm.guess(ind, interactive = FALSE)
  
  print("fitting model")
  
  FIT1_pHREML <- ctmm.select(ind, GUESS1, method = "pHREML", verbose = TRUE)
  
  print("filling the gaps")
  
  filled <- predict(ind, CTMM = FIT1_pHREML[[1]], dt = 3600)
  
  cn <- colnames(filled)
  
  filled <- as.data.frame(filled@.Data)
  
  colnames(filled) <- cn
  
  filled$track_id <- unique(x$track_id)
  
  return(filled)
}) #%>% reduce(rbind)
Sys.time()-start_time
#sk <- interpolated_tracks
#save(interpolated_tracks, file = "interpolated_tracks_100inds_21.09.22.RData")

interpolated_tracks$t <- as.POSIXct(interpolated_tracks$t, tz = "UTC", origin = '1970-01-01')

wgs <- "+proj=tpeqd +lat_1=38.4747235276509 +lon_1=-4.10793186111652 +lat_2=45.0136355260082 +lon_2=4.8721292939578 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

i_tracks_mv <- lapply(split(interpolated_tracks, f = interpolated_tracks$track_id), function(df){
  
  mv <- move(x = df$x, y = df$y, time = df$t, proj = wgs, data = df)
  
})

## make stack from filtered data
i_ms <- moveStack(i_tracks_mv, force_tz = TRUE)
#i_ms <- spTransform(fall_ms, center = T)

# # begin with a single id to check run
# ind <- fall_tracks %>% 
#   filter(track_id == unique(fall_tracks$track_id)[1])
# 
# inds <- as.telemetry(ind)
# 
# plot(inds,col=rainbow(length(inds)))
# 
# GUESS1 <- ctmm.guess(inds, interactive = FALSE)
# 
# start_time <- Sys.time()
# FIT1_pHREML <- ctmm.select(inds, GUESS1, method = "pHREML", verbose = TRUE)
# Sys.time() - start_time # Time difference of 1.630057 mins
# summary(FIT1_pHREML)

# SVF <- variogram(inds)
# plot(SVF, CTMM = FIT1_pHREML[[1]],
#      units = TRUE, fraction = 0.5, level = c(0.95, 0.50), 
#      col = "black", col.CTMM = "red")

# AKDE1_pHREML <- akde(inds, FIT1_pHREML, debias = TRUE)
# summary(AKDE1_pHREML, level.UD = 0.95)$CI

# newEXT <- extent(AKDE1_pHREML)
# plot(inds, UD = AKDE1_pHREML, ext = newEXT)
# title(expression("AKDEc"))
# 
# check <- predict(inds, CTMM = FIT1_pHREML[[1]])

# OD <- occurrence(DATA, FITS[[1]])
# SIM <- simulate(DATA, FITS[[1]], dt = 5 %#% "min")
# # should allow an occurrence distribution which sort of smooths over the tracks, then can be 
# # related to underlying environmental variables to see which covariates were sampled by the animal
# # (by multiplying PDF OD (probability of the animal in that cell) with a get value raster layer)
# # over the sampling time regardless of gappy data and not overweighting oversampled locations


## build kdes for each hour
check <- interpolated_tracks 
check$datestamp <- check$t
year(check$datestamp) <- 1995
check$datestamp <- round_date(datestamp)

length(unique(check$datestamp))

ind <- as.telemetry(check)
GUESS <- ctmm.guess(ind, interactive = FALSE)
FIT_pHREML <- ctmm.select(ind, GUESS, method = "pHREML", verbose = TRUE)
AKDE <- akde(ind, FIT_pHREML)








# fresh data
ind_locs <- getMovebankData(study = 24442409, animalName =  24450590, sensorID = "GPS", 
                            login = loginStored, removeDuplicatedTimestamps = T)
unique(year(ind_locs$timestamp))

# the distance covered in each day for the full data (bursts and all)
locs_df <- ind_locs %>% 
  as.data.frame() %>% 
  filter(year(timestamp) == 2015) %>% 
  group_by(date(timestamp)) %>% 
  mutate(distance = distVincentyEllipsoid(c(head(location_long,1), head(location_lat, 1)), 
                                          c(tail(location_long,1), tail(location_lat, 1)))) %>% 
  ungroup()

# the first time after August that it moved at least 100 km . day
on <- locs_df %>% 
  filter(timestamp > as.POSIXct("2015-08-01 00:00:00", tz = "UTC", origin = "1970-01-01") & distance > 100*1000) %>% 
  slice(1) %>% 
  select(timestamp)

# the last time after August that it moved at least 100 km / day
off <- locs_df %>% 
  filter(timestamp > as.POSIXct("2015-08-01 00:00:00", tz = "UTC", origin = "1970-01-01") & distance > 100*1000) %>% 
  slice(n()) %>% 
  select(timestamp)

# filter the data down to between on and off of migration
locs_df <- locs_df %>% 
  filter(timestamp > on & timestamp < off)



load("C:/Users/hbronnvik/Documents/loginStored.rdata")
studies <- c(24442409, 212096177, 76367850, 21231406, 1176017658, 173641633)

## determine the identities of the nestling birds (remove any care center adults)
nestlings <- lapply(studies, function(x){
  
  print(x)
  
  birds <- getMovebankReferenceTable(study = x, login = loginStored) %>% 
    drop_na(animal_id) %>% 
    filter(sensor_type_id == 653)
  
  if(length(unique(birds$animal_id[grep("juv|chick", birds$animal_comments, fixed = F)])) > 0){
    nesties1 <- unique(birds$animal_id[grep("juv|chick", birds$animal_comments, fixed = F)])
  } else {nesties1 <- NA}
  
  if("animal_life_stage" %in% colnames(birds) & length(unique(birds$animal_id[which(birds$animal_life_stage == "nestling")])) > 0){
    nesties2 <- unique(birds$animal_id[which(birds$animal_life_stage == "nestling")])
  } else {nesties2 <-  NA}
  
  nesties <- na.omit(unique(c(unique(nesties1), unique(nesties2))))
  
  return(nesties)
  
}) %>% unlist()


birds <- lapply(studies, function(x){
  birds <- getMovebankAnimals(study =  x, login = loginStored) %>% 
    filter(sensor_type_id == 653) %>% # even birds that die provide information before that (and possibly during)
    dplyr::select(individual_id) %>% 
    mutate(study_id = x)
  return(birds)
}) %>% reduce(rbind)

birds <- birds[which(birds$individual_id %in% nestlings),]

subset <- birds[1:100,]

start_time <- Sys.time()
fall_tracks <- lapply(split(subset, f = subset$individual_id), function(x){
  
  # fresh data
  ind_locs <- getMovebankData(study = x$study_id, animalName =  x$individual_id, sensorID = "GPS", 
                              login = loginStored, removeDuplicatedTimestamps = T)
  # unique(year(ind_locs$timestamp))
  
  # the distance covered in each day for the full data (bursts and all)
  locs_df <- ind_locs %>% 
    as.data.frame() %>% 
    filter(year(timestamp) == unique(year(timestamp))[1]) %>% 
    group_by(date(timestamp)) %>% 
    mutate(distance = distVincentyEllipsoid(c(head(location_long,1), head(location_lat, 1)), 
                                            c(tail(location_long,1), tail(location_lat, 1)))) %>% 
    ungroup()
  
  # the first time after August that it moved at least 100 km/day
  on <- locs_df %>% 
    filter(timestamp > as.POSIXct("2015-07-01 00:00:00", tz = "UTC", origin = "1970-01-01") & distance > 100*1000) %>% 
    slice(1) %>% 
    dplyr::select(timestamp)
  
  # the last time after August that it moved at least 100 km/day
  off <- locs_df %>% 
    filter(timestamp > as.POSIXct("2015-07-01 00:00:00", tz = "UTC", origin = "1970-01-01") & distance > 100*1000) %>% 
    slice(n()) %>% 
    dplyr::select(timestamp)
  
  # filter the data down to between on and off of migration
  locs_df <- locs_df %>% 
    filter(timestamp > on[,1] & timestamp < off[,1]) %>% 
    mutate(seq15 = round_date(timestamp, "15 minutes")) %>% 
    group_by(seq15) %>% 
    slice(1)
  
  return(locs_df)
}) #%>% reduce(rbind)
Sys.time() - start_time
#save(fall_tracks, file = "fall_tracks_byOnOff_15min100ind.RData")
#load("fall_tracks_byOnOff_15min100ind.RData")

fall_tracks <- lapply(fall_tracks, function(x){
  x <- x[, c("location_long","location_lat", "ground_speed", "height_above_ellipsoid",
             "timestamp", "individual_id", "distance", "seq15")]
  
  return(x)
}) %>% reduce(rbind)

# the first location in a given 15 minute interval (rm the bursts)
thin_tracks <- fall_tracks %>% 
  #   mutate(seq15 = round_date(timestamp, "15 minutes")) %>% 
  #   group_by(individual_id, seq15) %>% 
  #   slice(1) %>% 
  #   ungroup() %>% 
  # remove erroneous first fix of the day (only 5 minutes earlier than the next location)
  group_by(individual_id, date(timestamp)) %>% 
  slice(-1)

# define birds that took eastern routes as ones that are ever east of 12 longitude
eastern_birds <- unique(thin_tracks$individual_id[thin_tracks$location_long > 12])

# remove the eastern birds and add a phase column
thin_tracks <- thin_tracks %>% 
  filter(!individual_id %in% eastern_birds)

countries <- c("Niger", "Denmark", "Western Sahara", "Germany", "Spain", "Morocco", "Portugal", "Belgium", "Netherlands",
               "Algeria", "Liechtenstein", "Luxembourg", "Mauritania", "Mali", "Senegal", 
               "Gambia", "Ivory Coast", "Ghana", "Burkina Faso", "Switzerland", "France", "Benin")
ggplot(thin_tracks, aes(x=location_long, y=location_lat)) + 
  borders(region = countries, fill = "gray80") + geom_point()+ facet_wrap(~individual_id)+theme_minimal()

thin_tracks <- thin_tracks %>% 
  mutate(alignment = paste(month(seq15), day(seq15), hour(seq15), minute(seq15), sep = "_"))

thin_tracks$alignment <- thin_tracks$seq15
year(thin_tracks$alignment) <- 1995

check <- thin_tracks %>% 
  group_by(date(alignment)) %>% 
  summarise(rown = n()) %>% 
  group_by(rown) %>% 
  summarise(stamps = n())
hist(check$rown, breaks = nrow(check))


# the number of locations per day (one ID at each time stamp) is much higher than per 15 mins
# this may improve with the addition of the remaining 269 birds?

# fit the akde for each individual
# predict the missing locations using the total time from start to end in hourly intervals
# fit akde for each hour using one location for each individual
thin_tracks <- thin_tracks %>% 
  mutate(track_id = paste(individual_id, year(timestamp), sep = "_")) %>% 
  arrange(track_id)
start_time <- Sys.time()
interpolated_tracks <- lapply(split(thin_tracks, f = factor(thin_tracks$track_id)), function(x){
  
  print(paste0("Interpolating the track: ", unique(x$track_id), ", ", which(unique(thin_tracks$track_id) == unique(x$track_id)),
               " of ", length(unique(thin_tracks$track_id)), "."))
  
  # sdf <- split(thin_tracks, f = factor(thin_tracks$track_id))
  # x <- sdf[[3]]
  
  ind <- as.telemetry(x)#, keep = T)
  
  GUESS1 <- ctmm.guess(ind, interactive = FALSE)
  
  print("fitting model")
  
  FIT1_pHREML <- ctmm.select(ind, GUESS1, method = "pHREML", verbose = TRUE)
  
  print("filling the gaps")
  
  filled <- predict(ind, CTMM = FIT1_pHREML[[1]], dt = 3600)
  
  cn <- colnames(filled)
  
  filled <- as.data.frame(filled@.Data)
  
  colnames(filled) <- cn
  
  filled$track_id <- unique(x$track_id)
  
  return(filled)
}) #%>% reduce(rbind)
Sys.time()-start_time

# install.packages("maptools")
library(maptools)
data(wrld_simpl)
ws <- crop(wrld_simpl, extent(-20,20,0,60))
# plot(ws)
tmax <- raster::getData('worldclim', var = "tmax", res = 10)
mask <- crop(raster(tmax, 1), extent(ws))
raster::plot(mask)
mm <- buffer(mask, 30000)
raster::plot(mm)
# writeRaster(mm, "buffered_water_ras.grd", overwrite = T)
mm <- raster("~/storkSSFs/buffered_water_ras.grd")
outlines <- rasterToPolygons(mm, dissolve=TRUE)

# each time step gets its own kde
thin_tracks <- tracks %>% 
  mutate(hourly = round_date(timestamp, "hour")) %>% 
  group_by(individual.id, hourly) %>% 
  slice(1) %>% 
  ungroup() %>% 
  mutate(alignment = round_date(timestamp, "minute"))
year(thin_tracks$alignment) <- 1995
kde_tracks <- split(thin_tracks, thin_tracks$alignment)
kde_tsub <- kde_tracks[sapply(kde_tracks, function(x) nrow(x) > 30)]
start_time <- Sys.time()
kdes <- lapply(kde_tsub, function(ind){
  print(unique(ind$alignment))
  #ind <- kde_tsub[[1]] # a testing line
  
  ind$individual.id <- "a_name"
  ind$individual.local.identifier <- "a_name"
  
  ind <- as.telemetry(ind)#, keep = T)
  
  projection(ind) <-  "ESRI:54009" # mollwiede projektion to keep all the data map-able and comparable
  
  fit <- ctmm.fit(ind) # IID
  KDE1 <- akde(ind, fit, SP = outlines)
  #plot(ind, UD = KDE1, ext = extent(KDE1))
  return(KDE1)
})
Sys.time() - start_time
#save(kdes, file = "kdes_ls_22.9.RData")
raster::plot(kdes[[200]])

library(rgdal)
library(raster)
library(terra)

srtm_mosaic <- raster::raster("srtm_mosaic.grd")
projection(srtm_mosaic) <- "ESRI:54009"
binary <- srtm_mosaic
values(binary[!is.na(binary)]) <- 1
binary[is.na(srtm_mosaic)] <- 0 
plot(binary) #"+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs"
sp_land <- rasterToPoints(srtm_mosaic, spatial = TRUE)

## download DEM data for the countries along the flyway 
countries <- c("Niger", "Denmark", "Western Sahara", "Germany", "Spain", "Morocco", "Portugal", "Belgium", "Netherlands",
               "Algeria", "Liechtenstein", "Luxembourg", "Mauritania", "Mali", "Senegal", 
               "Gambia", "Ivory Coast", "Ghana", "Burkina Faso", "Switzerland", "France", "Benin") %>% 
  countrycode(origin = 'country.name', destination = 'iso3c')

shp <- shapefile("C:/Users/heste/Desktop/HB_storks/srtm/tiles.shp")
plot(shp)

#Get country geometry first
outlines <- lapply(countries, function(x){
  print(x)
  temp <- getData("GADM", country = x, level=0)
  crs(temp) <- crs(shp)
  
  return(temp)
}) %>% reduce(rbind)
ws <- crop(wrld_simpl, extent(outlines))
regions <- raster(extend(extent(ws), c(1000, 1000)), res = 1000, crs = CRS("+proj=longlat +datum=WGS84 +no_defs"), vals = 0)
#ws@data[,1] <- runif(nrow(ws))
plot(ws, col = ws$REGION)
check <- raster::rasterize(ws, regions, field = "FIPS", background = 0, fun = "first")
plot(check)


tmax <- raster::getData('worldclim', var = "tmax", res = 10)
mask <- crop(raster(tmax, 1), extent(ws))
raster::plot(mask)
mm <- buffer(mask, 60000)
raster::plot(mm)
#writeRaster(mm, "C:/Users/heste/Desktop/HB_storks/buffered_water_ras.grd")

outlines <- rasterToPolygons(mm, dissolve=TRUE)



uds <- lapply(kdes, function(x){
  hr <- ctmm::SpatialPolygonsDataFrame.UD(x, proj4string = datproj) %>%
    sp::spTransform(sp::CRS("+proj=longlat")) %>%
    fortify() %>%
    bind_rows()
  return(hr)
})

# plug in the data and the column names
ggplot(uds, aes(x = location.long, y = location.lat)) + 
  # if you want borders, change the fill shade
  borders(regions = countries, fill = "gray50") +
  # this is an alternative to listing all the countries, but I find it behaves oddly:
  #borders(database = "world", xlim = c(-10, 20), ylim = c(10, 60), fill = "gray50") + 
  # I don't know whether you would want to color by ID, year, season, etc.
  geom_polygon(aes(color = year(timestamp))) +  
  # optional to change colors, but requires the viridis library:
  #scale_fill_viridis(discrete=TRUE, option="A") + 
  # fix the labels
  labs(x = " Longitude", y = "Latitude")+
  # clean up the gray background
  theme_classic() +
  # remove the legend
  theme(legend.position="none", axis.text = element_text(color = "black")) # +
# and if you want to animate that:
#transition_time(time = timestamp)
#transition_states(states = year(timestamp))

datproj <- sp::CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs")
# Build images -> save them at .png format
png(file="%02d.png", width=480, height=480)

for (i in 1:length(kdes)){
  hr <- ctmm::SpatialPolygonsDataFrame.UD(kdes[[i]], proj4string = datproj) %>%
    sp::spTransform(sp::CRS("+proj=longlat"))
  #hr <- uds[[i]]
  par(mar = c(0, 0, 0, 0))  
  sp::plot(flyway, col="grey80", border= F, ylim = c(-11,65)) 
  sp::plot(hr, add=T)
  #mtext(paste(day(dates[i]), "-", month(dates[i])), side=3)
}
dev.off()

library(magick)
imgs <- list.files("C:/Users/heste/Desktop/HB_storks/kdegif/", full.names = TRUE)
img_list <- lapply(gtools::mixedsort(imgs), image_read)

## join the images together
img_joined <- image_join(img_list)

## animate at 2 frames per second
img_animated <- image_animate(img_joined, fps = 100)

image_write(image = img_animated,
            path = "C:/Users/heste/Desktop/HB_storks/kernels_ctmm_15min.gif")




Lines_ws <- lapply(split(fall_tracks,fall_tracks$track_id), function(x){
  points <- sp::SpatialPoints(cbind(x$long, x$lat), proj4string = wgs)
  sp_line <- as(points,"SpatialLines")
})


# bltr
par(mar = c(0, 0, 0, 0))  
sp::plot(flyway, col="grey70", border= F, ylim = c(-11,65)) 

#add latitude lines
lines(x = c(-20,45.9), y = c(0,0),lty = 2,lwd = 1, col = "grey50")
lines(x = c(-20,45.9), y = c(30,30),lty = 2,lwd = 1, col = "grey50")
lines(x = c(-20,45.9), y = c(60,60),lty = 2,lwd = 1, col = "grey50")
text(x = -17, y = 30.7, "30° N", col = "grey50", cex = 1.1, font = 3)
text(x = -17, y = 60.7, "60° N", col = "grey50", cex = 1.1, font = 3)
text(x = -17, y = 0.70, "0° N", col = "grey50", cex = 1.1, font = 3)

invisible(lapply(1:length(Lines_ws), function(x){
  l <- Lines_hb[[x]]
  lines(l, lty= 1, lwd =  2.5, col= rainbow)
}))

### 11.05.2022
library(maptools)
data("wrld_simpl")
ws <- crop(wrld_simpl, extent(-20,20,0,60))
raster::plot(ws)
raster::plot(raster(kdes[[1]]))
ws <- spTransform(ws, crs(raster(kdes[[1]])))
raster::plot(ws, add = T)

build <- raster(kdes[[1]])
build[build > 0.95] <- NA
raster::plot(build, col = heat.colors(8, alpha = 1), legend = F)
raster::plot(ws, add = T)

tracks <- tracks %>% 
  mutate(datestamp = timestamp) 

year(tracks$datestamp) <- 1995

kde_tsub[[1]]$alignment[1]

tsub_stamps <- lapply(1:length(kde_tsub), function(x){
  stamp <- kde_tsub[[x]]$alignment[1]
  return(stamp)
}) %>% 
  unlist() %>% 
  as.POSIXct(tz = "UTC", origin = "1970-01-01")


bt1 <- read.csv("bt1.csv")
bt2 <- read.csv("bt2.csv")
bt3 <- read.csv("bt3.csv")

burst_tracks <- rbind(bt1, bt2, bt3)

burst_tracks <- burst_tracks %>% 
  mutate(timestamp = as.POSIXct(gsub("\\.000", "", timestamp), tz = "UTC"), 
         alignment = round_date(timestamp, "minute"),
         location.long = as.numeric(location.long),
         location.lat = as.numeric(location.lat))
year(burst_tracks$alignment) <- 1995

t <- lapply(tsub_stamps, function(x){
  # print(x)
  df <- burst_tracks %>% 
    filter(alignment == x)
  
  if(nrow(df > 3)){
    # for each unique hour in the data, call the list item
    id <- which(tsub_stamps == x)
    kd <- raster(kdes[[id]])
    
    # make the data an SP object and transform to the same CRS as the AKDE
    spdf <- df
    coordinates(spdf) <- ~location.long + location.lat
    crs(spdf) <- crs("+proj=longlat +datum=WGS84 +no_defs +type=crs")
    spdf <- spTransform(spdf, crs(kd))
    
    # extract the value for each location
    vals <- raster::extract(kd, spdf)#[, c("location.long", "location.lat")])
    
    # append
    df <- df %>% 
      mutate(UD = vals)
    return(df)
  }
  
  
}) %>% reduce(rbind)

t <- t %>% 
  mutate(season = str_extract(track_id, "fall|spring"),
         individual_id = str_extract(track_id, "[^_]+"))

# number the migrations according to phase
sorting <- t  %>% 
  group_by(individual_id, season) %>% 
  count(track_id) %>% 
  mutate(no. = row_number(),
         number = ifelse(no. == 1, "first", 
                         ifelse(no. == 2, "second", 
                                ifelse(no. == 3, "third", 
                                       ifelse(no. == 4, "fourth", 
                                              ifelse(no. == 5, "fifth", 
                                                     ifelse(no. == 6, "sixth", 
                                                            ifelse(no. == 7, "seventh", 
                                                                   ifelse(no. == 8, "eighth", 
                                                                          ifelse(no. == 9, "ninth", "tenth")))))))))) %>% 
  ungroup()

# add those numbers onto the full data set
t <- full_join(t, sorting) %>% 
  mutate(migration = paste0(number, "_", season)) 
t$number <- factor(t$number, levels = c("first", "second", 'third', "fourth", "fifth", "sixth", "seventh", "eighth", "ninth"))

ggplot(t, aes(number, UD)) + 
  geom_boxplot() +
  theme_classic() +
  facet_wrap(~case_)

library(survival)
modelUD <- function(df) {
  clogit(case_ ~ scaled_UD + scaled_sl + scaled_ta + strata(step_id_), data = df)
}

SSF_results <- t %>%
  mutate(scaled_UD = scale(UD),
         scaled_sl = scale(sl_),
         scaled_ta = scale(ta_)) %>% 
  group_by(track_id) %>%
  filter(length(unique(UD)) > 3) %>% 
  nest() %>%
  mutate(ssf_modelUD = purrr::map(data, modelUD),
         ssf_coefsUD = purrr::map(ssf_modelUD, coef),
         AIC_TV = map_dbl(ssf_modelUD, ~AIC(.)))

# flatten the coefficient column
ssf_coefs <- unnest(SSF_results, ssf_coefsUD) %>% 
  group_by(track_id) %>% 
  # take the UD coefs
  slice(1) %>% 
  ungroup() %>% 
  mutate(individual_id = str_extract(track_id, "[^_]+")) %>% 
  full_join(sorting) %>% 
  group_by(number) %>% 
  mutate(samp = length(unique(track_id))) %>% 
  ungroup() %>% 
  mutate(migration = paste(number, season, sep = "_"))

library(EnvStats)
ggplot(ssf_coefs[which(abs(ssf_coefs$ssf_coefsUD) < 40),], aes(as.factor(migration), ssf_coefsUD, group = migration, fill = season, label = number)) +
  geom_boxplot() +
  labs(x = "Migration", y = "Coefficient") +
  theme_classic() +
  theme(text = element_text(size = 20),
        axis.text = element_text(colour = "black"))  +
  stat_n_text()

