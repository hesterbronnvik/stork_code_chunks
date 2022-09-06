### Kernel density estimate for storks migrating across western Europe
### Hester Bronnvik
### 18.07.2022

library(move)
library(scales)
library(ks)
library(adehabitatHR)
library(lubridate)
library(tidyverse)
library(gganimate)
library(leaflet)
library(viridis)


### retrieve the data
## this requires downloading large studies off movebank

# set wd and load login data
setwd("C:/Users/heste/Desktop/HB_storks")
load("loginStored.rdata")

# collect the names of the stork studies to use
studies <- c(24442409, 212096177, 76367850, 21231406,1176017658, 173641633)

start_time <- Sys.time()
locations <- lapply(studies, function(x){
  
  print(paste0("Study ", x, "."))
  # get the names of the animals to download
  birds <- getMovebankAnimals(study =  x, login = loginStored)
  # remove birds that had no data 
  bird_names <- unique(birds$individual_id[which(birds$number_of_events > 0)])
  
  # remove the names of birds that have already been downloaded
  #bird_names <- bird_names[!bird_names %in% sub(".rds", "", list.files("C:/Users/heste/Desktop/HB_storks/ind_locations/"))]
  
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
    
    # if internet is unstable, saving per individual may be necessary
    #saveRDS(locations_thin, file = paste0("C:/Users/heste/Desktop/HB_storks/ind_locations/", y, ".rds"))
    
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

# if internet is unstable, read in each saved individual file and bind
# ind_locations <- list.files("C:/Users/heste/Desktop/HB_storks/ind_locations/", full.names = T)
# 
# locations <- lapply(ind_locations, function(x){
#   ind <- readRDS(x)
# }) %>% reduce(rbind) # 46 million rows of data, a lengthy process

# study_locations <- list.files("C:/Users/heste/Desktop/HB_storks/study_locations/", full.names = T)
# 
# locations2 <- lapply(study_locations, function(x){
#    ind <- readRDS(x)
# }) %>% reduce(rbind) # 7.3 million rows of data

# locations_full <- rbind(locations, locations2)

# split the full data set up by month. They are all in the same year now, 
# so this provides 5 files just to reduce time reading in the 53 million rows of data. 
# If the whole workflow is done at once, this is unnecessary.
lapply(unique(month(locations_full$datestamp)), function(x){
  mo_locs <- locations_full %>% 
    dplyr::filter(month(locations_full$datestamp) == x)
  
  #save(mo_locs, file = paste0("C:/Users/heste/Desktop/HB_storks/months_6_studies/month_", x, "_locations.RData"))
})

# if the data were split, read them in and thin them
months_files <- list.files("C:/Users/heste/Desktop/HB_storks/months_6_studies/", full.names = T)

locations_thin <- lapply(months_files, function(x){
  
  load(x)
  
  # thin the data to 5 minute intervals
  locations_thin <- mo_locs %>% 
    group_by(individual.id) %>% 
    mutate(dt_5min = round_date(timestamp, "5 minutes")) %>% 
    group_by(dt_5min) %>% 
    slice(1) %>%  
    ungroup()
  
  return(locations_thin)
  
}) %>% reduce(rbind)

### build a kernel estimate for each subset

# assign each unique day a group so that kernels are built for every three days
dates <- data.frame(date = unique(date(locations_thin$datestamp)),
                    group = rep(1:51, times = c(3, 3, 3, 3, 3, 3, 3, 3, 
                                                   3, 3, 3, 3, 3, 3, 3, 3, 
                                                   3, 3, 3, 3, 3, 3, 3, 3, 
                                                   3, 3, 3, 3, 3, 3, 3, 3, 
                                                   3, 3, 3, 3, 3, 3, 3, 3, 
                                                   3, 3, 3, 3, 3, 3, 3, 3, 
                                                   3, 3, 3)))


dat_2 <- SpatialPointsDataFrame(locations_thin[,c("location.long", "location.lat")], locations_thin[,c(1,4,9,10,16)])
kus <- kernelUD(dat_2[, 4], same4all = TRUE)

data(puechabonsp)

## adehabitat library
kernels <- lapply(unique(date(locations_thin$datestamp)), function(f){
  print(f)
  data_sp <- locations_thin %>%  # store the data, then create a spatial object for plotting
    dplyr::filter(date(datestamp) == f) %>% 
    dplyr::select(location.long, location.lat)
  
  coordinates(data_sp) <- ~ location.long + location.lat
  
  #kern <- Hpi(x = coordinates(data_sp))
  
  kern <- kernelUD(data_sp, h = "href")# the ad hoc method is used for the smoothing parameter
  
  # fhat <- kde(x = data_sp, H = get("kern"), compute.cont = T)
  # 
  # png(filename= paste0("C:/Users/heste/Desktop/HB_storks/gif/", f, ".png"))
  # plot(fhat, xlab = "Longitude", ylab = "Latitude", xlim=c(-15,10), ylim = c(25, 50))
  # dev.off()
  
  return(kern)
})
image(getvolumeUD(kernels[[1]]), col = rev(viridis(15)))
#df <- as.data.frame(getvolumeUD(kernels[[1]]))
image(kernels[[1]])
plot(getverticeshr(kernels[[1]], standardize = TRUE, 95), add = TRUE)

plot(getverticeshr(kernels[[1]]))
data_sp <- locations_thin[which(date(locations_thin$datestamp) == unique(date(locations_thin$datestamp))[1]),]
coordinates(data_sp) <- ~ location.long + location.lat
points(data_sp, cex = 1, pch = 16)

# ud1 <- as(kernels[[1]], "SpatialPolygonsDataFrame")
# leaflet() %>% addTiles() %>% addPolygons(data = ud1)



ud95 <- lapply(kernels, function(x){
  getverticeshr(x, 95)
  }) #%>% reduce(rbind)

sapply(1:length(ud95), function(i) {
  
  #row.names(ud[[i]]) <- unique(date(locations_thin$datestamp))[i]
  row.names(ud[[i]]@data) <<- unique(date(locations_thin$datestamp))[i]
  
}) 

sdf_poly95 <- Reduce(rbind, ud95)
df95 <- fortify(sdf_poly95)

leaflet(sdf_poly95) %>% addTiles() %>%
  addPolygons(weight = 1, fillOpacity = .2)

ud50 <- lapply(kernels, function(x){
  getverticeshr(x, 50)
}) #%>% reduce(rbind)

sapply(1:length(ud50), function(i) {
  
  #row.names(ud[[i]]) <- unique(date(locations_thin$datestamp))[i]
  row.names(ud[[i]]@data) <<- unique(date(locations_thin$datestamp))[i]
  
})

sdf_poly50 <- Reduce(rbind, ud50)
df50 <- fortify(sdf_poly50)

g <- ggplot(df, aes(x = long, y = lat, fill = id, group = group)) +
  geom_polygon(alpha = .4) +
  coord_equal() +
  theme_void()+
  theme(legend.position="none", axis.text = element_text(color = "black"))
g
g + facet_wrap(~id)

mv <- df %>% 
  group_by(id) %>% 
  mutate(timing = 1:153)

mv95 <- data.frame()

for (i in 1:length(unique(df95$id))) {
  x <- df95[df95$id == unique(df95$id)[i],]
  x$date <- unique(date(locations_thin$datestamp))[i]
  mv95 <- rbind(mv95, x)
}

mv50 <- data.frame()

for (i in 1:length(unique(df50$id))) {
  x <- df50[df50$id == unique(df50$id)[i],]
  x$date <- unique(date(locations_thin$datestamp))[i]
  mv50 <- rbind(mv50, x)
}

ggplot(mv, aes(x = long, y = lat, fill = id, group = group)) + 
  geom_polygon(alpha = 0.4) +
  labs(x = " Longitude", y = "Latitude", title = "Date: {previous_state}") +
  transition_states(states = group)+
  theme_classic() +
  theme(legend.position="none", axis.text = element_text(color = "black"))

countries <- c("Germany", "Switzerland", "France", "Spain", "Morocco", "Italy", "Portugal", "Belgium", "Netherlands",
               "Tunisia", "Algeria", "Liechtenstein", "Luxembourg", "Mauritania", "Mali", "Libya", "Austria", "Slovenia",
               "Czech Republic", "Western Sahara", "Senegal", "Gambia", "Guinea Bissau", "Guinea", "Mali", "Sierra Leone",
               "Liberia", "Ivory Coast", "Togo", "UK", "Ghana", "Burkina Faso", "Benin", "Niger", "Nigeria", "Chad", "Cameroon")

ggplot(mv, aes(x = long, y = lat, fill = id, group = group)) + 
  borders(regions = countries, fill = "gray50") +
  geom_polygon(alpha = 0.4) +
  labs(x = " Longitude", y = "Latitude") +
  transition_time(time = date)+
  theme_classic() +
  theme(legend.position="none", axis.text = element_text(color = "black"))

img_animated <- ggplot(mv95, aes(x = long, y = lat, fill = group, group = group)) + 
  borders(regions = countries, fill = "gray70", colour = "gray70") +
  #borders(database = "world", xlim = c(-10, 20), ylim = c(10, 60), fill = "gray80", colour = "gray80") +
  geom_polygon(alpha = 0.4) +
  geom_polygon(data = mv50, alpha = 0.8) +
  #geom_point(data = locations_thin, aes(x = location.long, y = location.lat, group = date(datestamp))) +
  scale_fill_viridis(discrete=TRUE, option="A") +
  labs(x = " Longitude", y = "Latitude", title = 'Date: {gsub("1995-", "", current_frame)}')+
  transition_manual(frames = date)+
  theme_classic() +
  theme(legend.position="none", axis.text = element_text(color = "black"))
anim_save("C:/Users/heste/Desktop/HB_storks/kernels_multi.gif", img_animated)

library("rnaturalearth")
library("rnaturalearthdata")

world <- ne_countries(scale = "medium", returnclass = "sf")

ggplot(data = world) +
  geom_sf(fill = "gray70", colour = "gray70") +
  coord_sf(xlim = c(-30, 30), ylim = c(0, 60), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") + 
  theme(panel.grid.major = element_blank(), 
        panel.background = element_rect(fill = "aliceblue"))

bbox <- c(-30,0,30,60)
googlemap <- get_map(location = bbox)
ggmap(googlemap) +
  geom_polygon(data = mv95, aes(long, lat), alpha = 0.4) +
  geom_polygon(data = mv50, aes(long, lat), alpha = 0.8) +
  scale_fill_viridis(discrete=TRUE, option="A") +
  labs(x = " Longitude", y = "Latitude", title = 'Date: {gsub("1995-", "", current_frame)}')+
  transition_manual(frames = date)+
  theme_classic() +
  theme(legend.position="none", axis.text = element_text(color = "black"))

## KS library
kernels <- lapply(unique(date(locations_thin$datestamp))[1:3], function(f){
  print(f)
  data_sp <- locations_thin %>%  # store the data, then create a spatial object for plotting
    dplyr::filter(date(datestamp) == f) %>% 
    dplyr::select(location.long, location.lat)
  
  kern <- Hpi(x = coordinates(data_sp), binned = F)
  
  fhat <- kde(x = data_sp, H = get("kern"), compute.cont = T, binned = F)
  
  #png(filename= paste0("C:/Users/heste/Desktop/HB_storks/gif/", f, ".png"))
  #plot(fhat, xlab = "Longitude", ylab = "Latitude", xlim=c(-15,10), ylim = c(25, 50))
  #dev.off()
  
  return(fhat)
})

# rhr library (http://jmsigner.github.io/rhrman/cliHRE.html)
library(rhr)
kernels <- lapply(unique(date(locations_thin$datestamp)), function(f){
  print(f)
  
  data_sp <- locations_thin %>%  # store the data, then create a spatial object for plotting
    dplyr::filter(date(datestamp) == f) %>% 
    dplyr::select(location.long, location.lat) %>% 
    as.data.frame()
  
  #kd1 <- rhrKDE(data_sp)
  bw <- rhrHref(data_sp)
  #bw <- rhrHlscv(data_sp)
   
  kd <- rhrKDE(data_sp, h = bw$h)
  
  #png(filename= paste0("C:/Users/heste/Desktop/HB_storks/gif/", f, ".png"))
  #plot(fhat, xlab = "Longitude", ylab = "Latitude", xlim=c(-15,10), ylim = c(25, 50))
  #dev.off()
  
  return(kd)
})
iso <- rhrIsopleths(kernels[[10]], levels = seq(10, 90, 5))
plot(rhrUD(kernels[[153]]), xlim = c(-20, 20), ylim = c(10,60))
lines(iso, add = T)
maps::map("world", add = TRUE)

# make vector of the unique dates that each kernel belongs to
dates <- unique(date(locations_thin$datestamp))
# use the number in that to link to the item in the list

tracks <- tracks %>% 
  mutate(datestamp = timestamp) 

year(tracks$datestamp) <- 1995

t <- lapply(unique(date(tracks$datestamp)), function(x){

  df <- tracks %>% 
    filter(date(datestamp) == x)

  # for each unique date in the data, call the list item
  id <- which(dates == x)
  kd <- rhrUD(kernels[[125]])
  
  # extract the value for each location
  vals <- raster::extract(kd, df[, c("location.long", "location.lat")])
  
  # append
  df <- df %>% 
    mutate(kernel_dens = vals)
  
}) %>% reduce(rbind)

# Build 10 images -> save them at .png format
png(file="gif/%02d.png", width=480, height=480)

for (i in 1:length(dates)){
  maps::map("world", xlim = c(-40,20), ylim = c(10,60))
  plot(rhrUD(kernels[[i]]), add = T)
  maps::map("world", xlim = c(-40,20), ylim = c(10,60), add = T)
  mtext(paste(day(dates[i]), "-", month(dates[i])), side=3)
}
dev.off()

# Use image magick
system("*.png animated_count_down.gif")

# Remove png files
file.remove(list.files(pattern=".png"))

### save


dummydata <- data.frame(x = c(10, 6, 10), y = c("Wind support", "Uplift", "Uplift:Migration year"), sh = c(2, 2, 2), sl = c(-2, -2, -2))
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
        axis.line=element_line(size=8, color = "black"),
        axis.ticks.length=unit(8, "cm"),
        axis.ticks = element_blank(),
        axis.text.x.bottom = element_blank())

library(magick)
imgs <- list.files("C:/Users/heste/Desktop/HB_storks/gif/", full.names = TRUE)
img_list <- lapply(gtools::mixedsort(imgs), image_read)

## join the images together
img_joined <- image_join(img_list)

## animate at 2 frames per second
img_animated <- image_animate(img_joined, fps = 5)

image_write(image = img_animated,
            path = "C:/Users/heste/Desktop/HB_storks/kernels_rhr.gif")
