### Prepare rasters to annotate the observed stork movement and then 
### to create resistance layers
### 2025-05-21

library(sf)
library(lubridate)
library(terra)
library(dplyr)
library(ggplot2)
theme_set(theme_classic()+
            theme(axis.text = element_text(color = "black", size = 12), 
                  text = element_text(size = 15)))
colfunc <- colorRampPalette(c("#0081A7", "#0098b0", "#00AFB9", "#7fc4b8", "#FED9B7", "#f7a58f", "#f27e71", "#F07167", "#ee5e53"))

# 1. access elevation data (1 km)
# cnames <- maps::map("world", xlim = c(-30, 50), ylim = c(-30,55))$names
# countries <- geodata::country_codes() %>%
#   filter(NAME %in% cnames) %>%
#   dplyr::select(ISO3) %>%
#   deframe()
# add countries whose codes don't match, in addition Kosovo must be downloaded manually from https://geodata.ucdavis.edu/geodata/elv/ under the name "KO-_elv.zip"
# countries <- c(countries, "CIV", "GBR", "GNQ", "COG", "MWI", "PSE", "MKD", "SWZ")
# srtm <- lapply(countries, function(l){
#   elevation_30s(country = l, path = "/home/hbronnvik/Documents/chapter3/srtms/")
# })
srtm <- list.files("/home/hbronnvik/Documents/chapter3/srtms/elevation", full.names = T, pattern = ".tif")
srtm <- sprc(srtm)
srtm <- merge(srtm)
srtm <- crop(srtm, ext(-30, 50, -30, 60))
# get meters
srtm_m <- project(srtm, "ESRI:53009")
names(srtm_m) <- "DEM"
plot(srtm_m)
# writeRaster(srtm_m, file = "/home/hbronnvik/Documents/chapter3/env_rasts/moll_dem.tif", overwrite = T)
# srtm_m <- rast("/home/hbronnvik/Documents/chapter3/env_rasts/moll_dem.tif")

# this can be reduced to 90 m if needed:
# lons <- seq(from = -20, to = 15, by = 5)
# lats <- seq(5, 55, 5)
# locs <- expand.grid(lons, lats) %>% 
#   group_by(row_number()) %>% 
#   group_split()
# 
# srtm <- lapply(locs, function(l){
#   elevation_3s(l$Var1, l$Var2, path = "/home/hbronnvik/Documents/chapter3/srtms/")
# })
# srtm <- list.files("/home/hbronnvik/Documents/chapter3/srtms/elevation", full.names = T, pattern = ".tif")
# srtm <- sprc(srtm)
# srtm <- merge(srtm)
# plot(srtm)
# writeRaster(srtm, file = "/home/hbronnvik/Documents/chapter3/dem_3s.tif")
# srtm <- rast("/home/hbronnvik/Documents/chapter3/dem_3s.tif")

# 2. process DEM into products on original data (https://doi.org/10.1098/rsos.181440)
aspect <- terrain(srtm, "aspect")
slope <- terrain(srtm, "slope")
roughness <- terrain(srtm, "roughness")
tri <- terrain(srtm, "TRI")
tpi <- terrain(srtm, "TPI")

# reclassify aspects based on slope so that flat areas do not face any way
# if slope is less than or equal to 10 cm, no aspect
aspect[values(slope) <= 0.1] <- -1
plot(aspect)

# aspect binned to compass directions
# the bins
compass_points <- c("N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE",
                    "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW")
compass_breaks <- seq(0, 360, length.out = 17)
# define the bins as needed for classification
bins <- cbind(from = c(-1, compass_breaks[1:(length(compass_breaks)-1)]),
              to = c(-1, compass_breaks[-1]),
              becomes = c(-1, 1:16))
# reclassify the continuous aspects those using numeric codes
aspect_class <- classify(aspect, rcl = bins, include.lowest = T, right = F)
# also reclassify the NA values
aspect_class <- classify(aspect_class, cbind(NA, -1))
# convert to discrete compass points 
levels(aspect_class) <- data.frame(ID = c(-1, 1:16), compass = c("Flat", compass_points))
plot(aspect_class)

# a look at how these variables relate
c(srtm, aspect, slope, roughness, tri, tpi) %>% 
  as.data.frame(xy = T) %>% 
  tidyr::drop_na(TPI) %>%  
  tidyr::drop_na(TRI) %>% 
  tidyr::drop_na(roughness) %>% 
  tidyr::drop_na(slope) %>%
  tidyr::drop_na(DEM) %>% 
  cor() %>% 
  corrplot::corrplot.mixed(order = 'AOE')
# roughness and slope have a 0.98 correlation
# TRI and slope have 0.96 correlation
# slope and DEM have a 0.44 correlation

#                        x           y          DEM        aspect      slope   roughness
# x          1.0000000000 -0.166360758  0.283633469 -0.0622033818 0.08416346 0.085630492
# y         -0.1663607584  1.000000000 -0.313033341 -0.0197836092 0.09694945 0.102934627
# DEM        0.2836334693 -0.313033341  1.000000000  0.0026612383 0.43987567 0.460030930
# aspect    -0.0622033818 -0.019783609  0.002661238  1.0000000000 0.00579898 0.007279351
# slope      0.0841634635  0.096949453  0.439875668  0.0057989798 1.00000000 0.977239540
# roughness  0.0856304921  0.102934627  0.460030930  0.0072793510 0.97723954 1.000000000
# TRI        0.0825207649  0.102137205  0.457334817  0.0082525487 0.95996236 0.981044249
# TPI       -0.0003217311  0.001035777  0.087819358  0.0009445971 0.03812404 0.037519625
#                   TRI           TPI
# x         0.082520765 -0.0003217311
# y         0.102137205  0.0010357770
# DEM       0.457334817  0.0878193579
# aspect    0.008252549  0.0009445971
# slope     0.959962361  0.0381240363
# roughness 0.981044249  0.0375196254
# TRI       1.000000000  0.0738892401
# TPI       0.073889240  1.0000000000

# get meters projections
aspect <- project(aspect, srtm_m)
aspect_class <- project(aspect_class, srtm_m)
# writeRaster(aspect, file = "/home/hbronnvik/Documents/chapter3/env_rasts/aspect.tif")
# writeRaster(aspect_class, file = "/home/hbronnvik/Documents/chapter3/env_rasts/aspect_class.tif")

# 3. access high-res land cover data from https://esa-worldcover.org/en
# lapply(c("built", "grassland", "cropland", "wetland", "water"), function(v){
#   temp <- landcover(v, "/home/hbronnvik/Documents/chapter3")
#   temp <- crop(temp, ext(srtm))
#   terra::writeRaster(temp, paste0("/home/hbronnvik/Documents/chapter3/landuse/WorldCover_",
#                                  v, "_30s.tif"),
#                      overwrite = T)
# })
# plot(rast("/home/hbronnvik/Documents/chapter3/landuse/WorldCover_water_30s.tif"))

# 4. BLH as a static variable (365 days per year = 3791 layers)
blhs <- list.files("/home/hbronnvik/Documents/chapter3/PBLH", full.names = T)
# take the maximum per cell
blhs <- rast(blhs[!grepl("2013|2014", blhs)])
max_blhs <- quantile(blhs, 0.95)
names(max_blhs) <- "max_blh"
# take the median per cell
med_blhs <- median(blhs)
names(med_blhs) <- "med_blh"
# visualize
# png("/home/hbronnvik/Documents/chapter3/figures/med_95_blh.png",
#     width = 8.2, height = 11.7, units = "in", res = 400)
# library(rnaturalearth)
# world <- ne_countries(scale = "large", returnclass = "sf")
# lyr_names <- as_labeller(
#   c("max_blh" = "95th percentile PBL height",
#     "med_blh" = "Median PBL height"))
# ggplot() +
#   tidyterra::geom_spatraster(data = blhs)+
#   geom_sf(data=world, fill = NA, alpha = 0.5) +
#   coord_sf(xlim = c(-17, 50), ylim = c(-30, 60), expand = FALSE)+
#   tidyterra::scale_fill_whitebox_c(
#     palette = "muted",
#     n.breaks = 13,
#     guide = guide_legend(reverse = TRUE,
#                          title = "Boundary\nlayer height(m)")) +
#   labs(subtitle = "Boundary layer heights at 14:00 UTC 2015-2025") +
#   facet_wrap(~lyr, labeller = lyr_names)
# dev.off()
# get meters projections
med_blhs <- project(med_blhs, srtm_m)
max_blhs <- project(max_blhs, srtm_m)
# interpolate to smaller raster cell values, this won't change the outcome, just better match the elevation data
med_blhs <- terra::resample(med_blhs, srtm_m, method = "bilinear")
max_blhs <- terra::resample(max_blhs, srtm_m, method = "bilinear")
blhs <- c(med_blhs, max_blhs)
plot(blhs)
# blhs <- crop(blhs, srtm_m, extend = T)

# save out as a rough estimate of where boundary layer heights tend to be highest (uplift strongest)
# writeRaster(blhs, file = "/home/hbronnvik/Documents/chapter3/env_rasts/med_95_blh.tif")

max_blh_z <- scale(sqrt(max_blh))
plot(max_blh_z)

# 4. social information as a static variable
# access the full migration data
# 02_segment tracks & 09_eastern_segmentation
locsW <- readRDS("/home/hbronnvik/Documents/chapter3/migration_locations_40km70km15daySpeed_2025-03-08.rds")
locsE <- readRDS("/home/hbronnvik/Documents/chapter3/eastern_migration_locations_40km70km15daySpeed_2025-05-16.rds")

# the days classified as migratory, start and end migration dates, and outcome
metadW <- readRDS("/home/hbronnvik/Documents/chapter3/migration_dates_2025-03-08.rds") %>% 
  filter(journey_number == 1 & season == "fall")
metadE <- readRDS("/home/hbronnvik/Documents/chapter3/eastern_migration_dates_2025-05-16.rds") %>% 
  filter(journey_number == 1 & season == "fall")
metad <- rbind(metadW, metadE)

# recombine east and west individuals
locs <- locsW %>% 
  dplyr::select(trackID, individual.id, location.long, location.lat, timestamp) %>% 
  rbind(locsE %>% 
          dplyr::select(trackID, individual.id, location.long, location.lat, timestamp)) %>% 
  filter(trackID %in% unique(metad$trackID)) %>%
  # make an sf object to annotate
  st_as_sf(coords = c("location.long", "location.lat"), crs = crs(srtm)) %>%
  # and transform to the CRS of the background data
  st_transform(crs = crs(srtm_m)) %>%
  # add on the cell each observation fell in
  mutate(cell = terra::extract(srtm_m, ., cells = T)$cell,
         # and the coordinates
         location.long = st_coordinates(.)[,1],
         location.lat = st_coordinates(.)[,2])
# count the observations in each cell and vectorize
locs <- locs %>% 
  group_by(cell) %>% 
  # count total observations of different migrations in that cell
  mutate(total_obs = length(unique(trackID))) %>% 
  slice(1) %>% 
  ungroup() %>% 
  vect()
# rasterize the vector using DEM in moll proj as a template and total per-cell obs as value
locs_r <- rasterize(locs, srtm_m, field = "total_obs")
# if birds were not observed in a cell, place a 0
locs_r <- classify(locs_r, cbind(NA, 0))
names(locs_r) <- "obs"
plot(locs_r)
# save out as a rough estimate of where social information tends to be highest
# writeRaster(locs_r, file = "/home/hbronnvik/Documents/chapter3/env_rasts/total_obs_fall.tif")

# discretize the observations
# create bins with no (0-1], few (1-5], intermediate number (5-10], or many (> 10)
locs_b <- locs_r
locs_b <- classify(locs_b, cbind(c(-1,1,5,10), c(1,5,10,140), c(1,2,3,4)))
levels(locs_b) <- data.frame(ID = 1:4, Category = c("none", "few", "some", "many"))
plot(locs_b)
names(locs_b) <- "obs"
cons <- c(blhs[[2]], srtm_m, locs_b)

filled_3x <- lapply(cons, function(r){
  # fill the empty spaces with 0s 
  r <- classify(r, cbind(NA, 0))
  # then remove them except for sea crossing points
  r <- mask(r, crossing_coasts)
  return(r)
})
filled_3x <- rast(filled_3x)
plot(filled_3x)
filled_3x[["DEM"]] <- scale(log(filled_3x[["DEM"]]+420))
filled_3x[["max_blh"]] <- scale(filled_3x[["max_blh"]])
plot(filled_3x)

# library(MASS)
# x <- sample(values(filled_3x[[1]], na.rm = T), length(values(filled_3x[[1]], na.rm = T))/2)
# b <- boxcox(lm(x ~ 1))
# # Exact lambda
# lambda <- b$x[which.max(b$y)]
# lambda
# 0.8686869

# writeRaster(filled_3x, file = "/home/hbronnvik/Documents/chapter3/env_rasts/terr_z_3x.tif", overwrite = T)

# foraging areas more or less
# access high-res land cover data from https://esa-worldcover.org/en
lapply(c("grassland", "cropland", "wetland"), function(v){
  temp <- landcover(v, "/home/hbronnvik/Documents/chapter3")
  temp <- crop(temp, ext(srtm))
  temp <- project(temp, "ESRI:53009")
  terra::writeRaster(temp, paste0("/home/hbronnvik/Documents/chapter3/landuse/WorldCover_",
                                  v, "_30s.tif"),
                     overwrite = T)
})
crops_r <- rast("/home/hbronnvik/Documents/chapter3/landuse/WorldCover_cropland_30s.tif")
grass_r <- rast("/home/hbronnvik/Documents/chapter3/landuse/WorldCover_grassland_30s.tif")
wet_r <- rast("/home/hbronnvik/Documents/chapter3/landuse/WorldCover_wetland_30s.tif")
green_r <- crops_r+grass_r+wet_r
plot(green_r)
names(green_r) <- "forage"
# writeRaster(green_r, "/home/hbronnvik/Documents/chapter3/food_layer_mollweide.tif", overwrite = T)

# scale the cropland/grass/wetland layer 
forage <- rast("/home/hbronnvik/Documents/chapter3/food_layer_mollweide.tif")
# set oceans to 0 foraging
forage <- classify(forage, cbind(NA, 0))
forage <- mask(forage, crossing_coasts)
plot(forage)
# also pull in the binned aspects
filled_3x <- rast("/home/hbronnvik/Documents/chapter3/env_rasts/terr_z_3x.tif")
aspect_class <- rast("/home/hbronnvik/Documents/chapter3/env_rasts/aspect_class.tif")
aspect_class <- mask(aspect_class, crossing_coasts)
# add these together
filled4x <- c(filled_3x, forage, aspect_class)
filled4x[["obs"]] <- NULL
plot(filled4x)
# save
# writeRaster(filled4x, "/home/hbronnvik/Documents/chapter3/filled4x.tif", overwrite = T)

# this leads to 15 unique combinations excluding observations of conspecifics
items <- c("blh", "dem", "foo", "asp")
combs <- unlist(lapply(1:length(items), function(k) {
  combn(items, k, simplify = FALSE)
}), recursive = FALSE)

print(combs)



slope <- project(slope, crs(srtm_m))
roughness <- project(roughness, crs(srtm_m))
tri <- project(tri, crs(srtm_m))
tpi <- project(tpi, crs(srtm_m))
# combine all the variables into one giant raster and save it out for annotations
lyrs <- c(srtm_m, aspect, aspect_class, slope, roughness, tri, tpi, locs_r, blhs)
plot(lyrs)
# writeRaster(lyrs, file = "/home/hbronnvik/Documents/chapter3/env_rasts/full_terr.tif", overwrite = T)

### also prepare the prediction layers by allowing for water crossings where DEM is NA

# storks need to go where DEM is NA. 
# from Ras Mohammed to Hurghada is roughly 34 km, Spain-Morocco ~37, Dover-Calais ~34, Italy-Tunisia ~140
# this can act as an insulator in the Circuitscape models https://docs.circuitscape.org/Circuitscape.jl/latest/options/#Read-raster-mask-file:~:text=map%20(i.e.%2C-,treated%20as%20complete%20barriers,-).%20Positive%20integer%20cells
# but for memory reasons, we will not first predict across the oceans and then insulate, we will just remove the oceans first

# create polygons at places where storks cross
# wgs_sf <- sf::st_crs("+proj=longlat +datum=WGS84 +no_defs")

# xs <- c(35.55862,35.58746,36.00357,36.24939,36.23016,35.78384)
# ys <- c(36.57032,36.76859,36.95648,36.79829,36.59789,36.31844)
# iskenderun <- data.frame(x = xs, y = ys, id = "iskenderun")

# suez <- data.frame(x = xs, y = ys, id = "suez")
# bosphorus <- data.frame(x = xs, y = ys, id = "bosphorus")
# dover <- data.frame(x = xs, y = ys, id = "dover")
# gibraltar <- data.frame(x = xs, y = ys, id = "gibraltar")

# poly_xy <- list(gibraltar, bosphorus, suez, iskenderun, dover)
# poly_xy <- poly_xy %>% 
#   reduce(rbind)
# saveRDS(poly_xy, file = "/home/hbronnvik/Documents/chapter3/sea_crossing_polygons.rds")
poly_xy <- readRDS("/home/hbronnvik/Documents/chapter3/sea_crossing_polygons.rds")
localpoly <- poly_xy %>% 
  group_by(id) %>% 
  group_split() %>% 
  lapply(function(x) rbind(x,x[1,])) %>%
  lapply(function(x) x[,1:2]) %>%
  lapply(function(x) list(as.matrix(x))) %>%
  lapply(function(x) sf::st_polygon(x)) 
#convert polygons to sf object and add id column
crossing_points <- localpoly %>% 
  sf::st_sfc(crs = "EPSG:4326") %>% 
  sf::st_sf(geom=.) %>% 
  mutate(individual.id=names(localpoly)) %>% 
  sf::st_transform(crs = "ESRI:53009")
# visualize
mapview::mapview(crossing_points, color = "black", alpha.regions = 0, lwd = 2)
# add to coastlines
build <- rasterize(crossing_points, srtm_m)
plot(build)
coasts <- srtm_m
coasts[!is.na(values(coasts))] <- 1
crossing_coasts <- merge(coasts, build)
plot(crossing_coasts)
# writeRaster(crossing_coasts, filename = "/home/hbronnvik/Documents/chapter3/barrier_mask.tif")

# now we have a coastline buffer allowing for crossings where storks actually cross
# essentially allowing us to include water bodies in the model as something they always avoid
# we can reclassify the terrain variables so that oceans are 0 and then buffer
crossing_coasts <- rast("/home/hbronnvik/Documents/chapter3/barrier_mask.tif")
filled_terr <- rast("/home/hbronnvik/Documents/chapter3/env_rasts/full_terr.tif")
filled_terr <- lapply(names(filled_terr)[!grepl("TRI|TPI|roughness|aspect", names(filled_terr))], function(x){
  print(paste0("Filling in the seas for ", x, "."), quote = F)
  if(x != "aspect|compass"){
    r <- filled_terr[[x]]
    r <- classify(r, cbind(NA, 0))
  }else{
    r <- filled_terr[[x]]
  }
  r <- mask(r, crossing_coasts)
  return(r)
}) # pared down number of layers due to crashing
gc()
filled_terr <- rast(filled_terr)
levels(filled_terr[["compass"]]) <- data.frame(ID = c(-1, 1:16), compass = c("Flat", compass_points))
plot(filled_terr)
# write out another huge raster with all the sea crossings
# writeRaster(filled_terr, file = "/home/hbronnvik/Documents/chapter3/env_rasts/full_terr_crossings.tif", overwrite = T)
filled_terr <- rast("/home/hbronnvik/Documents/chapter3/env_rasts/full_terr_crossings.tif")

# normalized and scaled
srtm_m_z <- scale(log(srtm_m+420))
hist(values(srtm_m_z))
plot(srtm_m_z)

# and a huge data frame with all the variables at every cell for model predicting
filled_df <- as.data.frame(filled_terr[[c("DEM", "compass", "obs", "med_blh", "max_blh")]], 
                           xy = T, na.rm = FALSE)
# scaled to match what we do when modeling
filled_df <- filled_df %>% 
  mutate(sqrt_dem_z = scale(log(DEM+424))[,1], 
         obs_z = scale(obs)[,1], 
         med_blh_z = scale(med_blh)[,1], 
         max_blh_z = scale(max_blh)[,1])
# saveRDS(filled_df, file = "/home/hbronnvik/Documents/chapter3/full_terr_df.rds")

# finally, split up the rasters into smaller pieces
# model predictions require a LOT of memory (> 2125.798687 GB using this raster & models),
# split up the data to make smaller chunks for R
dim(filled_terr) # rows, cols, layers, 11867  7654     6
# set splitting parameters
rc <- ceiling(dim(filled_terr)[1:2] / c(2,2))
# make four quarters
# x <- makeTiles(filled_terr, rc, filename = "/home/hbronnvik/Documents/chapter3/Raven/pred_tiles_static/static_var_tile_.tif", overwrite = T)


