### Include a comparison of both flyways for both species
### Hester Bronnvik
### 2025-08-19
### hbronnvik@ab.mpg.de

library(terra)
library(sf)
library(tidyverse)
library(cowplot)
theme_set(theme_classic()+
            theme(axis.text = element_text(color = "black", size = 12), 
                  text = element_text(size = 15),
                  strip.background = element_rect(colour="white", fill="white")))
colfunc <- colorRampPalette(c("#0081A7", "#0098b0", "#00AFB9", "#7fc4b8", "#FED9B7", "#f7a58f", "#f27e71", "#F07167", "#ee5e53"))

# the white stork data annotated with current as a list of two flyways (in the fall)
ws <- readRDS("/home/hbronnvik/Documents/chapter3/fall_a_data_current_cent_uds4x.rds")
ws <- lapply(ws, function(x){
  x <- x %>% 
    # simplify and match to the black storks
    st_drop_geometry() %>% 
    mutate(flyway = ifelse(flyway == "BMF", "BSM", flyway))
  # get the relevant columns
  x <- x[, -grep("uds", names(x))]
  timen <- which(names(x) == "timestamp")
  coln <- grep("percentile_res_blh_dem_foo|current_res_blh_dem_foo|location.|flyway|trackID|individual.id", names(x))
  x <- x[, c(timen, coln)]
  # rename them
  names(x) <- sub("EAF", "", 
                  sub("_res_blh_dem_foo_", "", 
                      sub("BMF", "", 
                          sub("trackID", "track_id",
                              names(x)))))
  x$species <- "white"
  return(x)
})
# the black stork data annotated with current as a list of two flyways (in the fall)
bs <- readRDS("/home/hbronnvik/Documents/chapter3/black_stork_demXblhXfood_current.rds")
bs <- lapply(bs, function(x){
  x <- x %>% 
    st_drop_geometry()
  timen <- which(names(x) == "timestamp")
  coln <- grep("percentile|current|location.|flyway|individual.id|track_id", names(x))
  x <- x[, c(timen, coln)]
  names(x) <- sub("E", "", names(x))
  x$species <- "black"
  return(x)
})

# do two comparisons: a model predicting percentile of the current using model ID,
# and simply finding the average percentile produced by each model
storks <- c(ws, bs) %>% 
  reduce(rbind) %>% 
  mutate(ID = as.factor(individual.id),
         trackID = as.factor(track_id),
         model = paste(species, flyway, sep = "_"),
         model = factor(model, levels = c("white_EAF", "white_BSM", "black_EAF", "black_BSM"))) %>% 
  # add on relocation to have a record of movements in line with current
  group_by(ID, trackID) %>% 
  arrange(timestamp) %>% 
  mutate(reloc = 1:n(),
         reloc = as.factor(reloc),
         len = n()) %>% 
  ungroup() %>% 
  filter(len >= 7) %>% 
  dplyr::select(ID, trackID, reloc, percentile, model, species)

model <- lme4::lmer(percentile~model+(1|ID/trackID), data = storks)
summary(model)

graph <- confint(model, method = "Wald") %>%
  as.data.frame() %>%
  rownames_to_column(var = "Model") %>%
  filter(!grepl("sig|cept", Model)) %>%
  left_join(summary(model)$coefficients %>%
              as.data.frame() %>%
              rownames_to_column(var = "Model")) %>%
  rename(lower = "2.5 %",
         upper = "97.5 %") %>% 
  mutate(model = sub("model", "", Model))

# reshape the data
p_graph <- storks %>% 
  group_by(model) %>% 
  summarize(Mean = round(mean(percentile, na.rm = T), 2),
            Std.Dev. = round(sd(percentile, na.rm = T), 2),
            species = ifelse(unique(species) == "white", "White stork", "Black stork")) %>% 
  left_join(graph) %>% 
  # add the reference level to the data
  mutate(Estimate = ifelse(model == "white_EAF", 0, Estimate)) %>% 
  # add the top row of the "table"
  rbind(data.frame(x = names(.), y = NA) %>%
          pivot_wider(names_from = x, values_from = y)) %>% 
  dplyr::select(model:Estimate) %>%
  # reorder, rename, and round as needed
  mutate(flyway = ifelse(grepl("BSM", Model), "Black Sea Mediterranean",
                         "East Atlantic"),
         group = ifelse(!is.na(species), paste0(species, "<br>", flyway), NA),
         Mean = sprintf("%.1f", Mean),
         Mean = ifelse(Mean == "NA", NA, Mean),
         Std.Dev. = sprintf("%.1f", Std.Dev.),
         Std.Dev. = ifelse(Std.Dev. == "NA", NA, Std.Dev.),
         group = factor(group, levels = c("Black stork<br>Black Sea Mediterranean",
                                          "Black stork<br>East Atlantic",
                                          "White stork<br>Black Sea Mediterranean",
                                          "White stork<br>East Atlantic",
                                          NA))) %>% 
  # align the Means and Std.Dev.s to the Estimates
  arrange(group)

# png("/home/hbronnvik/Documents/chapter3/figures/curr_comp_lmer_4x.png",
#     height = 8.2/2, width = 11.7, units = "in", res = 400)
ggplot(p_graph, aes(Estimate, group, color = flyway)) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_vline(xintercept = c(2.75, 4.85)) +
  geom_pointrange(aes(xmin = lower, xmax = upper), lwd = 1.5, size = 0.5) +
  lims(x = c(-15, 6.5)) + 
  annotate("text", x = c(3, 5), y = nrow(p_graph), label = c("Mean", "Std. Dev."), 
           fontface = "bold", hjust = 0, size = 5) + 
  annotate("text", x = 3, y = 1:nrow(p_graph), label = p_graph$Mean, hjust = 0, size = 5) + 
  annotate("text", x = 5, y = 1:nrow(p_graph), label = p_graph$Std.Dev., hjust = 0, size = 5) +
  scale_color_manual(values = c("#0098b0","#f27e71")) +
  labs(y = "", x = "Estimate") +
  theme(legend.position = "none",
        axis.text.y = ggtext::element_markdown(size = 15))
# dev.off()

x <- x %>% 
  arrange(Mean) %>%
  mutate(Mean = sprintf("%.1f", Mean),
         Std.Dev. = sprintf("%.1f", Std.Dev.)) %>% 
  rbind(x2) %>% 
  mutate(Model = factor(Model, levels = Model))

x %>% 
  ggplot(aes(Estimate, Model)) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_vline(xintercept = c(2.75, 3.85)) +
  geom_pointrange(aes(xmin = lower, xmax = upper, color = ), lwd = 1.5, size = 0.5) +
  scale_x_continuous(limits = c(-1, 5), breaks = c(0, 1, 2, 3)) +
  labs(y = "", subtitle = unique(x$flyway)) + 
  annotate("text", x = c(3, 4), y = nrow(x), label = c("Mean", "Std. Dev."), 
           fontface = "bold", hjust = 0, size = 5) + 
  annotate("text", x = 3, y = 1:nrow(x), label = x$Mean, hjust = 0, size = 5) + 
  annotate("text", x = 4, y = 1:nrow(x), label = x$Std.Dev., hjust = 0, size = 5) +
  theme(axis.text.y = ggtext::element_markdown(size = 15))


library(rnaturalearth)
world <- ne_countries(scale = "large", returnclass = "sf") %>% 
  st_transform(crs = "ESRI:53009")
# a plot of the tracks
tracks <- c(ws, bs) %>% 
  reduce(rbind)%>% 
  mutate(flyway = ifelse(grepl("BSM", flyway), "Black Sea\nMediterranean",
                         "East Atlantic"),
         species = ifelse(species == "white", "A) White stork", "B) Black stork"),
         species = factor(species, levels = c("A) White stork", "B) Black stork")))

png(filename = "/home/hbronnvik/Documents/chapter3/figures/tracks_4x.png",
    height = 8, width = 11, units = "in", res = 500)
ggplot() + 
  geom_sf(data = world, color = "gray80", fill = "gray80") +
  geom_path(data = tracks, aes(location.long, location.lat, group = track_id, color = flyway), 
            alpha = 0.5) +
  geom_vline(xintercept = 915782.8, lty = 2) +
  coord_sf(xlim = c(-3003318, 5005879.4100524),
           ylim = c(-3639482.55008231, 6869064.04589582)) +
  scale_x_continuous(expand = c(0, 0), labels = function(x) as.character(x), 
                     breaks = seq(-30, 50, 20)) +
  scale_y_continuous(expand = c(0, 0), labels = function(x) as.character(x)) +
  scale_color_manual(values = c("#0098b0","#f27e71")) +
  labs(x = "Longitude", y = "Latitude", color = "Flyway") +
  facet_wrap(~species) +
  theme(strip.text = element_text(hjust = 0))
dev.off()
#

# explore properties of the flyways in the tracks
###------------------------------

## get the actual locations used by white storks
# the EAF days classified as migratory, start and end migration dates, and outcome
w_metad <- readRDS("/home/hbronnvik/Documents/chapter3/migration_dates_2025-03-08.rds") %>% 
  filter(season == "fall")
# the EAF location data during migration
west_ml <- readRDS("/home/hbronnvik/Documents/chapter3/migration_locations_40km70km15daySpeed_2025-03-08.rds")
# only migratory flight locations
# filter ground speeds and daily distances (speed is scalar and a few missing data can affect it dramatically)
west_ml <- west_ml %>% 
  mutate(id_date = paste(individual.id, date, sep = "_")) %>% 
  filter(id_date %in% w_metad$id_date & ground_speed_15 > 2 & ground_speed_15 <= 21) %>% 
  st_as_sf(coords = c("location.long", "location.lat"), crs = 4326) %>% 
  st_transform("ESRI:53009")
# the BMF days classified as migratory, start and end migration dates, and outcome
e_metad <- readRDS("~/Documents/chapter3/eastern_migration_dates_2025-05-16.rds") %>% 
  filter(season == "fall")
# the BMF location data during migration
east_ml <- readRDS("~/Documents/chapter3/eastern_migration_locations_40km70km15daySpeed_2025-05-16.rds") %>% 
  mutate(id_date = paste(individual.id, date, sep = "_")) %>% 
  filter(id_date %in% e_metad$id_date & ground_speed_15 > 2 & ground_speed_15 <= 21) %>% 
  st_as_sf(coords = c("location.long", "location.lat"), crs = 4326) %>% 
  # st_crop(xmin = 36.44682347007568, xmax = 31.004823415381654, 
  #         ymin = 39.08019820411506, ymax = 31.80700524124481)
  st_transform("ESRI:53009")

# environmental raster
srtm_m <- rast("/home/hbronnvik/Documents/chapter3/env_rasts/moll_dem.tif")
blh_m <- rast("/home/hbronnvik/Documents/chapter3/env_rasts/med_95_blh.tif")[["max_blh"]]
forage <- rast("/home/hbronnvik/Documents/chapter3/food_layer_mollweide.tif")
lyrs <- c(srtm_m, blh_m, forage)

# build <- as.data.frame(lyrs)
# cor.test(build$DEM, build$max_blh, method = "spearman")
# data:  build$DEM and build$max_blh
# S = 1.5442e+22, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# 0.3042863 

# white stork data
both_df <- west_ml %>%
  dplyr::select(geometry, timestamp) %>%
  mutate(flyway = "EAF") %>%
  rbind(east_ml %>%
          dplyr::select(geometry, timestamp) %>%
          mutate(flyway = "BMF")) %>% 
  mutate(terra::extract(lyrs, ., method = "bilinear"),
         flyway = ifelse(grepl("BMF", flyway), "Black Sea Mediterranean",
                         "East Atlantic"),
         species = "White stork")
# the black stork data annotated with current as a list of two flyways (in the fall)
bs <- readRDS("/home/hbronnvik/Documents/chapter3/black_stork_demXblhXfood_current.rds")
bs <- lapply(bs, function(x){
  x <- x %>% 
    dplyr::select(geometry, timestamp, flyway) %>% 
    mutate(terra::extract(lyrs, ., method = "bilinear"),
           flyway = ifelse(grepl("BSM", flyway), "Black Sea Mediterranean",
                           "East Atlantic"),
           species = "Black stork")
  return(x)
}) %>% reduce(rbind)

both_df <- both_df %>% 
  rbind(bs) %>% 
  st_drop_geometry() %>% 
  mutate(group = paste0(species, "\n", flyway)) %>% 
  dplyr::select(group, DEM, max_blh, forage, flyway)

both_df %>% 
  drop_na(DEM) %>% 
  drop_na(forage) %>% 
  dplyr::select(DEM, max_blh, forage) %>% 
  cor()

both_df <- both_df %>% 
  pivot_longer(cols = c(DEM, max_blh, forage), names_to = "variable")

ps <- lapply(split(both_df, both_df$variable), function(x){
  l <- ifelse(unique(x$variable) == "DEM", "Elevation (m)", 
              ifelse(unique(x$variable) == "max_blh", "Boundary layer height (m)",
                     "Percent grassland/cropland/wetland"))
  x <- x %>% 
    arrange(flyway)
  ggplot(x, aes(value, group, fill = flyway)) +
    ggridges::stat_density_ridges(quantile_lines = T, quantiles = 2,
                                  panel_scaling = T, alpha = 0.85) +
    scale_fill_manual("Flyway", values = c("#0098b0", "#f27e71")) +
    labs(x = l, y = "") +
    theme(legend.position = "none")
})
# png(filename = "/home/hbronnvik/Documents/chapter3/figures/env_dists.png",
#     height = 8, width = 11, units = "in", res = 500)
ggpubr::ggarrange(ps[[1]], ps[[2]], ps[[3]], nrow = 3)
# dev.off()
#


# the layers are too big to do all at once, R will crash (31.0 GiB, 12th Gen Intel® Core™ i7-12700 × 20)
# full <- lyrs %>% 
#   as.data.frame(xy = T) %>% 
#   mutate(flyway = ifelse(x > 915782.8, "Black Sea Mediterranean", "East Atlantic")) %>% 
#   pivot_longer(cols = 1:3) %>% 
#   mutate(variable = ifelse(name == "DEM", "Elevation (m)", 
#                            ifelse(name == "max_blh", "Boundary layer height (m)",
#                                   "Percent grassland/cropland/wetland")))

# compare the two halves of the background layers
ext_EAF <- ext(c(xmin(lyrs), 915782.8, -3074066, 6589402))
ext_BSM <- ext(c(915782.8, 4675682, -3074066, 6589402))

EAF <- values(crop(lyrs, ext_EAF), na.rm = T) %>% 
  as.data.frame() %>% 
  mutate(flyway = "East Atlantic") %>% 
  pivot_longer(cols = 1:3, names_to = "variable")

BSM <- values(crop(lyrs, ext_BSM), na.rm = T) %>% 
  as.data.frame() %>% 
  mutate(flyway = "Black Sea Mediterranean") %>% 
  pivot_longer(cols = 1:3, names_to = "variable")

full <- rbind(EAF, BSM)

background_ps <- lapply(split(full, full$variable), function(x){
  l <- ifelse(unique(x$variable) == "DEM", "Elevation (m)", 
              ifelse(unique(x$variable) == "max_blh", "Boundary layer height (m)",
                     "Percent grassland/cropland/wetland"))
  x <- x %>% 
    slice_sample(prop = 0.2)
  ggplot(x, aes(value, flyway, fill = flyway)) +
    ggridges::stat_density_ridges(quantile_lines = T, quantiles = 2,
                                  panel_scaling = T, alpha = 0.85) +
    scale_fill_manual("Flyway", values = c("#D575AC", "#8FC87B")) +
    labs(x = l, y = "") +
    theme(legend.position = "none")
})
# png(filename = "/home/hbronnvik/Documents/chapter3/figures/background_env_dists.png",
#     height = 8, width = 11, units = "in", res = 500)
ggpubr::ggarrange(background_ps[[3]], background_ps[[1]], background_ps[[2]], nrow = 3)
# dev.off()
#

# make polygons around the two flyways
fly_polys <- c(ws, bs) %>% 
  reduce(rbind) %>% 
  # exclude the 6 birds that took apparently unique routes to get a more general polygon
  filter(!individual.id %in% c(individual.id[which.min(location.lat)], 11089898,
                               292304721, 2908772879, 1585123983, 930267063)) %>% 
  st_as_sf(coords = c("location.long", "location.lat"), crs = crs(srtm_m)) %>% 
  group_by(flyway) %>% 
  summarise() %>% 
  st_cast("POLYGON") %>% 
  geos::geos_concave_hull(allow_holes = F, ratio = 0.2) %>% 
  st_as_sf() %>% 
  st_set_crs(crs(lyrs)) %>% 
  mutate(flyway = c("Black Sea Mediterranean", "East Atlantic"))

library(rnaturalearth)
world <- ne_countries(scale = "large", returnclass = "sf") %>% 
  st_transform(crs = crs(lyrs))
crossing_coasts <- rast("/home/hbronnvik/Documents/chapter3/barrier_mask.tif")

fly_polys <- lapply(1:nrow(fly_polys), function(x){
  poly <- fly_polys[x,] %>% 
    vect()
  
  poly <- rasterize(poly, crossing_coasts, field = "flyway")
  poly <- mask(poly, crossing_coasts)
  as.polygons(poly, dissolve = TRUE)
}) %>% vect()

# polys <- ggplot() +
#   geom_sf(data = world, color = "gray80", fill = "gray80") +
#   tidyterra::geom_spatvector(data = fly_polys, aes(fill = flyway, color = flyway), 
#           alpha = 0.75, lwd = 0.5) +
#   coord_sf(xlim = c(-2112348, 4534867), ylim = c(-3639776, 6869064)) +
#   scale_fill_manual("Flyway", values = c("#0098b0", "#f27e71")) +
#   scale_color_manual("Flyway", values = c("#0081A7", "#ee5e53")) +
#   labs(x = "Longitude", y = "Latitude") +
#   theme(legend.position = "top")

# see flyway_env_plots.R for background plotting of the rasters

# get the full backgrounds
background_polys <- lapply(1:nlyr(lyrs), function(l){
  lyr <- lyrs[[l]]
  varn <- ifelse(names(lyr) == "DEM", "Elevation (m)", 
                 ifelse(names(lyr) == "max_blh", "Boundary layer height (m)",
                        "Percent grassland/cropland/wetland"))
  
  # fill the empty spaces with 0s 
  lyr <- classify(lyr, cbind(NA, 0))
  # then remove them except for sea crossing points
  lyr <- mask(lyr, crossing_coasts)
  
  # go through each polygon and pull out the values
  fly_vals <- lapply(1:nrow(fly_polys), function(f){
    x <- mask(lyr, fly_polys[f]) %>% 
      as.data.frame() %>% 
      mutate(flyway = fly_polys[f]$flyway) %>% 
      rename(value = 1)
    return(x)
  })
  lyr_vals <- data.table::rbindlist(fly_vals)
  
  lns <- lyr_vals %>% 
    group_by(flyway) %>% 
    summarize(mn = min(value),
              mx = max(value))
  
  ggplot(lyr_vals, aes(value, flyway, fill = flyway)) +
    ggridges::stat_density_ridges(quantile_lines = T, quantiles = 2,
                                  panel_scaling = T, alpha = 0.85) +
    geom_vline(xintercept = c(lns$mn[1], lns$mx[1]), color = "#0098b0", lty = 2) +
    geom_vline(xintercept = c(lns$mn[2], lns$mx[2]), color = "#f27e71", lty = 2) +
    scale_fill_manual("Flyway", values = c("#0098b0", "#f27e71")) +
    labs(x = varn, y = "") +
    theme(legend.position = "none",
          axis.ticks.y = element_blank())
})


bg <- ggpubr::ggarrange(background_polys[[3]], background_polys[[1]], background_polys[[2]], 
                  nrow = 3, labels = c("B", "C", "D"))
# png(filename = "/home/hbronnvik/Documents/chapter3/figures/poly_background_env_dists.png",
#     height = 8, width = 11, units = "in", res = 500)
ggpubr::ggarrange(polys, bg, labels = c("A", NA))
# dev.off()

flyway_rasts <- lapply(1:nlyr(lyrs), function(i){
  varn <- ifelse(names(lyrs[[i]]) == "DEM", "Elevation (m)", 
                 ifelse(names(lyrs[[i]]) == "max_blh", "Boundary layer height (m)",
                        "Percent grassland/cropland/wetland"))
  r <- project(lyrs[[i]], "EPSG:4326")
  ggplot() +
    tidyterra::geom_spatraster(data = r) +
    tidyterra::scale_fill_grass_c(palette = "grey") +
    geom_sf(data = fly_polys, aes(color = flyway), fill = NA, lwd = 0.5) +
    scale_x_continuous(expand = c(0, 0), limits = c(-17, 50)) +
    scale_y_continuous(expand = c(0, 0), limits = c(-30, 60)) +
    labs(fill = varn, color = "Flyway") +
    theme(legend.position = "top",
          legend.box = "vertical",
          legend.text = element_text(size = 10),
          legend.title=element_text(size=12),
          axis.text = element_blank(),
          axis.ticks = element_blank()) +
    guides(fill = guide_colorbar(order = 2, title.position = "top", title.hjust = 0.5,
                                 barheight = unit(0.5, "cm"), barwidth = unit(6, "cm")),
           color = guide_legend(order = 1, title.position = "top", title.hjust = 0.5,
                                barheight = unit(0.4, "cm"), barwidth = unit(.4, "cm")))
})
# png(filename = "/home/hbronnvik/Documents/chapter3/figures/poly_background_rasts.png",
#     height = 8, width = 11, units = "in", res = 500)
ggpubr::ggarrange(flyway_rasts[[1]], flyway_rasts[[2]], 
                  flyway_rasts[[3]], labels = "AUTO", nrow = 1)
# dev.off()
png(filename = "/home/hbronnvik/Documents/chapter3/figures/poly_background_rast_dists.png",
    height = 8, width = 11, units = "in", res = 500)
ggpubr::ggarrange(ggpubr::ggarrange(flyway_rasts[[1]], flyway_rasts[[2]], 
                            flyway_rasts[[3]], labels = "AUTO", nrow = 1),
          ggpubr::ggarrange(background_polys[[1]]+theme(axis.text.y = element_blank()), 
                            background_polys[[2]]+theme(axis.text.y = element_blank()), 
                            background_polys[[3]]+theme(axis.text.y = element_blank()), 
                            nrow = 1, labels = c("D", "E", "F")), 
          nrow = 2, heights = c(3,2))
dev.off()

# check correlations
cors <- lapply(1:nrow(fly_polys), function(x){
  mask(lyrs, fly_polys[x]) %>% 
    as.data.frame() %>% 
    drop_na() %>% 
    cor(method = "spearman") %>% 
    as.data.frame() %>% 
    mutate(flyway = fly_polys[x]$flyway)
}) %>% reduce(rbind)


# a plot of the current
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


# the eastern white storks burst using 03_step_generation
av_df <- readRDS("/home/hbronnvik/Documents/chapter3/eastern_used_av_df_60min.rds") %>% 
  st_as_sf(coords = c("location.long", "location.lat"), crs = 4326) %>% 
  st_transform(crs(lyrs)) %>% 
  mutate(terra::extract(lyrs, .)) %>% 
  st_drop_geometry() %>% 
  group_by(track, stratum) %>% 
  summarize(step_dem_var = stats::sd(DEM, na.rm = T),
            step_blh_var = stats::sd(max_blh, na.rm = T),
            step_foo_var = stats::sd(forage, na.rm = T),
            # min_dem = min(DEM, na.rm = T),
            # max_dem = max(DEM, na.rm = T),
            flyway = "BSM") %>% 
  ungroup() %>% rbind(readRDS("/home/hbronnvik/Documents/chapter3/used_av_df_60min.rds") %>% 
                        st_as_sf(coords = c("location.long", "location.lat"), crs = 4326) %>% 
                        st_transform(crs(lyrs)) %>% 
                        mutate(terra::extract(lyrs, .)) %>% 
                        st_drop_geometry() %>% 
                        group_by(track, stratum) %>% 
                        summarize(step_dem_var = stats::sd(DEM, na.rm = T),
                                  step_blh_var = stats::sd(max_blh, na.rm = T),
                                  step_foo_var = stats::sd(forage, na.rm = T),
                                  # min_dem = min(DEM, na.rm = T),
                                  # max_dem = max(DEM, na.rm = T),
                                  flyway = "EAF") %>% 
                        ungroup())

long_av_df <- av_df %>% 
  pivot_longer(cols = contains("step"), names_to = "variable") %>% 
  mutate(variable = ifelse(variable == "step_foo_var", "Land use for foraging (%)",
                           ifelse(variable == "step_blh_var", "Boundary layer heights (m)",
                                  "Elevation (m)")))

choices <- ggplot(long_av_df, aes(value, fill = flyway)) +
  geom_density(alpha = 0.85) +
  scale_fill_manual(values = c("#0098b0", "#f27e71")) +
  labs(x = "Per-hour standard deviation in conditions", y = "Density", fill = "Flyway") +
  facet_wrap(~variable, scales = "free", nrow = 3)

var_max <- ggplot(av_df, aes(max_dem, step_dem_var)) +
  geom_point(alpha = 0.5) +
  geom_smooth(aes(color = flyway)) +
  scale_color_manual(values = c("#0098b0", "#f27e71")) +
  labs(y = "Per-step sd in elevation (m)", x = "Per-step maximum elevation (m)", color = "Flyway") 

cors <- lapply(split(av_df, av_df$flyway), function(x){
  flyway <- unique(x$flyway)
  if(flyway == "BSM"){
    cols <- c("#0081A7", "#009CB1", "#98C8B7")
  }else{
    cols <- c("#EE655A", "#F07369", "#F69D89")
  }
  x %>% 
    dplyr::select(max_dem, min_dem, step_dem_var) %>% 
    drop_na(max_dem) %>% 
    drop_na(min_dem) %>% 
    drop_na(step_dem_var) %>% 
    rename("Max. DEM" = max_dem,
           "Min. DEM" = min_dem,
           "Std. dev." = step_dem_var) %>% 
    cor() %>% 
    ggcorrplot::ggcorrplot(hc.order = TRUE,
                           type = "lower",
                           lab = TRUE,
                           colors = cols,
                           tl.srt = 90) +
    labs(title = paste0(flyway, " correlations in per-step values")) +
    theme(panel.grid = element_blank(),
          axis.text = element_text(color = "black"))
})

# png(filename = "/home/hbronnvik/Documents/chapter3/figures/flyway_dem_variance.png",
#     height = 8, width = 11, units = "in", res = 500)
ggpubr::ggarrange(choices, var_max, 
                  cors[[1]], cors[[2]],
                  nrow = 2, ncol = 2)
# dev.off()

# png(filename = "/home/hbronnvik/Documents/chapter3/figures/flyway_variance.png",
#     height = 8, width = 11, units = "in", res = 500)
choices
# dev.off()

