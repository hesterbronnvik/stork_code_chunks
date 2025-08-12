### Comparing the outputs of Circuitscape using different resistance surfaces
### Hester Bronnvik
### hbronnvik@ab.mpg.de
### 2025-06-03

library(sf)
library(terra)
library(dplyr)
library(purrr)
library(ggplot2)
library(ggpmisc)
theme_set(theme_classic()+
            theme(axis.text = element_text(color = "black", size = 12), 
                  text = element_text(size = 15),
                  strip.background = element_rect(colour="white", fill="white")))
colfunc <- colorRampPalette(c("#0081A7", "#0098b0", "#00AFB9", "#7fc4b8", "#FED9B7", "#f7a58f", "#f27e71", "#F07167", "#ee5e53"))

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

library(rnaturalearth)
ortho_proj <- "+proj=ortho +lat_0=15 +lon_0=10"
world <- ne_coastline(scale = "large", returnclass = "sf") %>% 
  st_transform(ortho_proj)
# png(filename = "/home/hbronnvik/Documents/chapter3/figures/flyway_tracks.png",
#     height = 11, width = 11, units = "in", res = 500)
ggplot() +
  geom_sf(data = world, fill = "gray80", color = "gray80") +
  # ggforce::geom_circle(aes(x0 = 10, y0 = 15, r = 64e5), color = "gray80") +
  geom_sf(data = west_ml %>%
               dplyr::select(geometry) %>%
               mutate(flyway = "EAF") %>%
               rbind(east_ml %>%
                       dplyr::select(geometry) %>%
                       mutate(flyway = "BMF")),
          aes(color = flyway), size = 0.5, alpha = 0.5) +
  # coord_sf(xlim = c(-1831566, 5039191), ylim = c(-3939776, 6169064), expand = T) +
  scale_color_manual("Flyway", values = c("#0098b0", "#f27e71"))
# dev.off()
# background map
# srtm_m <- rast("/home/hbronnvik/Documents/chapter3/env_rasts/moll_dem.tif")
# ggplot() +
#   tidyterra::geom_spatraster(data = srtm_m) +
#   tidyterra::scale_fill_whitebox_c(palette = "high_relief", guide = "none") +
#   geom_sf(data = west_ml %>%
#             dplyr::select(geometry) %>%
#             mutate(flyway = "EAF") %>%
#             rbind(east_ml %>%
#                     dplyr::select(geometry) %>%
#                     mutate(flyway = "BMF")),
#           aes(color = flyway), size = 0.5) +
#   scale_x_continuous(expand = c(0,0)) +
#   scale_y_continuous(expand = c(0,0)) +
#   scale_color_manual("Flyway", values = c("#f27e71", "#0098b0"))

both_df <- west_ml %>%
  dplyr::select(geometry, timestamp, trackID) %>%
  mutate(flyway = "EAF") %>%
  rbind(east_ml %>%
          dplyr::select(geometry, timestamp,trackID) %>%
          mutate(flyway = "BMF")) %>%
  mutate(seq_60 = round_date(timestamp, "hour")) %>% 
  group_by(trackID, seq_60) %>% 
  slice(1) %>% 
  ungroup() %>% 
  dplyr::select(geometry, flyway) %>% 
  st_transform("EPSG:4326") %>%
  # add on the current values at each location
  mutate(lon = st_coordinates(.)[,1],
         lat = st_coordinates(.)[,2]) %>%
  st_drop_geometry()
# write.csv(both_df, file = "/home/hbronnvik/Documents/chapter3/mapboxable_df1.csv")
## predict current maps using migration (https://doi.org/10.1007/s10980-016-0347-0)


check <- both_df %>% 
  mutate(dem = terra::extract(srtm_m, ., method = "bilinear")$DEM,
         flyway = ifelse(flyway == "BMF", "East", "West"))
ms <- check %>% 
  group_by(flyway) %>% 
  summarize(m = median(sqrt(dem), na.rm = T)) %>% 
  st_drop_geometry()
ggplot(check, aes(flyway, sqrt(dem))) +
  geom_hline(data = ms, aes(yintercept = m, color = flyway), lty = 2) +
  geom_violin(aes(fill = flyway)) +
  geom_boxplot(width = 0.33, alpha = 0.75) +
  scale_fill_manual(values = c("#f27e71", "#0098b0")) +
  scale_color_manual(values = c("#f27e71", "#0098b0")) +
  labs(x = "Flyway", y = "sqrt elevation (m)", fill = "", color = "")

# ggplot(check, aes(flyway, sqrt(dem))) + 
#   ## add half-violin from {ggdist} package
#   ggdist::stat_halfeye(
#     aes(fill = flyway),
#     ## custom bandwidth
#     adjust = .5, 
#     ## adjust height
#     width = .6, 
#     ## move geom to the right
#     justification = -.2, 
#     ## remove slab interval
#     .width = 0, 
#     point_colour = NA
#   ) +
#   # geom_jitter(alpha = 0.02, width = 0.1) + 
#   geom_boxplot(
#     aes(color = flyway),
#     width = .15, 
#     fill = NA,
#     ## remove outliers
#     outlier.color = NA ## `outlier.shape = NA` or `outlier.alpha = 0` works as well
#   ) +
#   coord_flip(xlim = c(1, 2)) +
#   scale_fill_manual("Flyway", values = c("#f27e71", "#0098b0")) +
#   scale_color_manual("Flyway", values = c("#f27e71", "#0098b0")) +
#   labs(x = "", y = "sqrt elevation (m)")

# plot resistance maps for the different variables
resistance <- rast(list.files("/home/hbronnvik/Documents/chapter3/Raven/resist", full.names = T))
# project back to geographic coordinates
resistance <- project(resistance, "EPSG:4326")
resistance <- crop(resistance, ext(-20, 50, -30, 60))
ordered_names <- c("dem", "blh", "obs", "dem_blh", "dem_obs", "obs_blh", "blh_obs_blh")
resistance <- resistance[[ordered_names]]
plot(resistance)
# using gg will reduce the resolution, but allows some flexibility as far as paneling
# New facet label names for supp variable
fac_labs <- c("Elevation", "Boundary layer height", "Conspecifics", "Elev. & BLH", 
              "Elev. & Conspecifics", "BLH & Conspecifics", "Elev., BLH, & Conspecifics")
names(fac_labs) <- c("dem", "blh", "obs", "dem_blh", "dem_obs", "obs_blh", "blh_obs_blh")

# png("/home/hbronnvik/Documents/chapter3/figures/resistance_maps.png",
#     height = 8.2, width = 11.7, units = "in", res = 400)
ggplot() +
  tidyterra::geom_spatraster(data = resistance) +
  tidyterra::scale_fill_grass_c("Resistance", palette = "grey") +
  # scale_fill_gradientn("Resistance", colors = colfunc(100), 
  #                      na.value = "transparent") +
  scale_x_continuous(expand = c(0, 0), labels = function(x) as.character(x)) +
  scale_y_continuous(expand = c(0, 0), labels = function(x) as.character(x)) +
  labs(x = "Longitude", y = "Latitude") +
  facet_wrap(~lyr, labeller = labeller(lyr = fac_labs), ncol = 4) +
  theme(strip.background = element_rect(linewidth = 0),
        axis.text = element_text(size = 11),
        legend.position = c(0.9, 0.2),
        legend.key.height = unit(1, "cm"))
# dev.off()

### --------------------------------
curmaps <- list.files("/home/hbronnvik/Documents/chapter3/Raven/Current", pattern = ".asc", full.names = T)
curmaps <- curmaps[grepl("(?=.*curconn)(?=.*\\.asc)", curmaps, perl = T)]

# do some plots to check the appearance
check <- lapply(curmaps, rast)
check <- rast(check)

# fix the variable names
names(check) <- gsub("curconn|curmap|res", "", names(check))
names(check) <- gsub("_", " ",
                     sub("dem", "DEM", 
                         sub("asp", "Aspect", 
                             sub("food", "Forage",
                                 sub("E_", " East",
                                     sub("W_", " West",
                                         sub("blh", "Max. BLH", 
                                             sub("mdh", "Med. BLH",
                                                 sub("obs", "Conspecifics",
                                                     names(check))))))))))
ord <- data.frame(name = names(check),
                  len = nchar(names(check))) %>% 
  arrange(len)
check <- terra::subset(check, ord$name)

# re-project, mask, and crop
# the mask used to allow barrier crossings but save computation of the seas
world_b <- rast("/home/hbronnvik/Documents/chapter3/barrier_mask.tif")
check <- mask(check, world_b)
check <- project(check, "EPSG:4326")
check <- crop(check, ext(-20, 50, -30, 60))

# obtain a map for orientation
library(rnaturalearth)
world <- ne_coastline(scale = "large", returnclass = "sf")

fac_labs <- c("Elevation", "Elevation", "Boundary layer height", "Boundary layer height", 
              "Elev. & BLH", "Elev. & BLH", "Conspecifics", "Conspecifics", 
              "Elev. & Consp.", "Elev. & Consp.", "BLH & Consp.", "BLH & Consp.", 
              "Elev., BLH, & Consp.", "Elev., BLH, & Consp.")
names(fac_labs) <- names(check)

cur_p <- lapply(1:nlyr(check), function(x){
  p <- ggplot() +
    tidyterra::geom_spatraster(data = sqrt(check[[x]])) +
    geom_sf(data = world, color = "grey30") +
    scale_fill_gradientn("sqrt current", colors = colfunc(100),
                         na.value = "transparent") +
    coord_sf(xlim = c(-20, 40), ylim = c(-30, 60)) +
    scale_x_continuous(expand = c(0, 0), labels = function(x) as.character(x)) +
    scale_y_continuous(expand = c(0, 0), labels = function(x) as.character(x)) +
    labs(x = "Longitude", y = "Latitude") +
    facet_wrap(~lyr, labeller = labeller(lyr = fac_labs), ncol = 4) +
    theme_void() + 
    theme(legend.position = "top",
          legend.title = element_text(size = 10),
          strip.text = element_text(size = 12))
  return(p)
})
# png("/home/hbronnvik/Documents/chapter3/figures/current_maps_west.png",
#     height = 8.2, width = 11.7, units = "in", res = 400)
ggpubr::ggarrange(cur_p[[2]], cur_p[[4]], cur_p[[8]], cur_p[[6]], 
                  cur_p[[10]], cur_p[[12]], cur_p[[14]], nrow = 2, ncol = 4)
# dev.off()
# png("/home/hbronnvik/Documents/chapter3/figures/current_maps_east.png",
#     height = 8.2, width = 11.7, units = "in", res = 400)
ggpubr::ggarrange(cur_p[[1]], cur_p[[3]], cur_p[[7]], cur_p[[5]], 
                  cur_p[[9]], cur_p[[11]], cur_p[[13]], nrow = 2, ncol = 4)
# dev.off()

# plot
# png("/home/hbronnvik/Documents/chapter3/figures/current_maps_east.png",
#     height = 8.2, width = 11.7, units = "in", res = 400)
# # make space for each subplot
# par(mfrow = c(2, 3))
# for (i in grep("East", names(check))) {
#   # remove water before plotting
#   pc <- check[[i]]
#   # pc <- mask(check[[i]], world_b)
#   plot(pc, col = colfunc(1000), main = names(check[[i]]))
#   lines(world, col = "gray20")
# }
# dev.off()
# png("/home/hbronnvik/Documents/chapter3/figures/current_maps_west.png",
#     height = 8.2, width = 11.7, units = "in", res = 400)
# # make space for each subplot
# par(mfrow = c(2, 3))
# for (i in which(!grepl("East", names(check)))) {
#   # remove water before plotting
#   pc <- mask(check[[i]], world_b)
#   plot(pc, col = colfunc(1000), main = names(check[[i]]))
#   lines(world, col = "gray20")
# }
# dev.off()
# par(mfrow = c(1,1))

# ggplot() +
#   tidyterra::geom_spatraster(data = check) +
#   tidyterra::scale_fill_grass_c(palette = "haxby") +
#   labs(x = "Longitude", y = "Latitude", fill = "Current") +
#   scale_x_continuous(expand = c(0, 0)) +
#   scale_y_continuous(expand = c(0, 0)) +
#   facet_wrap(~factor(lyr, levels = c("DEM", "Aspect", "Forage", "DEM X Aspect",
#                                      "DEM X Forage", "DEM X Aspect X Forage"))) +
#   theme_classic() +
#   theme(axis.text = element_text(color = "black", size = 10),
#         text = element_text(size = 15),
#         strip.background = element_rect(color = "white"),
#         strip.text = element_text(size = 15))


# count number of observations in each raster cell on each flyway
west_ml <- west_ml %>% 
  mutate(cell = terra::extract(rast(curmaps[1]), ., cells = T)$cell) %>%
  group_by(cell) %>% 
  mutate(n_obs_cell = length(unique(trackID))) %>% 
  ungroup()
east_ml <- east_ml %>% 
  mutate(cell = terra::extract(rast(curmaps[1]), ., cells = T)$cell) %>%
  group_by(cell) %>% 
  mutate(n_obs_cell = length(unique(trackID))) %>% 
  ungroup()

# annotate each flyway separately and make two big data frames for modeling
annotations <- lapply(c("E_", "W_"), function(f){
  if(f == "E_"){
    # the eastern current maps
    rel_curmaps <- curmaps[grepl("E_", curmaps)]
    # the eastern data
    ad_half <- east_ml
  }else{
    # the eastern current maps
    rel_curmaps <- curmaps[grepl("W_", curmaps)]
    # the western data
    ad_half <- west_ml
  }
  # each flyway
  flyway <- lapply(rel_curmaps, function(x){
    curmap <- rast(x)
    # find all the values of the current greater than 0
    v1 <- as.numeric(values(curmap, na.rm = T))
    v1 <- v1[v1>0]
    ecdf_func <- ecdf(v1)
    # annotate
    ad <- ad_half %>% 
      # add on the current values at each location
      mutate(current = terra::extract(curmap, .)[,2],
             # add on the percentile at each location
             cent = ecdf_func(current)*100,
             location.long = st_coordinates(.)[,1],
             location.lat = st_coordinates(.)[,2]) %>% 
      st_drop_geometry() %>% 
      dplyr::select(location.long, location.lat, timestamp, individual.id,
                    cent, current)
    colnames(ad)[grepl("current", colnames(ad))] <- sub("_curmap", "", sub("curconn", "current", names(curmap)))
    colnames(ad)[grepl("cent", colnames(ad))] <- sub("_curmap", "", sub("curconn", "percentile", names(curmap)))
    return(ad)
  })
  # combine the annotations on the original data
  ad_half <- ad_half %>% 
    left_join(flyway[[1]]) %>% 
    left_join(flyway[[2]]) %>% 
    left_join(flyway[[3]]) %>% 
    left_join(flyway[[4]]) %>% 
    left_join(flyway[[5]]) %>% 
    left_join(flyway[[6]]) %>% 
    left_join(flyway[[7]])
  return(ad_half)
})

# save out
# saveRDS(annotations, file = "/home/hbronnvik/Documents/chapter3/fall_a_data_cells_current_cent.rds")
annotations <- readRDS("/home/hbronnvik/Documents/chapter3/fall_a_data_cells_current_cent.rds")

# do two comparisons: a model predicting percentile of the current using model ID,
# and simply finding the average percentile produced by each model
outputs <- lapply(annotations, function(m){
  if("battery.charge.percent" %in% names(m)){
    m <- m %>% 
      dplyr::select(-battery.charge.percent)
  }
  # reshape the data to fit a model
  m_data <- m %>% 
    mutate(ID = as.factor(individual.id),
           trackID = as.factor(trackID)) %>% 
    # add on relocation to have a record of movements in line with current
    group_by(trackID) %>% 
    arrange(timestamp) %>% 
    mutate(reloc = 1:n()) %>% 
    ungroup() %>% 
    dplyr::select(ID, trackID, reloc, contains("cent")) %>% 
    st_drop_geometry() %>% 
    tidyr::pivot_longer(cols = contains("cent"), names_to = "model", values_to = "percentile")
  # make a factor order to determine the reference level for the model
  ord <- data.frame(name = unique(m_data$model),
                    len = nchar(as.character(unique(m_data$model)))) %>% 
    arrange(len)
  # order the factor
  m_data <- m_data %>%
    mutate(model = factor(model, levels = ord$name))
  # testing the effect of model on the response variable current while controlling for 
  # the random effect of individual variability and repeated locations.
  model <- lme4::lmer(percentile~model+(1|trackID/reloc), data = m_data)
  graph <- confint(model, method = "Wald") %>%
    as.data.frame() %>%
    rownames_to_column(var = "Model") %>%
    filter(!grepl("sig|cept", Model)) %>%
    left_join(summary(model)$coefficients %>%
                as.data.frame() %>%
                rownames_to_column(var = "Model")) %>%
    rename(lower = "2.5 %",
           upper = "97.5 %") %>% 
    mutate(Model = gsub("model|percentile|res", "", Model),
           Model = gsub("_", " ",
                        sub("dem", "DEM", 
                            sub("asp", "Aspect", 
                                sub("food", "Forage",
                                    sub("E_", " East",
                                        sub("W_", " West",
                                            sub("blh", "BLH", 
                                                sub("mdh", "Med. BLH",
                                                    sub("obs", "Conspecifics",
                                                        Model))))))))))
  
  avgs <- m_data %>% 
    group_by(model) %>% 
    summarize(Mean = round(mean(percentile), 1),
              Std.Dev. = round(sd(percentile), 1)) %>% 
    mutate(Model = gsub("percentile|res", "", model),
           Model = gsub("_", " ",
                        sub("dem", "DEM", 
                            sub("asp", "Aspect", 
                                sub("food", "Forage",
                                    sub("E_", " East",
                                        sub("W_", " West",
                                            sub("blh", "BLH", 
                                                sub("mdh", "Med. BLH",
                                                    sub("obs", "Conspecifics",
                                                        Model)))))))))) %>% 
    arrange(desc(Mean)) %>% 
    dplyr::select(Model, Mean, Std.Dev.)
  out <- avgs %>% 
    left_join(graph)
  return(out)
}) %>% 
  reduce(rbind) %>%
  mutate(flyway = ifelse(grepl("W", Model), "EAF", "BMF"))


# plot the comparisons
outs <- lapply(split(outputs, outputs$flyway), function(x){
  x <- x %>% 
    mutate(Model = substring(Model, 1, nchar(Model)-1))
  if(unique(x$flyway) == "BMF"){
    col <- "#0098b0"
  }else{
    col <- "#f27e71"
  }
  o <- x %>% 
    ggplot(aes(Estimate, Model)) +
    geom_vline(xintercept = 0, lty = 2) +
    geom_pointrange(aes(xmin = lower, xmax = upper), color = col) +
    labs(y = "", subtitle = unique(x$flyway)) +
    annotate(geom = 'table',
             x=ifelse(unique(x$flyway) == "BMF", 5, -3),
             y=2,
             label=list(x %>% 
                          dplyr::select(Model, Mean, Std.Dev.) %>% 
                          arrange(desc(Mean))),
             table.theme = ttheme_gtlight, 
             size = 5) 
  return(o)
})
# png("/home/hbronnvik/Documents/chapter3/figures/curr_comp_lmer_12x.png",
#     height = 8.2, width = 11.7, units = "in", res = 400)
ggpubr::ggarrange(outs[[2]], outs[[1]], nrow = 2)
# dev.off()

# png("/home/hbronnvik/Documents/chapter3/figures/curr_EAF_lmer_6x.png",
#     height = 6.5/2, width = 8, units = "in", res = 400)
outputs %>% 
  filter(flyway == "EAF") %>% 
  mutate(Model = gsub("W", "", Model),
         Estimate = ifelse(Model == "  DEM", 0, Estimate)) %>% 
  ggplot(aes(Estimate, forcats::fct_reorder(Model, Estimate))) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_pointrange(aes(xmin = lower, xmax = upper), color = "#f27e71", linewidth = 1.75) +
  labs(y = "")
# dev.off()
