### Comparing individuals relative to resistance
### Hester Bronnvik
### hbronnvik@ab.mpg.de
### 2025-06-18

library(sf)
library(terra)
library(dplyr)
library(purrr)
library(tibble)
library(lubridate)
library(rnaturalearth)
library(ggplot2)
library(ggpmisc)
theme_set(theme_classic()+
            theme(axis.text = element_text(color = "black", size = 12), 
                  text = element_text(size = 15)))
colfunc <- colorRampPalette(c("#0081A7", "#0098b0", "#00AFB9", "#7fc4b8", "#FED9B7", "#f7a58f", "#f27e71", "#F07167", "#ee5e53"))

## Now that we know which resistance surface performs best for predicting the flyway,
## we can use that resistance to look into individual variation.

# get the resistance map
resist <- rast("~/Documents/chapter3/Raven/resist/res_dem.tif")
# get the current map of the best model
curmap <- list.files("/home/hbronnvik/Documents/chapter3/Raven/Current", pattern = ".asc", full.names = T)
curmap <- curmap[grepl("demW", curmap)]
curmap <- rast(curmap)
names(curmap) <- "current"
# find all the values of the current greater than 0
v1 <- as.numeric(values(curmap, na.rm = T))
v1 <- v1[v1>0]
hist(v1)
# make a cumulative density function to estimate percentile later
ecdf_func <- ecdf(v1)

# get the actual locations used by white storks
# the EAF days classified as migratory, start and end migration dates, and outcome
w_metad <- readRDS("/home/hbronnvik/Documents/chapter3/migration_dates_2025-03-08.rds") %>% 
  filter(season == "fall")
# the EAF location data during migration
west_ml <- readRDS("/home/hbronnvik/Documents/chapter3/migration_locations_40km70km15daySpeed_2025-03-08.rds") %>% 
  dplyr::select(individual.id, trackID, location.long, location.lat, timestamp, 
                ground.speed, height.above.ellipsoid, gps.dop, heading, ground_speed_15, 
                date, daily_dist, daily_direction, track_displacement, journey_number, 
                total_journeys, track_status)

# only migratory flight locations
# filter ground speeds and daily distances (speed is scalar and a few missing data can affect it dramatically)
west_ml <- west_ml %>% 
  mutate(id_date = paste(individual.id, date, sep = "_")) %>% 
  filter(id_date %in% w_metad$id_date & ground_speed_15 > 2) %>%
  # make it spatial
  st_as_sf(coords = c("location.long", "location.lat"), crs = 4326) %>% 
  mutate(location.long = st_coordinates(.)[,1],
         location.lat = st_coordinates(.)[,2]) %>% 
  # match to the current map
  st_transform(crs(curmap)) %>%
  # annotate the movement
  mutate(current = terra::extract(curmap, .)$current,
         # add on the percentile at each location
         cent = ecdf_func(current)*100,
         resistance = terra::extract(resist, .)$dem)

# explore some potential covariates
west_ml %>% 
  # is there individual variation?
  filter(journey_number == 1 & track_status == "complete") %>% 
  ggplot(aes(resistance, forcats::fct_reorder(as.factor(individual.id), track_displacement), fill = track_displacement)) +
  ggridges::geom_density_ridges2() +
  labs(x = "Percentile of the current", y = "ID", fill = "Distance")

west_ml %>% 
  # did juvies that flowed with current succeed more?
  filter(journey_number == 1) %>% 
  mutate(fac = ifelse(total_journeys>1, "1 migration", "2+ migrations")) %>% 
  ggplot(aes(cent, fill = fac)) +
  geom_boxplot(alpha = 0.5) +
  labs(x = "Percentile of the current", y = "Density", fill = "Success")

west_ml %>% 
  mutate(fac = ifelse(journey_number==1, "Migration 1", "Migration 2+")) %>% 
  ggplot(aes(cent, fill = track_status)) +
  geom_density(alpha = 0.5) +
  labs(x = "Percentile of the current", y = "Density", fill = "Track\nstatus") +
  facet_wrap(~fac)

west_ml %>% 
  # did adults flow with more current?
  ggplot(aes(log(current), fill = as.factor(journey_number))) +
  geom_density(alpha = 0.5) +
  scale_fill_manual("Age", values = colfunc(9)) +
  labs(x = "Cumulative current", y = "Density")

west_ml %>% 
  mutate(age = factor(journey_number, levels = c(1:9))) %>% 
  arrange(age) %>% 
  # where did the current flow?
  ggplot(aes(location.lat, log(current), color = as.factor(journey_number))) +
  geom_point(alpha = 0.25) +
  scale_color_manual("Age", values = colfunc(9)) +
  scale_x_continuous(breaks = seq(10, 50, 10), 
                     labels = c("10\nGambia", "20\nSahara", "30\nAtlas", "40\nCentral Spain", "50\nFrankfurt")) +
  labs(x = "Latitude", y = "log cumulative current")

west_ml %>% 
  filter(location.lat >= 36 & track_status == "complete") %>% 
  # did adults flow with more current in Europe?
  ggplot(aes(log(current), fill = as.factor(journey_number))) +
  geom_density(alpha = 0.5) +
  scale_fill_manual("Age", values = colfunc(9)) +
  labs(x = "log cumulative current", y = "Density")

res_dens <- west_ml %>% 
  filter(location.lat >= 36) %>% 
  ggplot(aes(resistance, fill = as.factor(journey_number))) +
  geom_density(alpha = 0.5) +
  scale_fill_manual("Age", values = colfunc(9)) +
  coord_flip() +
  labs(x = "Resistance", y = "Density") +
  facet_wrap(~track_status)
res_vios <- west_ml %>% 
  filter(location.lat >= 36) %>% 
  ggplot(aes(as.factor(journey_number), resistance, fill = as.factor(journey_number))) +
  geom_violin() +
  geom_boxplot(alpha = 0.5, fill = "white") +
  scale_fill_manual("Age", values = colfunc(9)) +
  labs(y = "Resistance", x = "Age") +
  facet_wrap(~track_status)
res_lats <- west_ml %>% 
  filter(location.lat >= 34 & track_status == "complete") %>% 
  mutate(age = factor(journey_number, levels = c(1:9))) %>% 
  arrange(age) %>% 
  # where is the resistance?
  ggplot(aes(location.lat, resistance, color = as.factor(journey_number))) +
  geom_point(alpha = 0.25) +
  geom_smooth(color = "gray50") +
  scale_color_manual("Age", values = colfunc(9)) +
  scale_x_reverse(breaks = c(36, 40, 45, 50),
                  labels = c("36\nGibraltar", "40\nCentral Spain", "45\nLyon", "50\nFrankfurt")) +
  labs(x = "Latitude", y = "Resistance")
ggpubr::ggarrange(res_lats, res_dens, nrow = 1)

# how are current and resistance related?
ggplot(west_ml, aes(resistance, log(current))) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", aes(color = as.factor(journey_number)))  +
  scale_color_manual("Age", values = colfunc(9))

# how are current and distance related?
pr <- west_ml %>% 
  st_drop_geometry() %>% 
  filter(current > 0) %>%  
  filter(track_status == "complete" & location.lat >= 36) %>% 
  group_by(trackID) %>% 
  mutate(total_current = sum(log(current))) %>% 
  slice(1) %>% 
  ungroup() %>% 
  dplyr::select(total_current, track_displacement) %>% 
  tidyr::drop_na(total_current) %>% 
  cor()
dist_curr <- west_ml %>% 
  st_drop_geometry() %>% 
  filter(current > 0) %>%  
  filter(track_status == "complete") %>% 
  group_by(trackID) %>% 
  mutate(total_current = sum(log(current))) %>% 
  slice(1) %>% 
  ungroup() %>% 
  ggplot(aes(track_displacement, total_current)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = F, aes(color = as.factor(journey_number)))  +
  scale_color_manual("Age", values = colfunc(9)) +
  labs(x = "Total distance", y = "Sum cumulative current", subtitle = paste0("Pearson's r = ", round(as.data.frame(pr)[2,1], 3)))

# west_ml %>% 
#   st_drop_geometry() %>% 
#   filter(current > 0) %>%  
#   filter(track_status == "complete") %>% 
#   group_by(trackID, date) %>% 
#   summarize(daily_current = sum(log(current)),
#             daily_dist = unique(daily_dist),
#             age = unique(journey_number)) %>%  
#   ungroup()

# does moving through Gibraltar determine current?
corr_curr <- west_ml %>% 
  st_drop_geometry() %>% 
  filter(current > 0 & track_status == "complete" & location.lat >= 35.5) %>% 
  group_by(trackID) %>%  
  summarize(corridor = ifelse(T %in% unique(between(location.lat, 34, 37.5)), T, F),
            total_current = sum(log(current)),
            age = unique(journey_number)) %>% 
  ggplot(aes(total_current, fill = corridor)) +
  geom_density(alpha = 0.75) +
  scale_fill_manual(values = c("#f7a58f", "#ee5e53")) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.0025)) +
  labs(x = "Sum cumulative current North of 35.5°", y = "Density", fill = "Passed south\nof 37.5°")
age_corr <- west_ml %>%
  st_drop_geometry() %>% 
  filter(current > 0 & track_status == "complete" & location.lat >= 35.5) %>% 
  group_by(trackID) %>%  
  summarize(corridor = ifelse(T %in% unique(between(location.lat, 34, 37.5)), T, F),
            total_current = sum(log(current)),
            age = unique(journey_number)) %>% 
  group_by(age) %>% 
  # the total number of tracks from each age
  mutate(tracks_per_age = length(unique(trackID))) %>% 
  group_by(corridor, age) %>% 
  # the proportions of individuals in each age group in each distance group
  summarize(inds = length(unique(trackID))/unique(tracks_per_age),
            n = unique(tracks_per_age)) 
age_corr <- age_corr %>% 
  ggplot(aes(as.factor(age), inds, fill = corridor)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = c("#f7a58f", "#ee5e53")) +
  scale_x_discrete(expand = c(0, 0), labels = paste0(age_corr$age, "\nn = ", age_corr$n)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Passed south of 37.5°", y = "Proportion", fill = "Age")
# add coastlines
world <- ne_countries(scale = "large", returnclass = "sf")[1]
world <- st_transform(world, crs = crs(west_ml))
data.frame(x = c(0, 0), y = c(34, 37.5)) %>% 
  st_as_sf(coords = c("x", "y"), crs = 4326) %>% 
  st_transform(crs(west_ml))
sub_map <- west_ml %>%
  arrange(journey_number) %>% 
  ggplot() +
  geom_sf(data = world, color = "gray50", fill = "gray50") +
  geom_sf(aes(color = as.factor(journey_number)), alpha = 0.5, size = 0.15) +
  scale_color_manual(values = colfunc(9)) +
  scale_y_continuous(breaks = c(30, 34, 37.5, 40, 45, 50, 55)) +
  coord_sf(xlim = c(-901856.37732783, 1545226.84696131), 
           ylim = c(3703706.1360406, 6479941.38506427)) +
  geom_hline(yintercept = c(4103741, 4503046), lty = 2) +
  labs(x = "Longitude", y = "Latitude", color = "Age")

# png("/home/hbronnvik/Documents/chapter3/figures/corridor_current.png",
#     width = 6.5, height = 8, units = "in", res = 400)
ggpubr::ggarrange(ggpubr::ggarrange(sub_map, dist_curr,
                                    labels = c("A", "B"),
                                    common.legend = T, legend = "right"),
                  corr_curr, age_corr,
                  labels = c(NA, "C", "D"),
                  nrow = 3)
# dev.off()



# if adults are less responsive to difficulty (Ch. 1 & Ch. 2), age should
# show a positive effect on resistance and a negative effect on current

df4_mod <- west_ml%>% 
  # only compare where both adults and juveniles are observed
  filter(location.lat >= 40 & track_status == "complete" & journey_number <= 5) %>% 
  # useful information
  dplyr::select(individual.id, location.long, location.lat, date, journey_number, track_status, 
                cent, resistance, current) %>% 
  st_drop_geometry() %>% 
  # reduce skew
  mutate(log_curr = scale(log(current+1)),
         log_res = scale(log(resistance)),
         log_prob = scale(log(1-resistance)),
         age = scale(journey_number),
         ID = as.factor(individual.id)) 

# check the sample size
df4_mod %>% 
  st_drop_geometry() %>% 
  group_by(journey_number) %>% 
  summarize(n = length(unique(individual.id)))

cur_mod <- lme4::lmer(log_curr~age+(1|ID), data = df4_mod)
res_mod <- lme4::lmer(log_res~age+(1|ID), data = df4_mod)
cent_mod <- lme4::lmer(cent~age+(1|ID), data = df4_mod)
prob_mod <- lme4::lmer(log_prob~age+(1|ID), data = df4_mod)

coefs <- lapply(c(cur_mod, res_mod, cent_mod, prob_mod), function(x){
  pred <- sub("(|)", "", formula(x)[2])
  pred <- ifelse(grepl("res", pred), "Yr. 1\nresistance",
                 ifelse(grepl("prob", pred), "Yr. 1\nselection",
                        ifelse(grepl("curr", pred), "Cumulative\ncurrent",
                               ifelse(grepl("cent", pred), "Percentile\ncurrent",
                                      pred))))
  graph <- confint(x, method = "Wald") %>%
    as.data.frame() %>%
    rownames_to_column(var = "Model") %>%
    filter(!grepl("sig|cept", Model)) %>%
    left_join(summary(x)$coefficients %>%
                as.data.frame() %>%
                rownames_to_column(var = "Model")) %>%
    rename(lower = "2.5 %",
           upper = "97.5 %") %>% 
    mutate(Model = pred)
  return(graph)
}) %>% reduce(rbind)

# In Europe, age has negative effects on current, percentile of the current, and resistance

# png("/home/hbronnvik/Documents/chapter3/figures/ageXpred.png",
#     height = 6.5/2, width = 8, units = "in", res = 400)
ggplot(coefs, aes(Estimate, Model)) +
  geom_vline(lty = 2, xintercept = 0) +
  geom_pointrange(aes(xmin = lower, xmax = upper), color = "#f27e71") +
  labs(x = "Fixed effect estimate and 95% CI", y = "")
# dev.off()

# visualize
# add color
age_cols <- data.frame(journey_number = 1:5, col = colfunc(5))
df4_mod <- df4_mod %>% 
  left_join(age_cols)
# add coastlines
world <- ne_coastline(scale = "large", returnclass = "sf")[1]
world <- st_crop(world, c(xmin = -15, xmax = 17, ymin = 33, ymax = 55))
# reproject the map
wgs_curmap <- project(curmap, "EPSG:4326")
wgs_curmap <- crop(wgs_curmap, c(-15, 17, 33, 55))
wgs_curmap <- classify(wgs_curmap, cbind(50, 600, 50))
plot(wgs_curmap, col = rev(terrain.colors(300)))
lines(world, col = "gray30")
points(df4_mod$location.long, df4_mod$location.lat, col = df4_mod$col, cex = 0.25)
# text(x = par("usr")[2] + 0.5, y = mean(par("usr")[3:4]), 
#      "Legend Title", xpd = TRUE)
text(x = 19, y = 56.60897, "Cumulative\ncurrent", xpd = TRUE)

# look at only the "best" model and include age
library(glmmTMB)
library(performance)
## 2. access the migration data that have been burst and their alternative points (annotated)
a_data <- readRDS("/home/hbronnvik/Documents/chapter3/ua_data_60min_20250526.rds")

# the days classified as migratory, start and end migration dates, and outcome
metad <- readRDS("/home/hbronnvik/Documents/chapter3/migration_dates_2025-03-08.rds")
# the ages by track
ages <- metad %>% 
  group_by(trackID) %>% 
  slice(1) %>% 
  filter(track_status == "complete") %>% 
  dplyr::select(trackID, journey_number) %>% 
  rename(track = trackID,
         age = journey_number)

# prep the data
a_data <- a_data %>% 
  filter(track %in% ages$track) %>% 
  tidyr::separate(track, c("id", "season", "year"), sep = "_", remove = F) %>% 
  dplyr::select(-id, -year) %>% 
  # add on age
  left_join(ages) %>% 
  # filter(age <= 5) %>% 
  # only consider one season for now and use juvenile decisions
  filter(season == "fall") %>% 
  mutate(sqrt_dem = log(DEM+424),
         sqrt_sl = sqrt(step_length),
         age_z = scale(age)[,1],
         age = as.factor(age),
         stratum_ID = as.factor(stratum),
         stratum_ID = as.numeric(stratum_ID),
         individual.id = as.numeric(individual.id)) %>% 
  tidyr::drop_na(sqrt_dem_z)

# set up the model
# tmb_struc <- glmmTMB(used ~ sqrt_dem_z*age_z + sqrt_sl_z + ta_z + (1|stratum_ID) + (0 + sqrt_dem_z|individual.id) + (0 + age_z|individual.id),
#                      family = poisson,
#                      data = a_data, doFit = FALSE,
#                      #Tell glmmTMB not to change the first standard deviation, all other values are freely estimated (and are different from each other)
#                      map = list(theta = factor(c(NA, 1, 1))), # 3 is the n of random slopes
#                      #Set the value of the standard deviation of the first random effect (here (1|startum_ID)):
#                      start = list(theta = c(log(1e3), 0, 0))) #add a 0 for each random slope. in this case, 2
# # run the model
# TMB <- glmmTMB:::fitTMB(tmb_struc)
# summary(TMB)

grd_2x <- expand.grid(x = 1:9,
                      y = seq(from = quantile(a_data$sqrt_dem, 0.025, na.rm = T), 
                              to = quantile(a_data$sqrt_dem, 0.975, na.rm = T), length.out = 250)) %>%
  rename(age = x,
         sqrt_dem = y) %>% 
  mutate(#set other variables to their mean
    sqrt_sl = mean(a_data$sqrt_sl, na.rm = T),
    turning_angle = mean(a_data$turning_angle, na.rm = T),
    interaction = "dem_age")

set.seed(770)
n <- nrow(grd_2x)

new_data_only <- a_data %>%
  group_by(stratum_ID) %>% 
  slice_sample(n = 1) %>% #randomly selects one row (from each stratum)
  ungroup() %>% 
  slice_sample(n = n, replace = F) %>% #randomly select n of these strata. Make sure n is less than n_distinct(data$stratum)
  # only keep the columns that I need
  dplyr::select(c("stratum_ID", "individual.id")) %>% 
  bind_cols(grd_2x) %>% 
  #calculate z-scores. get the mean(center) and sd(scale) from the previous z transformation. for consistency
  mutate(age_z = (age - mean(a_data$age))/(sd(a_data$age)),
         sqrt_dem_z = (sqrt_dem - mean(a_data$sqrt_dem, na.rm = T))/(sd(a_data$sqrt_dem, na.rm = T)),
         sqrt_sl_z = (sqrt_sl - mean(a_data$sqrt_sl, na.rm = T))/(sd(a_data$sqrt_sl, na.rm = T)), 
         ta_z = (turning_angle - mean(a_data$turning_angle, na.rm = T))/(sd(a_data$turning_angle, na.rm = T)))

new_data <- a_data %>% 
  mutate(interaction = "OG_data") %>% 
  dplyr::select(names(new_data_only)) %>%  #only keep the columns that are necessary for the model
  bind_rows(new_data_only) %>% 
  #accoring to the predict.glmmTMB help file: "To compute population-level predictions for a given grouping variable 
  #(i.e., setting all random effects for that grouping variable to zero), set the grouping variable values to NA."
  mutate(stratum_ID = NA)#,
# individual.id = NA)

# now that we have the values to predict, run the model on them
preds <- predict(TMB, newdata = new_data, type = "link")

preds_pr <- new_data %>% 
  mutate(preds = preds) %>% 
  rowwise() %>%
  # using poisson inverse link
  mutate(probs = exp(preds)) %>% 
  ungroup()
# mutate(probs = gtools::inv.logit(preds)) #https://rpubs.com/crossxwill/logistic-poisson-prob

inter_preds <- preds_pr %>% 
  mutate(rel_prob = probs/max(preds_pr$probs, na.rm = T)) %>% 
  filter(interaction != "OG_data") 

# the whole population:
graph <- confint(TMB, method = "Wald") %>%
  as.data.frame() %>%
  rownames_to_column(var = "pred") %>%
  filter(!grepl("sig|cept|Std.Dev", pred)) %>%
  rename(lower = "2.5 %",
         upper = "97.5 %") %>% 
  mutate(pred = gsub("sqrt_dem_z", "log Elevation", 
                     gsub("age_z", "Age", 
                          gsub("sqrt_sl_z", "Step length",
                               gsub("ta_z", "Turning angle",
                                    pred)))),
         pred_label = c(
           "<span style='color:black'>log Elevation</span>",
           "<span style='color:black'>Age</span>",
           "<span style='color:grey40'>Step length</span>",
           "<span style='color:grey40'>Turning angle</span>",
           "<span style='color:black'>log Elevation:Age</span>"
         ))
coefs <- graph %>% 
  # it doesn't actually make sense to have age as a predictor alone because a 
  # bird can only be one age at a time, so only the interaction makes sense 
  # hence the huge error on the age alone
  filter(pred != "Age") %>% 
  ggplot(aes(Estimate, pred_label)) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_pointrange(aes(xmin = lower, xmax = upper), color = "#BD3786") +
  labs(y = "") +
  # parse HTML
  theme(axis.text.y = ggtext::element_markdown())
interactions <- ggplot(inter_preds, aes(as.factor(age), sqrt_dem, fill = rel_prob)) +
  geom_tile() +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  viridis::scale_fill_viridis("Selection\nprobability",
                              option = "C",
                              values = scales::rescale(c(0, 0.05, 0.2, 0.5, 1)),
                              limits = c(0, 0.55)) +
  labs(x = "Age", y = "log Elevation (m)")
# interactions <- ggplot(inter_preds, aes(y = exp(sqrt_dem), x = age, z = rel_prob)) +
#   geom_contour_filled() +
#   scale_fill_viridis_d(option = "C") +
#   scale_x_discrete(expand = c(0, 0)) +
#   scale_y_continuous(expand = c(0, 0))  +
#   labs(x = "Age", y = "Elevation (m)", fill = "Level")


# png("/home/hbronnvik/Documents/chapter3/figures/age_dem_ssf_contour.png",
#     height = 6.5, width = 8, units = "in", res = 400)
ggpubr::ggarrange(coefs, interactions, nrow = 2, labels = "AUTO",
                  heights = c(1, 3), common.legend = T, legend = "right")
# dev.off()

tmb_struc <- glmmTMB(used ~ sqrt_dem_z + sqrt_sl_z + ta_z + (1|stratum_ID) + (0 + sqrt_dem_z|age),
                     family = poisson,
                     data = a_data, doFit = FALSE,
                     #Tell glmmTMB not to change the first standard deviation, all other values are freely estimated (and are different from each other)
                     map = list(theta = factor(c(NA, 1))), # 3 is the n of random slopes
                     #Set the value of the standard deviation of the first random effect (here (1|startum_ID)):
                     start = list(theta = c(log(1e3), 0))) #add a 0 for each random slope. in this case, 2
# run the model
TMB <- glmmTMB:::fitTMB(tmb_struc)
summary(TMB)

fixed_effects <- confint(TMB, method = "Wald") %>%
  as.data.frame() %>%
  rownames_to_column(var = "pred") %>%
  filter(!grepl("sig|cept|Std.Dev", pred)) %>%
  rename(lower = "2.5 %",
         upper = "97.5 %") %>% 
  mutate(pred = gsub("sqrt_dem_z", "log Elevation",  
                     gsub("sqrt_sl_z", "Step length",
                          gsub("ta_z", "Turning angle",
                               pred))),
         pred_label = c(
           "<span style='color:black'>log Elevation</span>",
           "<span style='color:grey40'>Step length</span>",
           "<span style='color:grey40'>Turning angle</span>"
         ))
fxd <- fixed_effects %>% 
  ggplot(aes(Estimate, pred_label)) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_pointrange(aes(xmin = lower, xmax = upper), color = "#F07167") +
  labs(y = "", x = "Fixed effect estimate +/- 95% CI") +
  # parse HTML
  theme(axis.text.y = ggtext::element_markdown())

random_effects <- as.data.frame(ranef(TMB))

random_effects <- random_effects %>% 
  rename(age = grp) %>% 
  mutate(age = as.numeric(as.character(age)),
         lower = condval - condsd,
         upper = condval + condsd,
         sig = sign(lower)==sign(upper)) %>% 
  filter(term == "sqrt_dem_z")

rnd <- ggplot(random_effects, aes(x = condval, y = as.factor(age))) +
  geom_vline(xintercept = 0, linetype="dashed", 
             color = "gray50", linewidth = 0.5) +  
  geom_pointrange(aes(xmin = lower, xmax = upper), size = 0.8, color = "#F07167") +
  labs(x = "Difference from fixed effect estimate", y = "Age")

# png("/home/hbronnvik/Documents/chapter3/figures/dem_ssf_age_as_ranef.png",
#     height = 6.5, width = 8, units = "in", res = 400)
ggpubr::ggarrange(fxd, rnd, nrow = 2, labels = "AUTO")
# dev.off()


# make new data to predict across
# pred_grid <- expand.grid(
#   # the middle 95% of the elevations 
#   sqrt_dem_z = seq(from = quantile(a_data$sqrt_dem_z, 0.025, na.rm = T), 
#                    to = quantile(a_data$sqrt_dem_z, 0.975, na.rm = T), length.out = 250),
#   # the ages in the data
#   age_z = unique(a_data$age_z)) %>% 
#   # set other variables
#   mutate(sqrt_sl_z = mean(a_data$sqrt_sl_z, na.rm = T),  # keep other covariates fixed
#          ta_z = mean(a_data$ta_z, na.rm = T),
#          stratum_ID = NA,                      # these are random slopes, leave as NA
#          individual.id = sample(a_data$individual.id, n(), replace = F))
# # get the predicted values
# pred_grid$fit <- predict(TMB, newdata = pred_grid, type = "link")
# pred_grid <- pred_grid %>% 
#   # unscale the predictors
#   mutate(sqrt_dem = sqrt_dem_z * sd(a_data$sqrt_dem) + mean(a_data$sqrt_dem),
#          age = age_z * sd(a_data$age) + mean(a_data$age)) %>% 
#   rowwise() %>% 
#   # turn predictions to probabilities
#   mutate(probs = gtools::inv.logit(fit)) %>% 
#   ungroup() %>% 
#   mutate(rel_prob = probs/max(probs, na.rm = T))
# 
# ggplot(pred_grid, aes(as.factor(age), sqrt_dem, fill = rel_prob)) +
#   geom_tile() +
#   labs(x = "Age", y = "log DEM", fill = "Selection\nprobability")


## if deaths occur in areas of high conspecific density, they should occure in areas of relatively high current

# required information
load("/home/hbronnvik/Documents/storkSSFs/loginStored.RData") # Movebank credentials
studies <- c(24442409, 212096177, 76367850, 21231406, 1176017658, 173641633) # Movebank study IDs

# movebank_download_deployment
# get the current names and deployment dates of all birds
metad <- lapply(studies, function(x){
  md <- getMovebankReferenceTable(study = x, login = loginStored) %>%
    tidyr::drop_na(animal_id) %>%
    filter(sensor_type_id == 653) %>%
    dplyr::select(animal_id, deploy_on_timestamp, animal_local_identifier, 
                  tag_local_identifier, tag_timestamp_end) %>%
    mutate(deploy_on_timestamp = date(deploy_on_timestamp), study = x) %>%
    rename(individual.id = animal_id)
  return(md)
}) %>% 
  reduce(rbind) %>% 
  # remove birds that were captured twice
  filter(!duplicated(individual.id)) %>% 
  mutate(tag_local_identifier = as.numeric(tag_local_identifier),
         # add a year of deployment to compare to metadata
         tag_year = year(deploy_on_timestamp)) 

# birds that died, their death locations, and death dates
deaths <- read.csv("/home/hbronnvik/Documents/chapter2/Copy of Metadata_storks_a.csv") %>% 
  # add a year of deployment to compare to metadata
  mutate(tag_year = year(as.Date(Capture.deployment.date))) %>% 
  # discard a study we don't need here
  filter(!Movebank.Study %in% c("LifeTrack White Stork Loburg", 
                                "White Stork Loburg 2014",
                                "LifeTrack White Stork Loburg 2022-23"))

# not all the birds actually have matching names
missing <- deaths %>% 
  filter(!Individual.ID...same.as.movebank %in% metad$animal_local_identifier)

# e.g. has parentheses around the + sign
missing$Individual.ID...same.as.movebank[2]
metad$animal_local_identifier[grepl("eobs 3017", metad$animal_local_identifier)]

# put the metadata individual IDs on the deaths using the combination of year of tagging and tag number
missing <- missing %>% 
  left_join(metad[, c("tag_local_identifier", "individual.id", "tag_year", "animal_local_identifier")], 
            join_by(Transmitter.ID...same.as.movebank == tag_local_identifier, tag_year)) %>% 
  dplyr::select(-animal_local_identifier)

# attach individual IDs to the death records using name
deaths <- deaths %>% 
  filter(Individual.ID...same.as.movebank %in% metad$animal_local_identifier) %>% 
  left_join(metad[, c("tag_local_identifier", "individual.id", "tag_year", "animal_local_identifier")], 
            join_by(Individual.ID...same.as.movebank == animal_local_identifier, tag_year)) %>% 
  dplyr::select(-tag_local_identifier) %>% 
  # add on the other records using tag and year
  rbind(missing)

# there are still a couple missing
table(is.na(deaths$individual.id))
# one was re-tagged and has two tag numbers
# the other is missing its tag in the metadata
deaths$individual.id[which(deaths$Individual.ID...same.as.movebank == "Sissi / DER AZ999 (eobs 8010)")] <- 930962485
deaths$individual.id[which(deaths$Individual.ID...same.as.movebank == "Langenau22-5 ABT75 (e-obs 9692)")]

# now that we have the IDs of all the dead birds, we can clean up the metadata a little
deaths <- deaths %>% 
  tidyr::drop_na(individual.id) %>% 
  dplyr::select(individual.id, sex, Captive.bred..Y.N.., Rehabilitated..Y.N., Certainty.in.transmitter.failure.vs.animal.death,
                when.died., Where.died..X.coord, Where.died..Y.coord, Cause.of.death.) %>% 
  rename(captive_bred = Captive.bred..Y.N..,
         rehabed = Rehabilitated..Y.N.,
         did_die = Certainty.in.transmitter.failure.vs.animal.death,
         dod = when.died.,
         death_lon = Where.died..X.coord,
         death_lat = Where.died..Y.coord,
         cod = Cause.of.death.) %>% 
  mutate(cod = ifelse(cod == "" | cod == " ", NA, cod),
         dod = ifelse(dod == "" | dod == " ", NA, dod),
         dod = as.Date(dod, format = "%m/%d/%Y"),
         sex = ifelse(sex == "" | sex == " ", NA, sex),
         did_die = ifelse(did_die == "" | did_die == " ", NA, did_die),
         # if death may be due to crowding, H1, if due to exhaustion, H2
         death_type = ifelse(cod %in% c("collision", "electrocution"), "H1",
                             ifelse(cod %in% c("exhaustion", "landfill", "starvation",
                                               "drowned"), "H2", "Other"))) %>% 
  tidyr::drop_na(dod)
# this leaves the location, time, and cause of death, the individual ID, and sex


# where did they die?
deaths %>% 
  # only birds dead of collisions and electrocutions
  # filter(death_type == "H1") %>% 
  ggplot(aes(death_lon, death_lat, color = cod)) +
  borders(xlim = c(-20, 20), ylim = c(10, 55), fill = "gray60") +
  geom_point()

death_sf <- deaths %>% 
  filter(death_type == "H1"|death_type == "H2") %>% 
  st_as_sf(coords = c("death_lon", "death_lat"), crs = 4326) %>% 
  st_transform(crs(curmap)) %>% 
  mutate(current = terra::extract(curmap, .)$current,
         # add on the percentile at each location
         cent = ecdf_func(current)*100,
         resistance = terra::extract(resist, .)$dem,
         death_lon = st_coordinates(.)[1], 
         death_lat = st_coordinates(.)[2])

death_sf %>% 
  st_drop_geometry() %>% 
  dplyr::select(death_lon, death_lat, cent, current, resistance) %>% 
  mutate(type = "death") %>% 
  rename(lon = death_lon,
         lat = death_lat) %>% 
  rbind(west_ml %>% 
          mutate(type = "life",
                 lon = st_coordinates(.)[1], 
                 lat = st_coordinates(.)[2]) %>% 
          st_drop_geometry() %>% 
          dplyr::select(lon, lat, cent, current, resistance, type)) %>% 
  ggplot(aes(type, log(current))) +
  geom_boxplot()

dead_4mod <- death_sf %>% 
  # st_drop_geometry() %>% 
  dplyr::select(death_lon, death_lat, cent, current, resistance) %>% 
  mutate(type = "death") %>% 
  rename(lon = death_lon,
         lat = death_lat) %>% 
  rbind(west_ml %>% 
          mutate(type = "life",
                 lon = st_coordinates(.)[1], 
                 lat = st_coordinates(.)[2]) %>% 
          # st_drop_geometry() %>% 
          dplyr::select(lon, lat, cent, current, resistance, type))

# randomization of the death date/ live date identifier
no.perm <- 5000
m_vals <- data.frame()
# look at medians of BLH
for (i in 1:no.perm) {
  rando <- mosaic::resample(dead_4mod, shuffled = "type")#, groups = ID)
  
  m <- median(rando$current[rando$type == "death"])
  m_vals[i, 1] <- i
  m_vals[i, 2] <- m
}

# true median of the death location current 
m <- median(dead_4mod$max_blh[dead_4mod$death_date == T])

ggplot(m_vals, aes(V2)) +
  geom_histogram(bins = 100, color = "gray40", fill = "gray50") +
  # geom_segment(data = segment_data, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_segment(x = m, y = 0, xend = m, yend = 700, color = "#EF6E64") +
  # scale_y_continuous(expand = c(0, 0), limits = c(0, 800)) +
  labs(x = "Randomized median", y = "Count")

# death_sf <- deaths %>% 
#   filter(death_type == "H1"|death_type == "H2") %>% 
#   st_as_sf(coords = c("death_lon", "death_lat"), crs = 4326) %>% 
#   st_transform(crs(curmap)) %>% 
#   st_buffer(40000) 
# 
# death_current <- terra::extract(curmap, death_sf)$current
# hist(death_current, breaks = 134)
# 
# obs_buffered <- west_ml  %>% 
#   dplyr::select(geometry) %>% 
#   st_buffer(1000) 
# life_current <- terra::extract(curmap, obs_buffered)$current

