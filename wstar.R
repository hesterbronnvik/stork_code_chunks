### How to estimate w* (the code)
### Hester Bronnvik  
### 2023-06-08

library(readr)
library(tidyverse)
library(ggplot2)
library(parallel)
library(terra)
library(rasterVis)

# a function to interpolate data
int <- function(x, x1, x2, y1, y2){
  num <- x - x1
  den <- x2 - x1
  mul <- y2 - y1
  y <- y1+((num/den)*(mul))
  return(y)
}
# a function to take the cube root
CubeRoot<-function(x){
  sign(x)*abs(x)^(1/3)
}
# a function to calculate w*
convective <- function(data) {
  
  df <- data # the data input
  z <- df[, "blh"] # the boundary layer height
  s_flux <- df[, "sflux"] # the surface sensible heat flux
  m_flux <- df[, "mflux"] # the moisture flux
  T_k_2m <- df[, "temp2m"] # the temperature at ground level
  p0 <- df[, "p0"] # the pressure at ground level
  q <- df[, "blh_humidity"] # the humidity at the boundary layer height
  p1 <- df[, "blh_pressure"] # the pressure at the boundary layer height
  T_K_blh <- df[, "blh_temperature"] # the temperature at the boundary layer height
  
  k <- 0.2854 # Poisson constant
  g <- 9.80665 # acceleration due to gravity
  T_c_2m <- T_k_2m - 273.15 # surface temperature in degrees Celsius
  c_p <- 1012 # the isobaric mass heat capacity of air at common conditions
  p <- 1.225 # the density of air at sea level and 273 K
  
  Theta_k_z <- T_K_blh*((p0/p1)^k)
  Thetav_k_z <- Theta_k_z*(1+(0.61*q))
  wT <- (s_flux * -1)/c_p/p # reverse the sign. ECMWF upward fluxes are negative
  wq <- (m_flux * -1)/p # reverse the sign
  
  wthetav <- wT + 0.61 * T_c_2m * wq
  
  w_star <- CubeRoot(g*z*(wthetav/Thetav_k_z))
  
  return(w_star)
  
}


# the location data
data <- readRDS("C:/Users/hbronnvik/Documents/storkSSFs/used_av_df_230612.rds") %>% 
  rename(long = "location.long",
         lat = "location.lat") %>% 
  group_by(stratum) %>% 
  slice(1:100) %>%
  ungroup() %>% 
  as.data.frame()


# loop through each individual point and annotate it with its environmental data

(start_time <- Sys.time())
annotated_data <- lapply(1:nrow(data), function(rown){
  print(paste0("Extracting environmental data for observation ", rown, " of ", nrow(data), "."))
  
  d_year <- year(data[rown, "timestamp"])
  stamp <- round_date(data[rown, "timestamp"], unit = "hour")
  
  ## retrieve data with ECMWF variables (downloaded with the script "ecmwf_data.R")
  # the height of the boundary layer
  file <- paste0("C:/Users/hbronnvik/Documents/storkSSFs/ecmwf_data/single/boundary_layer_height_land_", d_year, ".nc")
  blh <- try(terra::rast(file))
  # levelplot(r[[1]])
  # time(r)
  blh <- blh[[which(time(blh) == round_date(data[rown, "timestamp"], unit = "hour"))]]
  
  # the boundary layer height at the point of interest
  poi_blh <- terra::extract(blh, vect(data[rown, ], geom = c("long", "lat")))[,2]
  
  # the height of the pressure levels at the point of interest
  levels <- c(seq(from = 600, to = 1000, by = 25))
  levels <- levels[!levels %in% c(625, 675, 725)]

  heights <- lapply(levels, function(x){
    file <- paste0("C:/Users/hbronnvik/Documents/storkSSFs/ecmwf_data/pressure/Geopotential_", x, "_", d_year, ".nc")
    geo <- try(terra::rast(file))
    geo <- geo[[which(time(geo) == stamp)]]
    poi_geoH <- terra::extract(geo, vect(data[rown, ], geom = c("long", "lat")))[ ,2]/9.80665
    h <- data.frame(level = x, height = poi_geoH)
  }) %>% reduce(rbind) %>% as.data.frame() 
  
  # the levels above and below the PBL at the point of interest
  low_pressure <- heights %>%
    filter(height > poi_blh) %>%
    filter(height == min(height))
  high_pressure <- heights %>%
    filter(height < poi_blh) %>%
    filter(height == max(height))
  
  variables <- c("Specific humidity_", "Temperature_")
  levels <- c(low_pressure$level, high_pressure$level)
  
  at_pressure <- lapply(variables, function(v){
    lapply(levels, function(l){
      file <- paste0("C:/Users/hbronnvik/Documents/storkSSFs/ecmwf_data/pressure/", v, l, "_", d_year, ".nc")
      ras <- try(terra::rast(file))
      ras <- ras[[which(time(ras) == stamp)]]
      val <- terra::extract(ras, vect(data[rown, ], geom = c("long", "lat")))[ ,2]
      variable <- ifelse(v == "Specific humidity_", "q", "temp")
      df <- data.frame(variable = variable, level = l, value = val)
    }) %>% reduce(rbind)
  }) %>% reduce(rbind) 
  
  pressures <- data.frame(index = 1)
  pressures$blh_humidity <- int(poi_blh, high_pressure$height, low_pressure$height,
                                at_pressure$value[at_pressure$variable == "q" & at_pressure$level == high_pressure$level],
                                at_pressure$value[at_pressure$variable == "q" & at_pressure$level == low_pressure$level])
  
  
  pressures$blh_temperature <- int(poi_blh, high_pressure$height, low_pressure$height,
                                   at_pressure$value[at_pressure$variable == "temp" & at_pressure$level == high_pressure$level],
                                   at_pressure$value[at_pressure$variable == "temp" & at_pressure$level == low_pressure$level])
  
  
  pressures$blh_pressure <- int(poi_blh, high_pressure$height, low_pressure$height,
                                high_pressure$level, low_pressure$level)
  
  # now get the single level data for the point of interest
  variables <- c("2m_temperature_land_", "surface_pressure_land_", 
                 "instantaneous_surface_sensible_heat_flux_land_", "instantaneous_moisture_flux_land_")
  
  singles <- lapply(variables, function(v){
    file <- paste0("C:/Users/hbronnvik/Documents/storkSSFs/ecmwf_data/single/", v, d_year, ".nc")
    ras <- try(terra::rast(file))
    ras <- ras[[which(time(ras) == stamp)]]
    val <- terra::extract(ras, vect(data[rown, ], geom = c("long", "lat")))[ ,2]
    variable <- ifelse(v == "2m_temperature_land_", "temp2m", 
                       ifelse(v == "surface_pressure_land_", "p0", 
                              ifelse(v == "instantaneous_surface_sensible_heat_flux_land_", "sflux", "mflux")))
    df <- data.frame(variable = variable, value = val)
  }) %>% 
    reduce(rbind) %>% 
    pivot_wider(names_from = "variable", values_from = "value")
  
  df <- cbind(data[1, 2:ncol(data)], singles, pressures[, 2:4])
  df$blh <- poi_blh
  df
}) %>% reduce(rbind)
Sys.time() - start_time # 26.43 minutes for 100 rows

# estimate w*
annotated_data$w_star <- convective(data = annotated_data)

# Parallel processed, saves ~6 minutes on 100 row trial
# -------------------------------------------------------

file_list_s <- list.files("C:/Users/hbronnvik/Documents/storkSSFs/ecmwf_data/single/", full.names = T) 
file_list_p <- list.files("C:/Users/hbronnvik/Documents/storkSSFs/ecmwf_data/pressure/", full.names = T) 

cl <- makeCluster(16)
clusterExport(cl, c("file_list_s", "file_list_p", "int", "data"))
clusterEvalQ(cl, c(library(tidyverse), library(terra)))

# loop through each individual point and annotate it with its environmental data

(start_time <- Sys.time())
annotated_data <- parLapply(cl, 1:nrow(data), function(rown){
  print(paste0("Extracting environmental data for observation ", rown, " of ", nrow(data), "."))
  
  progress <- paste0("Extracting environmental data for observation ", rown, " of ", nrow(data), ".")
  write.table(progress, "C:/Users/hbronnvik/Documents/storkSSFs/output.txt", append = T, quote = F, col.names = F)
  
  d_year <- year(data[rown, "timestamp"])
  stamp <- round_date(data[rown, "timestamp"], unit = "hour")
  
  ## retrieve data with ECMWF variables (downloaded with the script "ecmwf_data.R")
  # the height of the boundary layer
  file <- grep(d_year, file_list_s[grepl("boundary_layer", file_list_s)], value = T)
  blh <- try(terra::rast(file))
  blh <- blh[[which(time(blh) == round_date(data[rown, "timestamp"], unit = "hour"))]]
  
  # the boundary layer height at the point of interest
  poi_blh <- terra::extract(blh, vect(data[rown, ], geom = c("long", "lat")))[,2]
  
  # the height of the pressure levels at the point of interest
  files <- grep(d_year, file_list_p[grepl("Geopotential", file_list_p)], value = T)
  
  heights <- lapply(files, function(file){
    geo <- try(terra::rast(file))
    geo <- geo[[which(time(geo) == stamp)]]
    poi_geoH <- terra::extract(geo, vect(data[rown, ], geom = c("long", "lat")))[ ,2]/9.80665
    lev <- str_split(file, "_")[[1]][3]
    h <- data.frame(level = lev, height = poi_geoH)
  }) %>% reduce(rbind) %>% arrange(height) %>% as.data.frame() 
  
  # the levels above and below the PBL at the point of interest
  low_pressure <- heights %>%
    filter(height > poi_blh) %>%
    filter(height == min(height)) %>% 
    mutate(level = as.numeric(level))
  high_pressure <- heights %>%
    filter(height < poi_blh) %>%
    filter(height == max(height)) %>% 
    mutate(level = as.numeric(level))
  
  # the values at the levels above and below the PBL for each pressure variable
  files <- file_list_p[grep("v_component_of_wind|u_component_of_wind|Geopotential", file_list_p, invert = T)]
  
  at_pressure <- lapply(files[grepl(d_year, files)], function(file){
    ras <- try(terra::rast(file))
    ras <- ras[[which(time(ras) == stamp)]]
    val <- terra::extract(ras, vect(data[rown, ], geom = c("long", "lat")))[ ,2]
    v <- names(ras)
    variable <- ifelse(grepl("q_", v), "q", "temp")
    lev <- str_split(file, "_")[[1]][3]
    df <- data.frame(variable = variable, level = lev, value = val)
  }) %>% reduce(rbind) 
  
  # interpolate to the BLH
  pressures <- data.frame(index = 1)
  pressures$blh_humidity <- int(poi_blh, high_pressure$height, low_pressure$height,
                                at_pressure$value[at_pressure$variable == "q" & at_pressure$level == high_pressure$level],
                                at_pressure$value[at_pressure$variable == "q" & at_pressure$level == low_pressure$level])
  
  
  pressures$blh_temperature <- int(poi_blh, high_pressure$height, low_pressure$height,
                                   at_pressure$value[at_pressure$variable == "temp" & at_pressure$level == high_pressure$level],
                                   at_pressure$value[at_pressure$variable == "temp" & at_pressure$level == low_pressure$level])
  
  
  pressures$blh_pressure <- int(poi_blh, high_pressure$height, low_pressure$height,
                                high_pressure$level, low_pressure$level)
  
  # now get the single level data for the point of interest
  files <- grep(d_year, file_list_s[!grepl("boundary", file_list_s)], value = T)
  
  singles <- lapply(files, function(file){
    ras <- try(terra::rast(file))
    ras <- ras[[which(time(ras) == stamp)]]
    val <- terra::extract(ras, vect(data[rown, ], geom = c("long", "lat")))[ ,2]
    v <- names(ras)[1]
    variable <- ifelse(grepl("t2m", v), "temp2m", 
                       ifelse(grepl("sp_", v), "p0", 
                              ifelse(grepl("ishf", v), "sflux", "mflux")))
    df <- data.frame(variable = variable, value = val)
  }) %>% 
    reduce(rbind) %>% 
    pivot_wider(names_from = "variable", values_from = "value")
  
  df <- cbind(data[1, 2:ncol(data)], singles, pressures[, 2:4])
  df$blh <- poi_blh
  df
}) %>% reduce(rbind)
Sys.time() - start_time # 20.59281 mins for 100 rows
# 292 rows in an hour

stopCluster(cl)

# estimate w*
annotated_data$w_star <- convective(data = annotated_data)

# Parallel processed, uses the layer below the PBL rather than interpolating
# -------------------------------------------------------

file_list_s <- list.files("C:/Users/hbronnvik/Documents/storkSSFs/ecmwf_data/single/", full.names = T) 
file_list_p <- list.files("C:/Users/hbronnvik/Documents/storkSSFs/ecmwf_data/pressure/", full.names = T) 

cl <- makeCluster(16)
clusterExport(cl, c("file_list_s", "file_list_p", "int", "data"))
clusterEvalQ(cl, c(library(tidyverse), library(terra)))

# loop through each individual point and annotate it with its environmental data

(start_time <- Sys.time())
annotated_data <- parLapply(cl, 1:100, function(rown){
  print(paste0("Extracting environmental data for observation ", rown, " of ", nrow(data), "."))
  
  progress <- paste0("Extracting environmental data for observation ", rown, " of ", nrow(data), ".")
  write.table(progress, "C:/Users/hbronnvik/Documents/storkSSFs/output2.txt", append = T, quote = F, col.names = F)
  
  d_year <- year(data[rown, "timestamp"])
  stamp <- round_date(data[rown, "timestamp"], unit = "hour")
  
  ## retrieve data with ECMWF variables (downloaded with the script "ecmwf_data.R")
  # the height of the boundary layer
  file <- grep(d_year, file_list_s[grepl("boundary_layer", file_list_s)], value = T)
  blh <- try(terra::rast(file))
  blh <- blh[[which(time(blh) == round_date(data[rown, "timestamp"], unit = "hour"))]]
  
  # the boundary layer height at the point of interest
  poi_blh <- terra::extract(blh, vect(data[rown, ], geom = c("long", "lat")))[,2]
  
  # the height of the pressure levels at the point of interest
  files <- grep(d_year, file_list_p[grepl("Geopotential", file_list_p)], value = T)
  
  heights <- lapply(files, function(file){
    geo <- try(terra::rast(file))
    geo <- geo[[which(time(geo) == stamp)]]
    poi_geoH <- terra::extract(geo, vect(data[rown, ], geom = c("long", "lat")))[ ,2]/9.80665
    lev <- str_split(file, "_")[[1]][3]
    h <- data.frame(level = lev, height = poi_geoH)
  }) %>% reduce(rbind) %>% arrange(height) %>% as.data.frame() 
  
  # the level below the PBL at the point of interest
  low_pressure <- heights %>%
    filter(height > poi_blh) %>%
    filter(height == min(height)) %>% 
    mutate(level = as.numeric(level))

  # the value at the level below the PBL for each pressure variable
  files <- file_list_p[grep("v_component_of_wind|u_component_of_wind|Geopotential", file_list_p, invert = T)]
  files <- grep(low_pressure$level, files[grepl(d_year, files)], value = T)
  
  pressures <- lapply(files, function(file){
    ras <- try(terra::rast(file))
    ras <- ras[[which(time(ras) == stamp)]]
    val <- terra::extract(ras, vect(data[rown, ], geom = c("long", "lat")))[ ,2]
    v <- names(ras)
    variable <- ifelse(grepl("q_", v), "q", "temp")
    lev <- str_split(file, "_")[[1]][3]
    df <- data.frame(variable = variable, level = lev, value = val)
  }) %>% 
    reduce(rbind) %>% 
    pivot_wider(names_from = "variable", values_from = "value")
  
  # now get the single level data for the point of interest
  files <- grep(d_year, file_list_s[!grepl("boundary", file_list_s)], value = T)
  
  singles <- lapply(files, function(file){
    ras <- try(terra::rast(file))
    ras <- ras[[which(time(ras) == stamp)]]
    val <- terra::extract(ras, vect(data[rown, ], geom = c("long", "lat")))[ ,2]
    v <- names(ras)[1]
    variable <- ifelse(grepl("t2m", v), "temp2m", 
                       ifelse(grepl("sp_", v), "p0", 
                              ifelse(grepl("ishf", v), "sflux", "mflux")))
    df <- data.frame(variable = variable, value = val)
  }) %>% 
    reduce(rbind) %>% 
    pivot_wider(names_from = "variable", values_from = "value")
  
  df <- cbind(data[1, 2:ncol(data)], singles, pressures)
  df$blh <- poi_blh
  df
}) %>% reduce(rbind)
Sys.time() - start_time # 8.181818 mins for 100 rows

stopCluster(cl)

# estimate w*
annotated_data$w_star <- convective(data = annotated_data)

# Parallel processed, uses the layer below the PBL rather than interpolating, extracts multiple points at once
# -------------------------------------------------------

file_list_s <- list.files("C:/Users/hbronnvik/Documents/storkSSFs/ecmwf_data/single/", full.names = T) 
file_list_p <- list.files("C:/Users/hbronnvik/Documents/storkSSFs/ecmwf_data/pressure/", full.names = T) 

cl <- makeCluster(16)
clusterExport(cl, c("file_list_s", "file_list_p", "int", "data"))
clusterEvalQ(cl, c(library(tidyverse), library(terra)))

# loop through each individual point and annotate it with its environmental data

(start_time <- Sys.time())
annotated_data <- parLapply(cl, split(data, data$timestamp), function(steps){
  # print(paste0("Extracting environmental data for timestamp ", which(unique(data$timestamp) == steps$timestamp[1]), " of ", length(unique(data$timestamp)), "."))
  
  progress <- paste0("Extracting environmental data for stratum ", which(unique(data$timestamp) == steps$timestamp[1]), " of ", length(unique(data$timestamp)), ".")
  write.table(progress, "C:/Users/hbronnvik/Documents/storkSSFs/output3.txt", append = T, quote = F, col.names = F)
  
  d_year <- unique(year(steps$timestamp))
  stamp <- round_date(unique(steps$timestamp), unit = "hour")
  
  ## retrieve data with ECMWF variables (downloaded with the script "ecmwf_data.R")
  # the height of the boundary layer
  file <- grep(d_year, file_list_s[grepl("boundary_layer", file_list_s)], value = T)
  blh <- try(terra::rast(file))
  blh <- blh[[which(time(blh) == stamp)]]
  
  # the boundary layer height at the point of interest
  poi_blh <- terra::extract(blh, vect(steps, geom = c("long", "lat")))[,2]
  
  # the height of the pressure levels at the point of interest
  files <- grep(d_year, file_list_p[grepl("Geopotential", file_list_p)], value = T)
  
  heights <- lapply(files, function(file){
    geo <- try(terra::rast(file))
    geo <- geo[[which(time(geo) == stamp)]]
    poi_geoH <- terra::extract(geo, vect(steps, geom = c("long", "lat")))[ ,2]/9.80665
    lev <- str_split(file, "_")[[1]][3]
    h <- data.frame(poi_geoH)
    colnames(h) <- lev
    h
  }) %>% 
    reduce(cbind) %>% 
    # mutate(blh = poi_blh) %>% 
    as.data.frame() 
  
  # take the lowest pressure level below the PBL for each point.
  # in the event of the BLH being below the lowest pressure level, take the lowest pressure level.
  low_pressures <- lapply(1:nrow(heights), function(x){
    v <- heights[x,] %>% 
      t() %>% 
      as.data.frame() %>% 
      rownames_to_column() %>% 
      rename(h = 2,
             level = rowname) %>%
      arrange(h)
    if(min(v$h) < poi_blh[x]){
      v <- v %>% 
        filter(h < poi_blh[x]) %>% 
        filter(h == max(h))
    }else{
      v <- v %>% 
        filter(h == min(h))
    }
    v
  }) %>% reduce(rbind)
  
  steps$level <- low_pressures$level
  
  # now we have the pressure levels we need to extract values from, name the files
  files <- file_list_p[grep("Specific humidity|Temperature", file_list_p)]
  file_levels <- paste(unique(low_pressures$level), collapse = "|")
  files <- grep(paste0(file_levels), files[grepl(d_year, files)], value = T)
  variables <- c("Specific humidity", "Temperature")

  # the value at the level below the PBL for each pressure variable
  pressures <- lapply(split(steps, steps$level), function(l){
    temp <- lapply(variables, function(var){
      file <- grep(var, grep(unique(l$level), files, value = T), value = T)
      ras <- try(terra::rast(file))
      ras <- ras[[which(time(ras) == stamp)]]
      val <- terra::extract(ras, vect(l, geom = c("long", "lat")))[ ,2]
      v <- names(ras)
      variable <- ifelse(grepl("q_", v), "q", "temp")
      lev <- str_split(file, "_")[[1]][3]
      df <- data.frame(level = lev, value = val, long = l$long, lat = l$lat, timestamp = l$timestamp)
      colnames(df)[which(colnames(df) == "value")] <- variable
      return(df)
    }) %>% 
      reduce(cbind)
  }) %>% reduce(rbind)
  
  # now get the ground level data for the point of interest
  files <- grep(d_year, file_list_s[!grepl("boundary", file_list_s)], value = T)
  
  singles <- lapply(files, function(file){
    ras <- try(terra::rast(file))
    ras <- ras[[which(time(ras) == stamp)]]
    val <- terra::extract(ras, vect(steps, geom = c("long", "lat")))[ ,2]
    v <- names(ras)[1]
    variable <- ifelse(grepl("t2m", v), "temp2m", 
                       ifelse(grepl("sp_", v), "p0", 
                              ifelse(grepl("ishf", v), "sflux", "mflux")))
    df <- data.frame(value = val, long = steps$long, lat = steps$lat, timestamp = steps$timestamp)
    colnames(df)[which(colnames(df) == "value")] <- variable
    return(df)
  }) %>% reduce(cbind)
  
  extracts <- full_join(singles[c(1:4, 5, 9, 13)], pressures[, 2:7])
  df <- full_join(steps, extracts)
  df$blh <- poi_blh
  df
}) %>% reduce(rbind)
Sys.time() - start_time # 4.847536 secs for 100 rows

stopCluster(cl)

# a function to calculate w*
convective <- function(data) {
  
  df <- data # the data input
  z <- df[, "blh"] # the boundary layer height
  s_flux <- df[, "sflux"] # the surface sensible heat flux
  m_flux <- df[, "mflux"] # the moisture flux
  T_k_2m <- df[, "temp2m"] # the temperature at ground level
  p0 <- df[, "p0"] # the pressure at ground level (Pascals)
  p0 <- p0/100  # the pressure at ground level (mbar)
  q <- df[, "q"] # the humidity below the boundary layer height
  p1 <- as.numeric(df[, "level"]) # the pressure below the boundary layer height (mbar)
  T_K_blh <- df[, "temp"] # the temperature below the boundary layer height
  
  k <- 0.2854 # Poisson constant
  g <- 9.80665 # acceleration due to gravity
  T_c_2m <- T_k_2m - 273.15 # surface temperature in degrees Celsius
  c_p <- 1012 # the isobaric mass heat capacity of air at common conditions
  p <- 1.225 # the density of air at sea level and 273 K
  
  Theta_k_z <- T_K_blh*((p0/p1)^k)
  Thetav_k_z <- Theta_k_z*(1+(0.61*q))
  wT <- (s_flux * -1)/c_p/p # reverse the sign. ECMWF upward fluxes are negative
  wq <- (m_flux * -1)/p # reverse the sign
  
  wthetav <- wT + 0.61 * T_c_2m * wq
  
  w_star <- CubeRoot(g*z*(wthetav/Thetav_k_z))
  
  return(w_star)
  
}

# estimate w*
annotated_data$w_star <- convective(data = annotated_data)



library(readr)
library(tidyverse)
library(ggplot2)
library(parallel)
library(terra)
library(rasterVis)

# a function to take the cube root
CubeRoot<-function(x){
  sign(x)*abs(x)^(1/3)
}

# the location data
data <- readRDS("C:/Users/hbronnvik/Documents/storkSSFs/used_av_df_230612.rds") %>% 
  rename(long = "location.long",
         lat = "location.lat") %>% 
  group_by(stratum) %>% 
  slice(1:100) %>%
  ungroup() %>% 
  # filter(year(timestamp) == 2013) %>% 
  as.data.frame()

# Parallel processed, uses the layer below the PBL rather than interpolating, extracts a year at once
# -------------------------------------------------------

file_list_s <- list.files("C:/Users/hbronnvik/Documents/storkSSFs/ecmwf_data/single/", full.names = T) 
file_list_p <- list.files("C:/Users/hbronnvik/Documents/storkSSFs/ecmwf_data/pressure/", full.names = T) 

steps <- data %>% filter(year(timestamp) %in% c(2018, 2019, 2020, 2021, 2022, 2023))
data_ls <- split(steps, year(steps$timestamp))

cl <- makeCluster(12)
clusterExport(cl, c("file_list_s", "file_list_p", "data_ls"))
clusterEvalQ(cl, c(library(tidyverse), library(terra)))

# loop through each year and annotate it with its environmental data
(start_time <- Sys.time())
annotated_data <- parLapply(cl, data_ls, function(steps){
  progress <- paste0("Extracting environmental data for ", unique(year(steps$timestamp[1])), ".")
  write.table(progress, "C:/Users/hbronnvik/Documents/storkSSFs/output3.txt", append = T, quote = F, col.names = F)
  
  d_year <- unique(year(steps$timestamp))
  stamp <- round_date(unique(steps$timestamp), unit = "hour")
  
  ## retrieve data with ECMWF variables (downloaded with the script "ecmwf_data.R")
  # the height of the boundary layer
  file <- grep(d_year, file_list_s[grepl("boundary_layer", file_list_s)], value = T)
  blh <- try(terra::rast(file))
  blh <- blh[[which(time(blh) %in% stamp)]]
  
  # create a list to go through each timestamp individually
  steps_ls <- split(steps, steps$timestamp)
  
  # the boundary layer height at the points of interest
  # start_time <- Sys.time()
  pois_blh <- lapply(steps_ls, function(df){
    ts <- round_date(unique(df$timestamp), unit = "hour")
    h <- blh[[which(time(blh) == ts)]]
    poi_blh <- terra::extract(h, vect(df, geom = c("long", "lat")))[,2]
  }) %>% unlist()
  # Sys.time()-start_time
  
  steps$blh <- as.numeric(pois_blh)
  
  # the height of the pressure levels at the points of interest
  files <- grep(d_year, file_list_p[grepl("Geopotential", file_list_p)], value = T)
  
  # start_time <- Sys.time()
  heights <- lapply(files, function(file){
    geo <- try(terra::rast(file))
    temp <- lapply(steps_ls, function(df){
      ts <- round_date(unique(df$timestamp), unit = "hour")
      geo <- geo[[which(time(geo) == ts)]]
      poi_geoH <- terra::extract(geo, vect(df, geom = c("long", "lat")))[ ,2]/9.80665
      lev <- str_split(file, "_")[[1]][3]
      h <- data.frame(poi_geoH)
      colnames(h) <- lev
      h <- h %>% 
        mutate(long = df$long,
               lat = df$lat)
      return(h)
    }) %>% reduce(rbind)
  }) %>% 
    reduce(cbind) %>% 
    # mutate(blh = poi_blh) %>% 
    as.data.frame() 
  # Sys.time()-start_time #Time difference of 26.82361 mins
  
  heights <- heights %>% 
    select(colnames(heights)[which(!colnames(heights) %in% c("long", "lat"))]) %>% 
    mutate(long = heights[ , 2],
           lat = heights[ , 3])
  
  steps <- full_join(steps, heights)
  
  # take the lowest pressure level below the PBL for each point.
  # in the event of the BLH being below the lowest pressure level, take the lowest pressure level.
  # start_time <- Sys.time()
  low_pressures <- lapply(1:nrow(steps), function(x){
    v <- steps[x, c(15:28)] %>% 
      t() %>% 
      as.data.frame() %>% 
      rownames_to_column(var = "level") %>% 
      rename(h = 2) %>%
      arrange(h)
    if(min(v$h) < steps$blh[x]){
      v <- v %>% 
        filter(h < steps$blh[x]) %>% 
        filter(h == max(h))
    }else{
      v <- v %>% 
        filter(h == min(h))
    }
    v <- v %>% 
      mutate(long = steps$long[x],
             lat = steps$lat[x])
  }) %>% reduce(rbind)
  # Sys.time()-start_time # Time difference of 2.84052 mins
  
  steps <- full_join(steps, low_pressures)
  
  # now we have the pressure levels we need to extract values from, name the files
  files <- file_list_p[grep("Specific humidity|Temperature", file_list_p)]
  file_levels <- paste(unique(low_pressures$level), collapse = "|")
  files <- grep(paste0(file_levels), files[grepl(d_year, files)], value = T)
  variables <- c("Specific humidity", "Temperature")
  
  # the value at the level below the PBL for each pressure variable
  # start_time <- Sys.time()
  # split the year into different pressure levels so that only one file need be called per group
  pressures <- lapply(split(steps, steps$level), function(l){
    # for each variable, 
    temp <- lapply(variables, function(var){
      # call the file at the level
      file <- grep(var, grep(unique(l$level), files, value = T), value = T)
      # subset it to have only the times when the PBL was at that pressure level
      ras <- try(terra::rast(file))
      # then go through each time step individually and extract the values at each location
      val <- lapply(split(l, l$timestamp), function(lt){
        ts <- round_date(unique(lt$timestamp), unit = "hour")
        r <- ras[[which(time(ras) == ts)]]
        val <- terra::extract(r, vect(lt, geom = c("long", "lat")))[ ,2] 
        val_df <- data.frame(timestamp = lt$timestamp, value = val, long = lt$long, lat = lt$lat)
      }) %>% reduce(rbind)
      # tidy and save
      v <- names(ras)[1]
      variable <- ifelse(grepl("q_", v), "q", "temp")
      lev <- str_split(file, "_")[[1]][3]
      df <- val %>% 
        mutate(level = lev)
      colnames(df)[which(colnames(df) == "value")] <- variable
      return(df)
    }) %>% 
      reduce(cbind)
  }) %>% reduce(rbind)
  # Sys.time()-start_time # Time difference of 14.89281 mins
  
  # now get the ground level data for the point of interest
  files <- grep(d_year, file_list_s[!grepl("boundary", file_list_s)], value = T)
  variables <- c("instantaneous_moisture_flux", "temperature", "instantaneous_surface_sensible", "surface_pressure")
  
  # start_time <- Sys.time()
  # split the year into different pressure levels so that only one file need be called per group
  singles <- lapply(variables, function(var){
    # call the file at the level
    file <- grep(var, files, value = T)
    # subset it to have only the times when the PBL was at that pressure level
    ras <- try(terra::rast(file))
    # then go through each time step individually and extract the values at each location
    val <- lapply(split(steps, steps$timestamp), function(lt){
      ts <- round_date(unique(lt$timestamp), unit = "hour")
      r <- ras[[which(time(ras) == ts)]]
      val <- terra::extract(r, vect(lt, geom = c("long", "lat")))[ ,2] 
      val_df <- data.frame(timestamp = lt$timestamp, value = val, long = lt$long, lat = lt$lat)
    }) %>% reduce(rbind)
    # tidy and save
    v <- names(ras)[1]
    variable <- ifelse(grepl("t2m_", v), "T2m", 
                       ifelse(grepl("ie_", v), "mflux", 
                              ifelse(grepl("ishf_", v), "sflux",
                                     ifelse(grepl("sp_", v), "p0", NA))))
    colnames(val)[2] <- variable
    return(val)
  }) %>% 
    reduce(cbind)
  # Sys.time()-start_time # 7.434673 mins
  
  singles <- singles %>% 
    select(colnames(singles)[which(!colnames(singles) %in% c("long", "lat", "timestamp"))]) %>% 
    mutate(long = singles[ , 3],
           lat = singles[ , 4], 
           timestamp = singles[ , 5])
  
  extracts <- full_join(singles, pressures[, c(1:5, 7)])
  df <- full_join(steps, extracts)
  df
}) %>% reduce(rbind)
Sys.time() - start_time 
# Time difference of 59.21467 mins for year 2013 (44981 rows)
# Time difference of 5.136102 hours for 2014
# Time difference of 9.617492 hours for 2015 & 2016 (731200 rows)
# Time difference of 5.531022 hours for 2017 (272900 rows)

stopCluster(cl)

annotated_data <- annotated_data %>% 
  as.data.frame()

# a function to calculate w*
convective <- function(data) {
  
  df <- data # the data input
  z <- df$blh # the boundary layer height
  s_flux <- df$sflux # the surface sensible heat flux
  m_flux <- df$mflux # the moisture flux
  T_k_2m <- df$T2m # the temperature at ground level
  p0 <- df$p0 # the pressure at ground level (Pascals)
  p0 <- p0/100  # the pressure at ground level (mbar)
  q <- df$q # the humidity below the boundary layer height
  p1 <- as.numeric(df$level) # the pressure below the boundary layer height (mbar)
  T_K_blh <- df$temp # the temperature below the boundary layer height
  
  k <- 0.2854 # Poisson constant
  g <- 9.80665 # acceleration due to gravity
  T_c_2m <- T_k_2m - 273.15 # surface temperature in degrees Celsius
  c_p <- 1012 # the isobaric mass heat capacity of air at common conditions
  p <- 1.225 # the density of air at sea level and 273 K
  
  Theta_k_z <- T_K_blh*((p0/p1)^k)
  Thetav_k_z <- Theta_k_z*(1+(0.61*q))
  wT <- (s_flux * -1)/c_p/p # reverse the sign. ECMWF upward fluxes are negative
  wq <- (m_flux * -1)/p # reverse the sign
  
  wthetav <- wT + 0.61 * T_c_2m * wq
  
  w_star <- CubeRoot(g*z*(wthetav/Thetav_k_z))
  
  return(w_star)
  
}

# estimate w*
annotated_data$w_star <- convective(data = annotated_data)

ggplot(annotated_data, aes(w_star)) +
  geom_histogram(color = "#970C11", fill = "firebrick") +
  labs(x = "w*", y = "Count") +
  theme_classic()

# saveRDS(annotated_data, file = "C:/Users/hbronnvik/Documents/storkSSFs/ecmwf_data/annotated/adata_2018_19_20.rds")

cor.test(annotated_data$w_star, annotated_data$blh)
ggplot(annotated_data, aes(w_star, blh)) +
  geom_point(alpha = 0.01) +
  geom_smooth(method = "gam", color = "firebrick") +
  labs(y = "Boundary layer height (m)", x = "w* (m/s)") +
  geom_vline(xintercept = 0, lty = 2) +
  theme_classic()

