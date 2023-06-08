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
convective <- function(blh, T2m, Tblh, s_flux, m_flux, qblh, p0, p1) {
  
  p <- 0.2854 # Poisson constant.
  g = 9.80665 # acceleration due to gravity.
  z <- blh # boundary layer height.
  T_k_2m <- T2m # surface temperature.
  T_K_blh <- Tblh # temperature at the boundary layer height.
  q <- qblh # humidity at the boundary layer height.
  T_c_2m <- T_k_2m - 273.15 # surface temperature in degrees Celsius.
  Theta_k_z <- T_K_blh*((p0/p1)^p)
  Thetav_k_z <- Theta_k_z*(1+(0.61*q))
  wT <- (s_flux * -1)/1012/1.2 #reverse the sign. ECMWF upward fluxes are negative
  wq <- (m_flux * -1)/1.2 #reverse the sign.
  
  wthetav <- wT + 0.61 * T_c_2m * wq
  
  w_star <- CubeRoot(g*z*(wthetav/Thetav_k_z))
  
  return(w_star)
  
}




# the location data
data <- read_csv("C:/Users/hbronnvik/Documents/storkSSFs/burst_data/speed40x80/A_2023-04-13.csv") %>% 
  rename(long = "location-long",
         lat = "location-lat") %>% 
  slice(1:100) %>% 
  as.data.frame()


# loop through each individual point and annotate it with its environmental data
(start_time <- Sys.time())
annotated_data <- lapply(1:nrow(data), function(rown){
  
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
Sys.time() - start_time

annotated_data$w_star <- convective(annotated_data$blh, annotated_data$temp2m, annotated_data$blh_temperature,
                                    annotated_data$hflux, annotated_data$mflux, annotated_data$blh_humidity,
                                    annotated_data$p0, annotated_data$blh_pressure)


