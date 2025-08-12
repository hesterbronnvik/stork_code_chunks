### Get environmental data from the Copernicus Climate Data Store
library(ecmwfr)
library(lubridate)
cds.key <- "0216b1a7-50e7-4799-8247-a2e644cfe518"
wf_set_key(user = "4c3ff79f-8b52-4928-8e25-99a90adf0a88", key = cds.key) # add user here


# define the pressure levels at which to download data (if needed)
# levels <- c(seq(from = 600, to = 1000, by = 25))
# levels <- levels[!levels %in% c(625, 675, 725)]

# define the variables to download
variables <- c("boundary_layer_height")

#list the years with data
years <- as.character(c(2015:2025))
# loop through each year 
lapply(years, function(z){
  # set up the request to the archive
  request <- list(
    # the repository
    "dataset_short_name" = "reanalysis-era5-single-levels",
    "product_type"   = "reanalysis",
    "variable"       = "boundary_layer_height",
    "year"           = z,
    # all months
    "month"          = c(paste0(0:12)),
    # all days
    "day"            = c(paste0(1:31)),
    # all hours
    "time"           = c("14:00"),
    # the top left and bottom right points in lat/long
    "area"           = "15/-90/-55/-30",
    "format"         = "netcdf",
    # the file name
    "target"         = paste0("boundary_layer_height_", z, ".nc")
  )
  
  file <- wf_request(user = "4c3ff79f-8b52-4928-8e25-99a90adf0a88",
                     request = request,
                     transfer = TRUE,
                     path = "/home/hbronnvik/Documents/chapter3/maguari_storks/PBLH/",
                     verbose = TRUE)
})
