### A function to segment white stork tracking data migration tracks for consecutive years
### Hester Bronnvik
### 2025
### hester.bronnvik@gmail.com

# library(data.table)
# library(dplyr)
# library(tibble)
# library(lubridate)

# The approach is to use ground speed thresholds to define rapid movement
# then, use time of the year to define season of the migration
# if the time is unclear, add in direction.

# This assumes an input of clean data split into a list of individuals

# required information
# high_d_thresh <- the number of meters required to define what is certainly migration
# low_d_thresh <- the number of meters required to define "not sedentary"
# g_thresh <- days between low_d_thresh days to define bursts of activity
# fall_months <- the months that are probably post-breeding season migrations (e.g. c(8:11))
# spring_months <- the months that are probably pre-breeding season migrations (e.g. c(1:6))
# distance_col <- the name of the column containing distance moved by the animal per day

# optional information
# visDODs is a data frame with individual ID and death date
# if use_ragged is true, days with speeds between low_d_thresh and high_d_thresh are kept
# if the bird died then (i.e. the reason the distance was not achieved is that the bird died)
# and if the bird also had traveled high_d_thresh in the preceding g_thresh days
# if use_ragged is false, low_d_thresh is used to define tracks but only high_d_thresh is maintained

seg_storks <- function(clean_locations, high_d_thresh, low_d_thresh, g_thresh,
         fall_months, spring_months, distance_col, use_ragged = F){
  # ensure required packages are loaded
  req_packs <- c("dplyr", "data.table", "lubridate", "tibble")
  invisible(lapply(req_packs, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste0("Package ", pkg, "' is required but not installed."))
    }
  }))
  # symbolize the column containing daily ground speeds
  dist_sym <- rlang::sym(distance_col)
  
  # isolate long-distance movement and label it
  ld_locs <- lapply(1:length(clean_locations), function(i){
    # one animal
    x <- clean_locations[[i]]
    # update progress
    print(paste0("Segmenting data for individual ", x$individual.id[1], ". Animal ",
                 i, " of ", length(clean_locations), "."), quote = F)
    df <- x %>% 
      dplyr::select(individual.id, timestamp, location.lat, location.long, !!dist_sym, daily_direction, ground_speed_15) %>% 
      # remove all days that the animal did not travel the pre-set distance
      filter(!!dist_sym > low_d_thresh) %>%
      arrange(timestamp) %>%
      mutate(# add a column to measure difference in time between consecutive points
        time_lag = as.numeric(timestamp - lag(timestamp), units = "secs"),
        # insert an enormous time lag for the first location to remove the NA
        time_lag = ifelse(is.na(time_lag), 1e5, time_lag),
        # add whether there was a gap between days of pre-set distance travel
        new_burst = ifelse(round(time_lag) <= days(g_thresh), F, T),
        # newCluster = ifelse(round(timeLag) <= weeks(6), F, T),
        # take the cumulative sum to act as a unique ID for each burst
        cumu_check_for_burst = cumsum(new_burst),
        # cumu_check_for_clust = cumsum(newCluster)
      ) %>%
      group_by(cumu_check_for_burst) %>% 
      # for the bursts, count number of days
      mutate(burst_days = length(unique(date(timestamp))),
             # add an ID to each burst
             burstID = as.character(cur_group_id())) %>%
      # remove the bursts that do not meet user-set criteria
      # filter(burstLength > MinBurstLength) %>% 
      ungroup() %>%
      # clean up the sorting columns
      dplyr::select(-"new_burst", -"cumu_check_for_burst")
    
    # if the animal met the pre-set distance criterion
    if(nrow(df) > 0){
      # use the above bursts to classify locations as belonging to spring or fall
      class_df <- df %>%
        group_by(burstID) %>% 
        # if the group of long-distance movements end at a lower latitude than they started,
        # this is a southward burst
        mutate(burst_angle = ifelse(location.lat[1] > location.lat[n()], "south", "north"),
               # if the group of long-distance movements occurred in months defined neither as fall 
               # nor as spring, then the season is inconclusive
               burst_season = ifelse(month(timestamp[n()]) %in% fall_months, "late", 
                                     ifelse(month(timestamp[n()]) %in% spring_months, "early", "change")),
               # for those inconclusive locations, use the burst angle to determine whether this was 
               # a movement towards the breeding grounds or away
               burst_class = ifelse(burst_season == "early", "spring", 
                                    ifelse(burst_season == "late", "fall",
                                           ifelse(burst_angle == "north" & burst_season == "change", "spring", 
                                                  ifelse(burst_angle == "south" & burst_season == "change", "fall", NA)))),
               # ensure that any bursts in december count as belonging to the next spring
               burst_year = ifelse(month(timestamp[n()]) == 12 & burst_class == "spring", year(years(1) + timestamp[n()]), year(timestamp[n()])),
               # and finally, paste together the ID, season, and year so that all bursts
               # matching those criteria are grouped under a single unique identifier
               trackID = paste(unique(individual.id), burst_class, burst_year, sep = "_")) %>% 
        ungroup()
      
      # a diagnostic plot
      # ggplot(class_df, aes(location.long, location.lat, color = trackID)) +
      #   borders(xlim = c(-20, 20), ylim = c(0, 50)) +
      #   geom_point()
      
      # if interested in using mortality to discard some movements,
      # remove ragged edges of tracks the bird survived
      if(use_ragged){
        # identify ragged edges as bursts that do not have high speed days in them
        ragged <- class_df %>% 
          group_by(burstID) %>% 
          mutate(high_speed = max(!!dist_sym > high_d_thresh)) %>%
          ungroup() %>% 
          filter(high_speed == F)
        # did the bird die?
        if(unique(x$individual.id) %in% visDODs$individual_id){
          dod <- visDODs %>% 
            filter(individual_id == unique(x$individual.id)) %>% 
            dplyr::select(dod) %>% 
            deframe()
        }else{
          dod <- NA
        }
        # using that death date, and whether the bird had travel days
        # between the low and high threshold speeds:
        if(nrow(ragged) > 0){
          if(!is.na(dod)){
            # did the bird die in a ragged edge?
            cutoff <- date(dod) %in% unique(date(ragged$timestamp))
            # if it did, 
            if(cutoff == T){
              # which burstID?
              death_burst <- class_df %>% 
                filter(date(timestamp) == date(dod)) %>% 
                dplyr::select(burstID) %>% 
                unique() %>% 
                deframe() %>% 
                as.numeric()
              # was the previous burst high speed?
              ifelse((death_burst-1) %in% as.numeric(unique(ragged$burstID)), 
                     "Ragged", "Rapid")
              # if the "burst" preceding the one with a death was high speed,
              if((death_burst-1) %in% as.numeric(unique(ragged$burstID))){
                # filter class_df to contain only high speed bursts and the one with a death
                class_df <- class_df %>% 
                  filter(!burstID %in% ragged$burstID | burstID == death_burst)
              }
              # if the bird did not die on a mid-speed day, keep only high speed bursts
            }else{
              class_df <- class_df %>% 
                filter(!burstID %in% ragged$burstID)
            }
            # if the bird did transmit a mid-speed day, keep only high speed bursts
          }else{
            class_df <- class_df %>% 
              filter(!burstID %in% ragged$burstID)
          }
        }
        # if the user is ignoring mortality, keep only high speed bursts
      }else{
        class_df <- class_df %>% 
          group_by(burstID) %>% 
          mutate(high_speed = max(!!dist_sym > high_d_thresh)) %>%
          ungroup() %>% 
          filter(high_speed == T)
      }
      
      # whether or not the mortality dates were used, 
      # add the classified tracks back on to the full data
      if(nrow(class_df) > 0){
        # just the first and last locations of each track
        tracks <- class_df %>%
          group_by(trackID) %>% 
          slice(1, n()) %>% 
          group_split()
        # for each track, define all days between start and end as part of the track, 
        # regardless of daily distance
        tracks <- lapply(tracks, function(track){
          ID <- unique(track$trackID)
          track <- x %>% 
            filter(between(timestamp, track$timestamp[1], track$timestamp[2])) %>% 
            mutate(trackID = ID)
          track
        }) %>% rbindlist()
        
        suppressMessages(x <- x %>% left_join(tracks))
        
        # another diagnostic plot
        # ggplot(x, aes(location.long, location.lat, color = trackID)) +
        #   borders(xlim = c(-20, 10), ylim = c(40, 50)) +
        #   geom_point()
        # 
        # x %>% 
        #   filter(!is.na(trackID)) %>% 
        #   sf::st_as_sf(coords = c("location.long", "location.lat"), crs = "EPSG:4326") %>% 
        #   mapview::mapview(zcol = "trackID")
        
        return(x)
      }
    }
  })
}

# # an example usage
# clocs <- readRDS("/home/hbronnvik/Documents/storkSSFs/clean_locations_2023-08-30.rds")
# build <- seg_storks(clocs, 70000, 40000, 7, c(8:11), c(1:6), "daily_dist")
# # remove empty items (birds that do not migrate)
# build <- build[lapply(build, length) > 1]

