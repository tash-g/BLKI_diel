## ---------------------------
##
## Script name: 2_prepare_attendance_data.R
##
## Purpose of script: Quantify visits to the colony across time 
##
## Author: Dr. Natasha Gillies
##
## Created: 2022-10-03
##
## Email: gilliesne@gmail.com
##
## ---------------------------


# Preamble ---------------------------------------------------------------------

# Define the packages
packages <- c("zoo", "stringr", "lubridate", "purrr", "dplyr", "suncalc", "magrittr",
              "tidyr", "lutz")

# Install packages not yet installed - change lib to library path
# installed_packages <- packages %in% rownames(installed.packages())

# if (any(installed_packages == FALSE)) {
#  install.packages(packages[!installed_packages])
# }

# Load packages and functions
invisible(lapply(packages, library, character.only = TRUE))

source("BLKI_functions.R")


# Load and process data --------------------------------------------------------

# Load the full trip data
load("Data_inputs/BLKI_gls-trips.RData")
gls_trips$year <- substr(gls_trips$start, 1, 4)

# Get data ready for processing
gls_trips <- rename(gls_trips, start_time = start, end_time = end)
gls_trips$start <- as.POSIXct(gls_trips$start_time, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
gls_trips$end_time <- as.POSIXct(gls_trips$end_time, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# Load the hourly behaviour data
load("Data_inputs/BLKI_gls-hourly-behaviour.RData")
gls_hourly$year <- substr(gls_hourly$date, 1, 4)

# Load temperature file
load("Data_inputs/temperature_data_2013-2022.RData")
col_temp_data %<>% dplyr::select(colony, datetime, temp)

# Process year by year, colony by colony, due to very large datasets -----------
years <- unique(gls_trips$year)
years <- years[order(years)]
kit.list <- vector(mode = "list", length = length(years))

for (i in 1:length(years)) {
  
  pb <- txtProgressBar(min = 0, max = length(years), style = 3)
  setTxtProgressBar(pb, i)
  
  trips_year <- subset(gls_trips, year == years[i])
  
  colonies <- unique(trips_year$colony)[order(unique(trips_year$colony))]
  col.list <- vector(mode = "list", length = length(colonies))
  
  for (j in 1:length(colonies)) {
  
  trips_col <- subset(trips_year, colony == colonies[j])
  behaviour_col <- subset(gls_hourly, year == years[i] & colony == colonies[j])
  
  ### Find number minutes attendance at colony by hour -------------------------
  
  kits.att <- trips_col %>% 
    group_by(ring, loc, colony, col_lat, col_lon) %>%
    mutate(time = purrr::map2(start_time, end_time, seq, by = 'min')) %>%
    tidyr::unnest(time) %>%
    mutate(hour = hour(time), date = as.Date(time)) %>%
    dplyr::count(date, hour) %>% 
    mutate(datetime = paste(date, 
                            format(strptime(hour, format = "%H"), "%H:%M:%S"))) %>%
    data.frame()
  
  kits.att <- kits.att[order(kits.att$ring, kits.att$datetime),]
  
  # Set trip attendance values to inverse
  kits.att$n[kits.att$loc == "trip"] <- 60 - kits.att$n[kits.att$loc == "trip"]
  
  # Where duplicates, remove the trip value
  kits.att <-  kits.att[order(kits.att$ring, kits.att$datetime),]
  kits.att <- kits.att[!duplicated(kits.att[c("ring","datetime")]),]
  
  # Convert date and hour column to posixct
  kits.att$datetime <- as.POSIXct(kits.att$datetime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  
  ## Add year and hour column
  kits.att$year <- substr(kits.att$date, 1, 4)
  kits.att$hour <- strftime(kits.att$datetime, format = "%H", tz = "UTC")
  
  
  ### Find number of visits per hour -------------------------------------------
  
  # Make another hour column
  trips_col$hour <- strftime(trips_col$start_time, format = "%H", tz = "UTC")
  
  # Count number of colony visits per hour
  kits.visits <- trips_col %>% 
    filter(loc != "trip") %>% 
    group_by(ring, date = as.Date(start_time), hour) %>%
    summarize(num_visits = n()) %>%
    data.frame()
  
  # Merge into kits.att dataframe
  kits.att <- merge(kits.att, kits.visits, by = c("ring", "date", "hour"), all.x = TRUE)
  kits.att <- kits.att[order(kits.att$ring, kits.att$datetime),]
  
  # Set NA num_visits to 0
  kits.att$num_visits[is.na(kits.att$num_visits)] <- 0
  
  # Remove loc column - no longer needed
  kits.att$loc <- NULL
  
  
  ### Extend hours column to include missing hours -----------------------------
  
  kits.att <- kits.att %>%
    group_by(ring) %>%
    tidyr::complete(datetime = seq(min(datetime), max(datetime), by = 'hour'), 
                    fill = list(n = 0, num_visits = 0)) %>%
    tidyr::fill(ring) %>%
    data.frame()
  
  # Give each timestamp an ID per bird; rename 'n' column 
  kits.att <- kits.att %>%
    group_by(ring) %>%
    mutate(time_ID = row_number()) %>%
    rename(att_time = n) %>%
    data.frame()
  
  # Add columns representing total possible attendance minutes
  kits.att$att_poss <- 60
  
  # Find proportional attendance at colony
  kits.att$prop_att <- kits.att$att_time/kits.att$att_poss
  
  # Ensure birds have at least 60 hours of data (lag time used for analysis)
  kits.att <- kits.att %>% 
    group_by(ring) %>% 
    filter(n_distinct(datetime) > 60)
  
  ### Add behaviour data --------------------------------------------------------
  
  behaviour_col %<>% dplyr::select(-c(n_forage, for_poss, year, colony))
  
  kits.att <- merge(kits.att, behaviour_col, 
                    by = c("ring", "date", "hour"), all.x = T)
  
  kits.att[,c("date", "hour")] <- NULL
  
  # Turn NA prop_forage to 0
  kits.att$prop_forage[is.na(kits.att$prop_forage)] <- 0
  
  
  ### Get sun elevation --------------------------------------------------------
  
  ## Get local times (useful for plotting later)
  kits.att %<>%
        # Add 1 second to times to circumvent POSIXct midnight bug
        mutate(datetime = datetime + 1,
              timezone = tz_lookup_coords(lat = col_lat[1], lon = col_lon[1], method = "fast", warn = FALSE),
              local_time = map2(datetime, timezone, ~format(.x, tz = .y))) %>% 
    unnest(local_time) %>%
    dplyr::select(-c(timezone))
  
  kits.att$local_time <- as.POSIXct(kits.att$local_time, format = "%Y-%m-%d %H:%M:%S") # this resets internal tz but not a problem
  kits.att$local_time <- kits.att$local_time - 1
  kits.att$datetime <- kits.att$datetime - 1

  # Calculate sun elevation for local time
  kits.att %<>%
  mutate(sun_altitude = getSunlightPosition(date = datetime, 
                                          col_lat[1], col_lon[1])[,4]) 
  
  
  # Convert sun elevation from radians to degrees
  kits.att$sun_altitude <- rad2deg(kits.att$sun_altitude)
 
  ### Add in temperature data --------------------------------------------------
  
  # Isolate columns and merge with kits.att
  col_temp <- subset(col_temp_data, colony == colonies[j])
  col_temp %<>% dplyr::select(-colony)
  kits.att <- merge(kits.att, col_temp, by = "datetime", all.x = TRUE)
  
  ## Add year and hour column (hour is LOCAL TIME)
  kits.att$year <- format(kits.att$datetime, "%Y")
  kits.att$local_hour <- strftime(kits.att$local_time, format = "%H")
  
  
  ### Clean up dataset ---------------------------------------------------------
  
  # Make sure only have data in May/June/July (breeding season)
  kits.att %<>% mutate(month = as.numeric(lubridate::month(datetime))) %>% 
    filter(month >= 5 & month <= 7) %>% dplyr::select(-month)
  
  # Change datetime to local time to avoid confusion
  kits.att$datetime <- kits.att$local_time
  kits.att$local_time <- NULL
  kits.att <- rename(kits.att, "hour" = "local_hour")
  
  # Add the dataset to the list
  col.list[[j]] <- kits.att
  
  }

  col_df <- do.call("rbind", col.list)
  kit.list[[i]] <- col_df
  
}
  
kits.att_all <- do.call("rbind", kit.list)
close(pb)

beepr::beep(8)

# Output the dataframe ---------------------------------------------------------
save(kits.att_all, file = "Data_inputs/BLKI_attendance_2013-2022.RData")
