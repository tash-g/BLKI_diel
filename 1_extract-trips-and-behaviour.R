
## ---------------------------
##
## Script name: 1_extract-trips-and-behaviour.R
##
## Purpose of script: Extract trips from GLS data using immersion values
##
## Author: Dr. Natasha Gillies
##
## Created: 2022-08-19
##
## Email: gilliesne@gmail.com
##
## ---------------------------


# Load functions, packages, & data ---------------------------------------------

# Functions for GPS processing and plotting
source("BLKI_functions.R")

# Define the packages
packages <- c("ggplot2", "cowplot", "zoo", "dplyr", "lubridate", "suncalc",
              "lutz", "data.table", "magrittr")

# Install packages not yet installed - change lib to library path
#installed_packages <- packages %in% rownames(installed.packages())

#if (any(installed_packages == FALSE)) {
#  install.packages(packages[!installed_packages], lib = "C:/Users/libraryPath")
#}

# Load packages
invisible(lapply(packages, library, character.only = TRUE))

# Suppress dplyr summarise warning
options(dplyr.summarise.inform = FALSE)


# Extract trips from GLS -------------------------------------------------------

### Set parameters for trip identification -------------------------------------
rolling.width <- 20
thresh.immersion <- 15

### Load datasets --------------------------------------------------------------

# Huge amount of data - so loop through previously subsetted files individually

# Load metadata file
load("Data_inputs/BLKI_metadata.RData")
meta %<>% distinct(ring, .keep_all = TRUE) 

# Identify relevant files 
gls_files <- list.files("Data_workshop/GLS Processed")[grepl("GLS",list.files("Data_workshop/GLS Processed"))][-1]

# Loop through each file
gls_trips_list <- vector(mode = "list", length = length(gls_files))
gls_daily_behaviour_list <- vector(mode = "list", length = length(gls_files))
gls_hourly_behaviour_list <- vector(mod = "list", length = length(gls_files))


for (i in 1:length(gls_files)) {
  
  print(paste0("Processing file ", i, " of ", length(gls_files)))
  
  ### FILE LOOP - LOAD FILES ---------------------------------------------------
  
  load(paste0("Data_inputs/", gls_files[i]))

  # Make date a posixct variable
  gls$datetime <- as.POSIXct(gls$datetime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  
  # Merge metadata
  gls <- merge(gls, meta, by = "ring", all.x = TRUE)
  
  # Get ringyr for processing
  gls$ringYr <- paste(gls$ring, format(gls$datetime, "%Y"), sep = "_")

    ### Calculate rolling mean and add threshold immersion parameters ------------
  glsDat <- gls %>%
    group_by(ringYr) %>%
    # Get rolling immersion for each bird
    mutate(
      immersion.rolling = rollmean(immersion, rolling.width, na.pad = T),
      rolling.percentage = immersion.rolling / max(immersion.rolling, na.rm = T) * 100
    ) %>%
    data.frame()
  
  # Set NA rolling immersion to 0 so function will work
  glsDat$rolling.percentage[is.na(glsDat$rolling.percentage)] <- 0
  
  ### TRIP LOOP - IDENTIFY TRIPS FROM GLS --------------------------------------
  
  myKits <- unique(glsDat$ringYr)
  tripsList <- vector(mode = "list", length = length(myKits))
  dailyList <- vector(mode = "list", length = length(myKits))
  hourlyList <- vector(mode = "list", length = length(myKits))

  for (k in 1:length(myKits)) {
    
    mydf <- subset(glsDat, ringYr == myKits[k])
    mydf <- mydf[order(mydf$datetime),]
    
    # Add second to deal with datetime formatting errors
    mydf$datetime <- as.POSIXct(mydf$datetime) + 1
    
    # Catch GLS that didn't go out/have errors
    if (sum(mydf$immersion) < 5) next

    tripsID <- tripFinder(mydf$datetime,
                          mydf$immersion,
                          mydf$rolling.percentage,
                          thresh.immersion)

    mytrips <- data.frame(ring = mydf$ring[1],
                          type = "GLS",
                          start = tripsID[,1],
                          loc = "trip",
                          end = tripsID[,2])

    # Error if no end to trip
    mytrips <- subset(mytrips, !is.na(end))

    # Combine trips with short gaps between them
    mytrips$tripGap <- difftime(lead(mytrips$start), mytrips$end, units = "mins")
    mytrips$end[as.numeric(mytrips$tripGap) <= 20] <- NA
    mytrips$end <- na.locf(mytrips$end, fromLast = TRUE)
    mytrips <- subset(mytrips, !duplicated(end))
    mytrips$tripGap <- NULL

    # Add in colony rows
    myvisits <- data.frame(ring = mydf$ring[1],
                           type = "GLS",
                           start = c(mydf$datetime[1], mytrips$end),
                           loc = "colony",
                           end = c(mytrips$start, as.character(mydf$datetime[nrow(mydf)])))

    myvisits$start <- as.character(myvisits$start)

    # Combine datasets
    mytrips <- rbind(mytrips, myvisits)
    mytrips$start <- fasttime::fastPOSIXct(mytrips$start)
    mytrips$end <- fasttime::fastPOSIXct(mytrips$end)
    mytrips <- mytrips[order(mytrips$start),]
    mytrips$duration_mins <- difftime(mytrips$end, mytrips$start,
                                      units = "mins")
    mytrips <- subset(mytrips, as.numeric(duration_mins) > 0)

    # Add metadata
    mymeta <- mydf %>% dplyr::select(c(ring, colony, col_lat, col_lon)) %>% distinct(.keep_all = TRUE)
    mytrips <- merge(mytrips, mymeta, by = "ring", all.x = TRUE)

    # Remove additional second
    mytrips$start <- as.POSIXct(mytrips$start, tz = "UTC") - 1
    mytrips$end <- as.POSIXct(mytrips$end, tz = "UTC") - 1

    mytrips <- mytrips[order(mytrips$start),]
    
    tripsList[[k]] <- mytrips

    ### HOURLY DATAFRAME - IDENTIFY HOURLY BEHAVIOUR FROM GLS ------------------
    
    mydf$hour <- format(mydf$datetime, "%H")
    mydf$date <- as.Date(mydf$datetime)
    
    hourly_behaviour <- mydf %>% 
      group_by(colony, ring, date, hour) %>%
      summarise(total_imm = sum(immersion),
                n_flight = length(behaviour[behaviour == "flight"]),
                n_rest = length(behaviour[behaviour == "rest"]),
                n_forage = length(behaviour[behaviour == "foraging"]),
                for_poss = length(hour),
                prop_forage = n_forage/for_poss,
                rest_flight_ratio = n_rest/n_flight) %>%
      data.frame()
    
    mydf$hour <- NULL
    mydf$date <- NULL
    
    hourlyList[[k]] <- hourly_behaviour
    
    
    ### TRIP LOOP/DAILY DATAFRAME - IDENTIFY DAILY BEHAVIOUR FROM GLS ----------
    
    # Expand location column (colony vs trip) 
    mytrips[,c("start", "end")] <- lapply(mytrips[,c("start", "end")], as.POSIXct)
    
    attDF <- mytrips %>% 
      mutate(time = purrr::map2(start, end, seq, by = 'min')) %>%
      tidyr::unnest(time)  %>%
      dplyr::select(loc, time) %>%
      rename(datetime = time)
    
    # Merge into mydf
    setkey(setDT(attDF), datetime)
    setkey(setDT(mydf), datetime)
    
    mydf <- attDF[mydf, roll = "nearest"]
    
    # Set behaviour to colony if at colony
    mydf$behaviour[mydf$loc == "colony"] <- "colony"
    
    # Get time spent in each behaviour - 24 hour
    behavTime <- mydf %>%
      mutate(date = as.Date(datetime)) %>%
      group_by(date) %>%
      summarise(n_recs24 = n(),
                prop_forage24 = sum(behaviour == "foraging")/n_recs24, 
                prop_rest24 = sum(behaviour == "rest")/n_recs24,
                prop_flight24 = sum(behaviour == "flight")/n_recs24,
                col_att24 = sum(behaviour == "colony")/n_recs24,
                total_imm24 = sum(immersion)/n_recs24)
    
    # Get proportion of daylight spend in behaviour
    dayNightTime <- mydf %>% 
      mutate(sun_altitude = getSunlightPosition(date = datetime, 
                                                col_lat[1], col_lon[1])[,4]) %>%
      mutate(sun_altitude = rad2deg(sun_altitude), 
             date = as.Date(datetime),
             LoD = ifelse(sun_altitude > 0, "day", "night"),
             daylight.mins = as.numeric(difftime(getSunlightTimes(as.Date(datetime), lat = col_lat[1], lon = col_lon[1], keep = "sunset")[,4],
                                                 getSunlightTimes(as.Date(datetime), lat = col_lat[1], lon = col_lon[1], keep = "sunrise")[,4],
                                                 units = "mins"))) %>%
      group_by(date, LoD, daylight.mins) %>%
      summarise(n_recs = n(),
                n_recs.noCol = sum(behaviour != "colony"),
                prop_forage = sum(behaviour == "foraging")/n_recs, 
                prop_rest = sum(behaviour == "rest")/n_recs,
                prop_flight = sum(behaviour == "flight")/n_recs,
                prop_forage.noCol = sum(behaviour == "foraging")/n_recs.noCol, 
                prop_rest.noCol = sum(behaviour == "rest")/n_recs.noCol,
                prop_flight.noCol = sum(behaviour == "flight")/n_recs.noCol,
                col_att = sum(behaviour == "colony")/n_recs,
                total_imm = sum(immersion)/n_recs) %>%
      mutate(daylight.mins = ifelse(is.na(daylight.mins), 0, daylight.mins)) %>%
      data.frame()
    
    # If daylight mins = 0, means full daylight (no sunrise/sunset identified)
    dayNightTime$daylight.mins[dayNightTime$daylight.mins == 0] <- 1440
    
    # Change to wide format
    dayNightWide <- reshape(dayNightTime, idvar = c("date", "daylight.mins"), timevar = "LoD", direction = "wide")
    
    # Add dummy night variables where missing
    if (ncol(dayNightWide) == 12) {
      
      new_cols <- c("n_recs.night", "n_recs.noCol.night", "prop_forage.night", 
                    "prop_rest.night", "prop_flight.night", "prop_forage.noCol.night", 
                    "prop_rest.noCol.night", "prop_flight.noCol.night", "col_att.night", 
                    "total_imm.night")
      dayNightWide[,new_cols] <- 0
      
    }
    
    ## Combine all datasets into single daily dataframe
    dailyBehaviour <- merge(behavTime, dayNightWide, by = "date")
    
    dailyBehaviour %<>% 
      mutate(colony = mydf$colony[1],
             col_lat = mydf$col_lat[1],
             col_lon = mydf$col_lon[1],
             ring = mydf$ring[1]) %>% 
      relocate("ring")
    
    dailyList[[k]] <- dailyBehaviour
    
  }
  
  glsTrips <- do.call("rbind", tripsList)
  dayBehaviour <- do.call("rbind", dailyList)
  hourBehaviour <- do.call("rbind", hourlyList)
  
  gls_trips_list[[i]] <- glsTrips
  gls_daily_behaviour_list[[i]] <- dayBehaviour
  gls_hourly_behaviour_list[[i]] <- hourBehaviour
  
}

gls_trips <- data.table::rbindlist(gls_trips_list)
gls_daily <- data.table::rbindlist(gls_daily_behaviour_list)
gls_hourly <- data.table::rbindlist(gls_hourly_behaviour_list)

### Output the data ----------------------------------------------------------
save(gls_trips, 
     file = "Data_inputs/BLKI_gls-trips.RData", 
     row.names = F)
save(gls_daily, 
     file = "Data_inputs/BLKI_gls-daily-behaviour.RData", 
     row.names = F)
save(gls_hourly, 
     file = "Data_inputs/BLKI_gls-hourly-behaviour.RData", 
     row.names = F)



