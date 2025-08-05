
# GCD.HF ------------------------------------------------------------------

gcd.hf <- function(long1, lat1, long2, lat2) { 
  R <- 6371 # Earth mean radius [km]
  deg2rad <- function(deg) return(deg*pi/180)
  long1 <- deg2rad(long1)
  long2 <- deg2rad(long2)
  lat1 <- deg2rad(lat1)
  lat2 <- deg2rad(lat2)
  delta.long <- (long2 - long1)
  delta.lat <- (lat2 - lat1)
  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
  c <- 2 * asin(min(1,sqrt(a)))
  
  coo2 <- function(x){2 * asin(min(1,sqrt(x)))}
  c <- lapply(a, coo2)
  c <- do.call("c", c)
  d = R * c
  return(d) # Distance in km
}


# SAMPLE INDIVIDUAL DATA BY TIME PERIOD -----------------------------------

# Function to randomly sample a period
subsample_by_hours <- function(data, sample_time) {
  
  # Sort the data by individual and datetime
  data <- data %>% arrange(ring, datetime)
  
  # Initialize an empty data frame to store the sampled periods
  sampled_data <- data.frame()
  
  # Iterate over each individual
  unique_individuals <- unique(data$ring)
  
  for (individual in unique_individuals) {
    
    print(paste0("Processing individual ", individual, "of ", length(unique_individuals)))
    
    # Filter the data for the current individual
    individual_data <- data %>% filter(ring == individual)
    
    # Determine the number of observations 
    n_obs <- nrow(individual_data)
    
    if (n_obs >= sample_time) {
      
      # Randomly select a starting index 
      start_index <- sample(1:(n_obs - (sample_time - 1)), 1)
      
      # Extract the period
      sampled_period <- individual_data[start_index:(start_index + (sample_time - 1)), ]
      
      # Add the sampled period to the final data frame
      sampled_data <- rbind(sampled_data, sampled_period)
    }
  }
  
  return(sampled_data)
}



# FIND TRIP DURATION FROM PROCESSED IMMERSION DATA ------------------------

tripFinder <- function(datetime, immersion, rolling.immersion, thresh) {
  
  stts <- c()
  ends <- c()
  
  i <- 1
  r <- 1
  
  while (r < length(immersion)) {
    
    if (rolling.immersion[r] > thresh) {  # Point at which cross threshold
  
      pointR <- r   
     
      while (immersion[r] > 0 & r > 1) { # Move backwards to point where raw immersion = 0
        
        r <- r - 1 }  
      
      # Bug in R means POSIXct at midnight excludes time component
      posix.startDate <- format(as.POSIXct(datetime[r]), format="%Y-%m-%d %H:%M:%S", tz = "UTC")
      
      stts[i] <- as.character(posix.startDate) #  Label start of trip
      
      r <- pointR + 1     # Advance r by 1 
    
      
      while (rolling.immersion[r] > thresh & r < length(immersion)-3) {   # While immersion > threshold, move forwards to find end of trip
        
        r <- r + 1 }    
      
      r <- r + 1
      
      while (sum(immersion[r:(r+3)]) > 1 & r < length(immersion)-3) {  # Move forward until raw immersion = 0 for 30 mins
        
          r <- r + 1 }
      
      r <- r + 1  # Add 10 min commuting time penalty for return 
      
      # Bug in R means POSIXct at midnight excludes time component
      posix.endDate <- format(as.POSIXct(datetime[r]), format="%Y-%m-%d %H:%M:%S", tz = "UTC")
      
      ends[i] <- as.character(posix.endDate) # Set end time 
      
      r <- r + 1
      
      i <- i + 1 
      
    } else { r <- r + 1} } 
  
  trips <- data.frame(starttime = stts, endtime = ends)
  
  return(trips)
  
}

# FIND BEHAVIOUR FROM ACTIVITY LOGGER -------------------------------------

find_behaviour <- function(immersionValues) {
  
  behaviour <- rep("NA", length = length(immersionValues))
  
  next_immersion <- c(immersionValues[2:length(immersionValues)], NA)
  prev_immersion <- c(NA, immersionValues[1:length(immersionValues) - 1])
  
  behaviour[immersionValues <= max(immersionValues*0.10) &
              ( next_immersion == min(immersionValues) |
                  prev_immersion == min(immersionValues) )] <- "flight"
  behaviour[immersionValues >= max(immersionValues*0.90) &
             ( next_immersion == max(immersionValues) |
              prev_immersion == max(immersionValues) )] <- "rest"
  behaviour[behaviour == "NA"] <- "foraging"
  
  return(behaviour)
}



# LOAD RAW GLS DATA -------------------------------------------------------

load_raw_gls <- function(myfilename) {
  
  mydata <- read.csv(myfilename, header = F, skip = 1)
  
  if(ncol(mydata) > 4) { 
    mydata <- mydata[,c(6, 7, 3),] } else {
        mydata <- mydata[,c(2,4)]
        mydata$ring <- str_sub(myfilename,-11,-5) }
  
  colnames(mydata) <- c("datetime", "value", "ring")
  
  mydata
}


# CONVERT RADIANS AND BACK AGAIN ------------------------------------------

rad2deg <- function(rad) {(rad * 180) / (pi)}
deg2rad <- function(deg) {(deg * pi) / (180)}


# BLKI THEME ------------------------------------------------------------

require(ggplot2)

BLKI_theme <- theme_bw() + theme(axis.text.x = element_text(size = 12), 
                    axis.text.y = element_text(size = 12), 
                    axis.title.x = element_text(size = 14),
                    axis.title.y = element_text(size = 14),
                    plot.tag = element_text(size = 15),
                    plot.title = element_text(size = 16, face = "bold"))


# ## TASH'S HACKY DATE TIME FORMATTER -------------------------------------

time_sorter <- function(x){
  date1 <- sapply(strsplit(x, " "), "[[", 1)[1]
  time1 <- sapply(strsplit(x, " "), "[[", 2)[1]
  seconds <- length(strsplit(time1, ":")[[1]])
  
  if(grepl("-", date1) == TRUE) { 
    return(as.POSIXct(x, format = "%Y-%m-%d %H:%M:%S"))
  } else {
    
    if(nchar(date1) == 10 & seconds == 3 & nchar(sapply(strsplit(date1, "/"), "[[", 1)) == 4){
      return(as.POSIXct(x, format = "%Y/%m/%d %H:%M:%S")) 
    } else {
      
      if(nchar(date1) == 10 & seconds == 3 & nchar(sapply(strsplit(date1, "/"), "[[", 1)) == 2){
        return(as.POSIXct(x, format = "%d/%m/%Y %H:%M:%S"))
      } else {
        
        if(nchar(date1) == 8 & seconds == 3){
          return(as.POSIXct(x, format = "%d/%m/%y %H:%M:%S"))
        } else {
          
          if(nchar(date1) == 10 & seconds == 2){
            return(as.POSIXct(x, format = "%d/%m/%Y %H:%M"))
          } else {
            
            if(nchar(date1) == 8 & seconds == 2){
              return(as.POSIXct(x, format = "%d/%m/%y %H:%M"))
            }
          }
        }
      }
      
    }
  }
}



# REMOVE LONG ZERO RUNS ---------------------------------------------------

# These represent pre- or post-recording intervals

rm_non_recording <- function(df, threshold_hours = 48) {
  
  # Calculate number of rows corresponding to threshold_hours
  threshold_rows <- threshold_hours * 60 / 10
  
  df %<>% 
    arrange(datetime) %>%
    mutate(zero_flag = (immersion == 0),
           zero_group = cumsum(zero_flag != lag(zero_flag, default = FALSE)) ) %>%
    group_by(zero_group) %>%
    mutate(run_length = n()) 

  last_group <- df$zero_group[nrow(df)]
  
  # Filter logic: remove rows if they are in start or end zero run and exceed threshold
  df_filtered <- df %>%
    filter(!( (zero_group == 1 | zero_group == last_group)
      & run_length > threshold_rows ) ) %>%
    ungroup() %>%
    select(-zero_flag, -zero_group, -run_length)
  
  return(df_filtered)

  }




# GET PEAK STRENGTH -------------------------------------------------------


get_peak_strength <- function(df) {
  
  # restrict to lags > 1 (exclude lag 1 autocorrelation of 1.0)
  df_use <- df %>% filter(Lag > 1)
  
  # Find first peak after lag 1 (maximum ACF)
  peak_lag <- df_use$Lag[which.max(df_use$ACF)]
  peak_val <- max(df_use$ACF)
  
  # Get first trough after the peak
  df_after_peak <- df_use %>% filter(Lag > peak_lag)
  if (nrow(df_after_peak) == 0) return(NA)
  trough_val <- min(df_after_peak$ACF)
  
  # Peak strength = peak - trough
  return(peak_val - trough_val)
}
