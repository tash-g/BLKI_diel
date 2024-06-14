
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

# INTERPOLATE  ----------------------------------------------------

interpolate_gps <- function(lat, lon, datetime, int) {
  
  require(geosphere); require(pracma)
  
  df <- data.frame(cbind(as.character(datetime), lon, lat))
  colnames(df)[1] <- "datetime"
  df$datetime <- as.POSIXct(df$datetime, format = "%Y-%m-%d %H:%M:%S")
  df[,2:3] <- lapply(df[,2:3], as.numeric)
  
  ## Resample GPS to exact fixes
  df$x1 <-  c(abs(df$lon[1:(nrow(df)-1)]-df$lon[2:(nrow(df))]),1)
  df$x2 <-  c(abs(df$lat[1:(nrow(df)-1)]-df$lat[2:(nrow(df))]),1)
  df$x3 <- c(as.numeric(df$datetime)[1:(length(as.numeric(df$datetime))-1)]- as.numeric(df$datetime)[2:(length(as.numeric(df$datetime)))],-1)
  df <- df[which(df$x1 != 0 & df$x2 != 0 & df$x3 < 0),]
  
  df$distance_next_metres <- c(0,distHaversine(data.frame(df$lon[1:(nrow(df)-1)],
                                                          df$lat[1:(nrow(df)-1)]),
                                               data.frame(df$lon[2:(nrow(df))],
                                                          df$lat[2:(nrow(df))]))) 
  
  df$cumdist <- cumsum(df$distance_next_metres)
  df <- df[order(df$datetime),]
  
  nrt <- nrow(df)
  idt <- seq(from=df$datetime[1], to=df$datetime[nrt], by=int)
  ix = pchip(as.numeric(df$datetime), df$lon, as.numeric(idt))
  iy = pchip(as.numeric(df$datetime), df$lat, as.numeric(idt))
  icumdist = pchip(as.numeric(df$datetime), df$cumdist, as.numeric(idt))
  
  df <- data.frame(datetime = idt,
                   Longitude = ix,
                   Latitude = iy)
  
  return(df)
  
}


# CALC_TEST_STATISTICS ----------------------------------------------------

calc_test_statistics <- function(coefs, vcovs, b_pos1, b_pos2, conf.level = NULL) {
  
  if (is.null(conf.level)) {
    critical_value <- qnorm((1 + 0.95) / 2) } else {
      critical_value <- qnorm((1 + conf.level) / 2)
    }
  
  beta <- coefs[b_pos1,1] + coefs[b_pos2,1]
  se <- sqrt(coefs[b_pos1,2]^2 + cond_coefs[b_pos2,2]^2 + 2*cond_vcov[b_pos1,b_pos2]) 
  upper_conf <- beta + (critical_value * se)
  lower_conf <- beta - (critical_value * se)
  z_score <- beta / se
  p_value <- 2 * (1 - pnorm(abs(z_score)))
  
  print(data.frame(cbind(beta, se, lower_conf, upper_conf, z_score, p_value))[1,])
  
}

# LOGIT2PROB FUNCTION -----------------------------------------------------

logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
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


# GET DECIMAL FLOOR AND CEILING -------------------------------------------

# From : https://stackoverflow.com/questions/35807523/r-decimal-ceiling

floor_dec <- function(x, level=1) round(x - 5*10^(-level-1), level)
ceiling_dec <- function(x, level=1) round(x + 5*10^(-level-1), level)


# FIND ATTENDANCE BY HOUR -------------------------------------------------

## This is a clunky function but serves to streamline the randomisation analysis

att_by_hour <- function(dataset) {
  
  # Time spent at colony
  kits.att <- dataset %>% 
    group_by(ring, loc) %>%
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

  # Remove loc column - no longer needed
  kits.att$loc <- NULL

  # Convert date and hour column to posixct
  kits.att$datetime <- as.POSIXct(kits.att$datetime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

  ## Extend hours column to include missing hours ----------------------------
  kits.att <- kits.att %>%
    group_by(ring) %>%
    tidyr::complete(datetime = seq(min(datetime), max(datetime), by = 'hour'), fill = list(n = 0)) %>%
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
  
  # Output dataset
  kits.att

}


# FIND TRIP DURATION FROM PROCESSED IMMERSION DATA ------------------------

# For troubleshooting
# mydf <- mydf
# rolling.immersion <- mydf$rolling.percentage
# datetime <- mydf$datetime
# immersion <- mydf$immersion
# thresh <- thresh.immersion


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
      posix.startDate <- format(as.POSIXct(datetime[r], tz="UTC"), format="%Y-%m-%d %H:%M:%S")
      
      stts[i] <- as.character(posix.startDate) #  Label start of trip
      
      r <- pointR + 1     # Advance r by 1 
    
      
      while (rolling.immersion[r] > thresh & r < length(immersion)-3) {   # While immersion > threshold, move forwards to find end of trip
        
        r <- r + 1 }    
      
      r <- r + 1
      
      while (sum(immersion[r:(r+3)]) > 1 & r < length(immersion)-3) {  # Move forward until raw immersion = 0 for 30 mins
        
          r <- r + 1 }
      
      r <- r + 1  # Add 10 min commuting time penalty for return 
      
      # Bug in R means POSIXct at midnight excludes time component
      posix.endDate <- format(as.POSIXct(datetime[r], tz="UTC"), format="%Y-%m-%d %H:%M:%S")
      
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




# OVERDISPERSION FUNCTION (Ben Bolker) -----------------------------------------

overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}


