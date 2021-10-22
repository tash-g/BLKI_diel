
# INTERPOLATE GPS ---------------------------------------------------------

interpolate_gps <- function(lat, lon, datetime, int) {
  
  df <- data.frame(cbind(as.character(datetime), lon, lat))
  colnames(df)[1] <- "datetime"
  df$datetime <- as.POSIXct(df$datetime, format = "%Y-%m-%d %H:%M:%S")
  df[,2:3] <- lapply(df[,2:3], as.numeric)
  
  ## Resample GPS to exact 15-minute fixes
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
  
  ## Calculate speed
  n <- nrow(df)
  df$dist_next <- c(gcd.hf(df$Longitude[2:n], df$Latitude[2:n], df$Longitude[1:(n-1)], df$Latitude[1:(n-1)]),NA)
  
  dt = c(as.numeric(difftime(df$datetime[2:n], df$datetime[1:(n-1)], units = "secs")), NA)
  df$calc_speed <- df$dist_next*1000/dt
  
  return(df)
  
}



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



# ## GREAT CIRCLE DISTANCE CALC -------------------------------------------

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




# PLOTTING THEME ----------------------------------------------------------

theme_tash <- function() {
  
  theme_bw(base_size =15) +
    theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           axis.text.x = element_text(colour = "black", angle = 0, hjust = 0.5, size=17),
           axis.title.y = element_text(colour = "black", size = 17, margin = margin(t = 0, r = 12, b = 0, l = 0)),
           axis.text.y = element_text(colour = "black", angle = 0, hjust = 0.5, size=17),
           axis.title.x = element_text(colour = "black", size = 17,margin = margin(t = 12, r = 0, b =0, l = 0)),
           legend.background = element_blank())
}

