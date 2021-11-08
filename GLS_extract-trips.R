### PURPOSE: This script extracts foraging trips from the GLS data (within folder: GLS recovered 2018)

### GLS files are in 10 minute bins with following columns:
### X = row number; fileID = device ID, year, ID; ring; season; datetime; wetdry; date; year

### Notes: for light sub-sampling stage, would be better to average every 2 rows

library(dplyr); library(cowplot)

# LOAD DATA ---------------------------------------------------------------

### Load immersion files

# Assign folder
gls_folder.imm <- "./GLS-GPS Comparison/GLS recovered 2018/Activity files/"

filenames.imm <- paste0(gls_folder.imm, list.files(gls_folder.imm, pattern="*.csv"))
glsFiles.imm <- lapply(filenames.imm, read.csv)
glsData.imm <- do.call("rbind", glsFiles.imm)


### Load light files

# Assign folder
gls_folder.lig <- "./GLS-GPS Comparison/GLS recovered 2018/Light files/"

filenames.lig <- paste0(gls_folder.lig, list.files(gls_folder.lig, pattern="*.csv"))
glsFiles.lig <- lapply(filenames.lig, read.csv)
glsData.lig <- do.call("rbind", glsFiles.lig)

## OUTSTANDING: Average every 2 rows to bring into same sampling frequency as immersion (10 min vs 5 min bins)


### Merge datasets

## Isolate variables
glsData.lig <- glsData.lig[,c(3,6,7),] # ring, datetime, light
glsData.imm <- glsData.imm[,c(3,6,7),] # ring, datetime, wetdry

glsData <- merge(glsData.imm, glsData.lig, by = c("ring", "datetime"), all.x = TRUE)

## Process dataset
glsData$ring <- as.character(glsData$ring)
glsData$datetime <- as.POSIXct(glsData$datetime, format = "%Y-%m-%d %H:%M:%S")


### Write into csv
write.csv(glsData, "gls_summarised.csv", row.names = F)


# ALIGNED PLOTS FOR EACH INDIVIDUAL ---------------------------------------

### Plot immersion, light, and GPS attendance for each individual

myKits <- unique(glsData$ring)

for (i in 1:length(myKits)) {
  
  gpsDat <- subset(tripSummary, ring == myKits[i])
  
  # Get min and max times
  mn <- min(gpsDat$start)
  mx <- max(na.omit(gpsDat$end))
  
  # Subset GLS data
  glsDat <- subset(glsData, ring == myKits[i] & datetime > mn & datetime < mx)
  
  
  # Build the plots
  wetdryPlot <- ggplot(glsDat, aes(x = datetime, y = wetdry)) + geom_line() + xlab("") + ylab("Immersion")
  lightPlot <- ggplot(glsDat, aes(x = datetime, y = light)) + geom_line() + xlab("") + ylab("Light")
  gpsPlot <- ggplot(gpsDat, aes(x = start, y = ring, colour = loc)) + 
    geom_segment(aes(x = start, xend = end, yend = ring), size = 6) +
    theme_bw() + theme(legend.position = c(0.95,0.8), legend.key = element_rect(fill = "transparent"),
                       legend.background = element_rect(fill="transparent"),
          axis.text.y=element_blank(), axis.title.y=element_blank()) +
    xlab("Datetime") + 
    scale_colour_manual(values = c("#FFC107", "#1E88E5", "white"))
  gpsPlot <- addSmallLegend(gpsPlot)
  
  alignedPlot <- plot_grid(lightPlot + scale_x_datetime(limits=c(mn, mx)),
            wetdryPlot + scale_x_datetime(limits=c(mn, mx)),
            gpsPlot + scale_x_datetime(limits=c(mn, mx)), ncol = 1, align = "v")
  
  
  png(filename = paste0("Figures/glsPlot_", glsDat$ring[i], ".png"),
      width = 21, height = 18, units = "cm", res = 300)
  print(alignedPlot)
  dev.off()

  
}




# ROLLING MEAN PLOTS ------------------------------------------------------

### AIM: To get visitation to nearest 30 mins

### Rolling mean

library(zoo)

glsData <- glsData %>% group_by(ring) %>%
                mutate(wetdry.rolling = rollmean(wetdry, 3, fill = NA),
                       light.rolling = rollmean(light, 3, fill = NA)) %>%
                data.frame()


myKits <- unique(glsData$ring)

for (i in 1:length(myKits)) {
  
  gpsDat <- subset(tripSummary, ring == myKits[i])
  
  # Get min and max times
  mn <- min(gpsDat$start)
  mx <- max(na.omit(gpsDat$end))
  
  # Subset GLS data
  glsDat <- subset(glsData, ring == myKits[i] & datetime > mn & datetime < mx)
  
  
  # Build the plots
  wetdryPlot <- ggplot(glsDat, aes(x = datetime, y = wetdry.rolling)) + geom_line() + xlab("") + ylab("Immersion")
  lightPlot <- ggplot(glsDat, aes(x = datetime, y = light.rolling)) + geom_line() + xlab("") + ylab("Light")
  gpsPlot <- ggplot(gpsDat, aes(x = start, y = ring, colour = loc)) + 
    geom_segment(aes(x = start, xend = end, yend = ring), size = 6) +
    theme_bw() + theme(legend.position = c(0.95,0.8), legend.key = element_rect(fill = "transparent"),
                       legend.background = element_rect(fill="transparent"),
                       axis.text.y=element_blank(), axis.title.y=element_blank()) +
    xlab("Datetime") + 
    scale_colour_manual(values = c("#FFC107", "#1E88E5", "white"))
  gpsPlot <- addSmallLegend(gpsPlot)
  
  alignedPlot <- plot_grid(lightPlot + scale_x_datetime(limits=c(mn, mx)),
                           wetdryPlot + scale_x_datetime(limits=c(mn, mx)),
                           gpsPlot + scale_x_datetime(limits=c(mn, mx)), ncol = 1, align = "v")
  
  
  png(filename = paste0("Figures/glsPlot_", glsDat$ring[i], "_rolling-30mins.png"),
      width = 21, height = 18, units = "cm", res = 300)
  print(alignedPlot)
  dev.off()
  
  
}


### Summed every 30 minutes (i.e 30 minute bins)

n <- 3

glsData.dates <- glsData %>% group_by(ring) %>%
                         summarise(datetime2 = datetime[seq(1, length(datetime), 3)]) %>%
                         data.frame()

glsData.summed <- glsData %>%
                    group_by(ring, indx = gl(ceiling(length(ring)/n), n, length(ring))) %>%
                    summarise(light.sum = sum(light),
                              wetdry.sum = sum(wetdry)) %>%
                   data.frame()

glsData.summed$datetime <- glsData.dates$datetime2
glsData.summed$indx <- NULL


myKits <- unique(glsData.summed$ring)

for (i in 1:length(myKits)) {
  
  gpsDat <- subset(tripSummary, ring == myKits[i])
  
  # Get min and max times
  mn <- min(gpsDat$start)
  mx <- max(na.omit(gpsDat$end))
  
  # Subset GLS data
  glsDat <- subset(glsData.summed, ring == myKits[i] & datetime > mn & datetime < mx)
  
  
  # Build the plots
  wetdryPlot <- ggplot(glsDat, aes(x = datetime, y = wetdry.sum)) + geom_line() + xlab("") + ylab("Immersion")
  lightPlot <- ggplot(glsDat, aes(x = datetime, y = light.sum)) + geom_line() + xlab("") + ylab("Light")
  gpsPlot <- ggplot(gpsDat, aes(x = start, y = ring, colour = loc)) + 
    geom_segment(aes(x = start, xend = end, yend = ring), size = 6) +
    theme_bw() + theme(legend.position = c(0.95,0.8), legend.key = element_rect(fill = "transparent"),
                       legend.background = element_rect(fill="transparent"),
                       axis.text.y=element_blank(), axis.title.y=element_blank()) +
    xlab("Datetime") + 
    scale_colour_manual(values = c("#FFC107", "#1E88E5", "white"))
  gpsPlot <- addSmallLegend(gpsPlot)
  
  alignedPlot <- plot_grid(lightPlot + scale_x_datetime(limits=c(mn, mx)),
                           wetdryPlot + scale_x_datetime(limits=c(mn, mx)),
                           gpsPlot + scale_x_datetime(limits=c(mn, mx)), ncol = 1, align = "v")
  
  
  png(filename = paste0("Figures/glsPlot_", glsDat$ring[i], "_summed-30mins.png"),
      width = 21, height = 18, units = "cm", res = 300)
  print(alignedPlot)
  dev.off()
  
  
}




# LABEL START/END OF TRIPS USING GLS --------------------------------------

## Extract indices satisfying some condition 
## e.g. light == 60 & wetdry > 0
## label plots with these

## Set parameters outside the function for easy adjustment
thresh.wetdry <- 0
thresh.light <- 64

## Plotting approach to examine how worked

myKits <- unique(glsData$ring)

for (i in 1:length(myKits)) {
 
  gpsDat <- subset(tripSummary, ring == myKits[i])
  gpsDat <- gpsDat[,c(1,3,7,8,9),] # ring, start, loc, end, dur.mins
  
  
  # Get min and max times
  mn <- min(gpsDat$start)
  mx <- max(na.omit(gpsDat$end))
  
  # Subset GLS data
  glsDat <- subset(glsData, ring == myKits[i] & datetime > mn & datetime < mx)
  
  
  ## Characterise location based on light
  glsDat$loc <- ifelse(glsDat$light.rolling == thresh.light, "trip", "colony")
  
  # REFINE: Any immersion means not at nest
  glsDat$loc[glsDat$wetdry.rolling > 0] <- "trip"
  
  # Get start and end times
  runs <- rle(glsDat$loc == "trip")
  tripstarts <- which(runs$values == TRUE)
  
  tripfinder <- function(x) {sum(runs$lengths[1:(x-1)])}  ## Ask Ollie to explain how this works 
  endfinder <- function(x) {sum(runs$lengths[1:(x)])}
  
  tripst <- do.call("c", lapply(tripstarts, tripfinder))
  tripend <- do.call("c", lapply(tripstarts, endfinder))
  
  starttime <- glsDat$datetime[tripst]
  endtime <- glsDat$datetime[tripend]
  
  tripdf <- data.frame(ring = "GLS", start = starttime, loc = "trip", end = endtime)
  tripdf$dur.mins <- as.numeric(difftime(tripdf$end, tripdf$start, units = "mins"))
  
  # Trips must be at least 1 hour
  tripdf$loc[tripdf$dur.mins < 60] <- "unknown"

  gpsDat$ring <- "GPS"
  gpsDat <- rbind(gpsDat, tripdf)
  
  # Build the plots
  wetdryPlot <- ggplot(glsDat, aes(x = datetime, y = wetdry.rolling)) + geom_line() +
    xlab("") + ylab("Immersion") + theme_bw()
  
  lightPlot <- ggplot(glsDat, aes(x = datetime, y = light.rolling)) + geom_line() +
    xlab("") + ylab("Light") + theme_bw()
 
   gpsPlot <- ggplot(gpsDat, aes(x = start, y = ring, colour = loc)) + 
    geom_segment(aes(x = start, xend = end, yend = ring), size = 6) +
    theme_bw() + theme(legend.position = "none") +
    xlab("Datetime") + ylab("Data Source") + 
    scale_colour_manual(values = c("#FFC107", "#1E88E5", "white"))
  
   
  alignedPlot <- plot_grid(lightPlot + scale_x_datetime(limits=c(mn, mx)),
                           wetdryPlot + scale_x_datetime(limits=c(mn, mx)),
                           gpsPlot + scale_x_datetime(limits=c(mn, mx)), ncol = 1, align = "v")
  
  
  png(filename = paste0("Figures/glsPlot_", glsDat$ring[i], "_validation.png"),
      width = 21, height = 18, units = "cm", res = 300)
  print(alignedPlot)
  dev.off()
  
  
}



#### Diving code to adapt

runs <- rle(divedat$depth < -0.5)
divestarts <- which(runs$values == TRUE) ## RLE has split data into chunks where TRUE = depth is more than -0.5 
## and FALSE = depth is less than -0.5. We only want the true values. 
divefinder <- function(x) {sum(runs$lengths[1:(x-1)])}  ## Ask Ollie to explain how this works 
endfinder <- function(x) {sum(runs$lengths[1:(x)])}

divest <- do.call("c", lapply(divestarts, divefinder))
divend <- do.call("c", lapply(divestarts, endfinder))

starttime <- divedat$dates[divest]
endtime <- divedat$dates[divend]

outdf <- data.frame(starttime =  divedat$dates[divest], endtime = divedat$dates[divend])
maxdepth <- c()

for(dive in 1:length(divest)) {
  
  maxdepth[dive] <- min(divedat$depth[divest[dive]:divend[dive]])
  
}

outdf$maxdepth <- maxdepth




# FUNCTIONS ---------------------------------------------------------------

### Small legend function

addSmallLegend <- function(myPlot, pointSize = 4, textSize = 5, spaceLegend = 1) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.text  = element_text(size = textSize),
          legend.title = element_blank(),
          legend.key.size = unit(spaceLegend, "lines"))
}

addSmallLegend(gpsPlot)
