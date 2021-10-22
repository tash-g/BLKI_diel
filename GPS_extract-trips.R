### PURPOSE: This script extracts foraging trips from the GPS data (summary file: Grumant2017_GPS.csv)
### Summary file already contains distance from colony and calculated trip versus colony

### Outstanding questions:
### 1) Is colony-or-trip column accurate? E.g. can see 20 minute trips?

### Notes:
### - Not all GPS are on the same fix schedule


# Preamble ----------------------------------------------------------------

### Load functions/packages

library(dplyr); library(ggplot2)

### Load and clean the data

gpsDat <- read.csv("GLS-GPS Comparison/GPS from 2017/Grumant2017_GPS.csv")
gpsDat$datetime <- as.POSIXct(gpsDat$datetime, format = "%Y-%m-%d %H:%M:%S")



# Get trip metrics for each individual -----------------------------------------

### 1) Label individual colony versus trip bouts

gpsDat <- gpsDat %>% group_by(ring) %>%
                      mutate(bout = rep(1:length(rle(ColonyorTrip)$lengths), rle(ColonyorTrip)$lengths),
                             boutID = paste(ColonyorTrip, ceiling(bout/2), sep = "_"),
                             fixtime.mins = as.numeric(difftime(lead(datetime),datetime, units = "mins"))) %>%
                      data.frame()

gpsDat$loc <- ifelse(gpsDat$fixtime.mins > 1440, "unknown", as.character(gpsDat$ColonyorTrip))
gpsDat$loc[is.na(gpsDat$loc)] <- gpsDat$ColonyorTrip[is.na(gpsDat$loc)]

### 2) Get start and end time of each bout 

tripSummary <- gpsDat %>% group_by(ring, boutID) %>%
                          summarise(start = datetime[1], 
                                    colony = colony[1],
                                    stage = stage[1],
                                    year = year[1],
                                    loc = loc[1]) %>%
                          data.frame()

tripSummary <- tripSummary[order(tripSummary$ring, tripSummary$start),]

tripSummary <- tripSummary %>% group_by(ring) %>%
                               mutate(end = lead(start),
                                     dur.mins = as.numeric(difftime(end, start))) %>%
                               data.frame()


# Plot visitation timeline for each individual ---------------------------------------

tripSummary$ring <- as.factor(as.character(tripSummary$ring))

## Still some unrealistically long trips
tripSummary$loc[tripSummary$dur.mins > 5000] <- "unknown"

ggplot(tripSummary, aes(x = start, y = ring, colour = loc)) + 
  geom_segment(aes(x = start, xend = end, yend = ring), size = 6) +
  theme_bw() +
  ylab("Ring") + xlab("Time") + 
  scale_colour_manual(values = c("#FFC107", "#1E88E5", "white"))

