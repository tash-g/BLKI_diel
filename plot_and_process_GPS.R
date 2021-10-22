### PURPOSE: This script processes the raw GPS data and plots them

source("gps_functions.R")

library(dplyr); library(ggplot2); library(raster); library(rasterVis); library(scales); library(rgeos)


# Load in data ------------------------------------------------------------

allfiles <- list.files(path = "./GPS Pyramiden 2018", pattern = '.csv')

### Load and process each file 

gps_list <- vector(mode = "list", length = length(allfiles))

for(i in 1:length(allfiles)){
  
  bird <- data.table::fread(file = paste0("./GPS Pyramiden 2018/",allfiles[i], sep = ""))
  
  ## Clean the data
  bird$datetime <- paste(bird$Date, bird$Time)
  bird$datetime <- fasttime::fastPOSIXct(bird$datetime)

  
  ## Calculate speed and distance from home
  
  # Get home coords
  bird$home_lat <- 78.65562
  bird$home_lon <- 16.37982
  
  # Calc distance to home
  bird$dist_home <- gcd.hf(bird$Longitude, bird$Latitude, bird$home_lon, bird$home_lat)
  
  # Calc speed
  n <- nrow(bird)
  
  bird$dist_next <- c(gcd.hf(bird$Longitude[2:n], bird$Latitude[2:n], bird$Longitude[1:(n-1)], bird$Latitude[1:(n-1)]),NA)
  dt = c(as.numeric(difftime(bird$datetime[2:n], bird$datetime[1:(n-1)], units = "secs")), NA)
  bird$calc_speed <- bird$dist_next*1000/dt
  
  
  ## Add filename metadata to the lon/lat/date info
  gps_list[[i]] <- bird
  gps_list[[i]]$Ring <- rep(strsplit(allfiles[i], "\\W")[[1]][2], nrow(bird))
  
  print(i)
}

# bind em up!
alltracks <- data.table::rbindlist(gps_list)


### Add meta data

meta <- read.csv("GPSmetafile_Pyramiden2018.csv")

# Isolate relevant columns
meta <- meta[,c(2,3,4,5,7,8,10),]
meta$deployment <- paste(meta$Deployment.date, meta$Deployment.time)
meta$deployment <- fasttime::fastPOSIXct(meta$deployment)

meta$retrieval <- paste(meta$Retrieval.date, meta$Retrieval.time)
meta$retrieval <- fasttime::fastPOSIXct(meta$retrieval)

meta[,c(1,2,5,6)] <- NULL
meta$Ring <- as.character(meta$Ring)

## Merge with GPS
alltracks <- merge(alltracks, meta, by = "Ring", all.x = T)

# Export
write.csv(alltracks, "gps_processed.csv", row.names = F)


# Split into trips --------------------------------------------------------

## Simple plot to determine threshold for colony departure
# Need more refined condition e.g. speed? 

gpsDat.df$threshold <- ifelse(gpsDat.df$dist_home > 3, "away", "home")

ggplot(data = sval3) + 
  geom_sf(fill = "cadetblue", colour = "grey") + 
  coord_sf(crs = proj.utm, xlim = c(-15377.9, 13524.07), ylim = c(-33192.5, 17386.9),
           label_axes = list(top = "E", left = "N", bottom = "E", right = "N")) +
  geom_path(aes(x = Longitude, y = Latitude, col = threshold), dat = gpsDat.df) +
  theme(legend.position = "none")



## Mixture model to see if trimodal distribution in distance from home

mixmd <- mixtools::normalmixEM(alltracks$dist_home, k = 3)
plot(mixmd, which = 2)

# Plot the data -----------------------------------------------------------

library(raster); library(rgdal)

## SET PROJECTIONS 

# Load shape file and make sure CRS matches coordinate data and projection
proj.dec <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
proj.utm <- "+proj=laea +lat_0=78.65562 +lon_0=16.37982 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

## Set projection of GPS data
gpsDat.proj <- alltracks
coordinates(gpsDat.proj) <- ~ Longitude+Latitude
proj4string(gpsDat.proj) <- proj.dec

gpsDat.sp <- spTransform(gpsDat.proj, CRS(proj.utm))
gpsDat.df <- as.data.frame(gpsDat.sp)

## Colony projection (for plotting)
colony <- data.frame(Lon = 16.37982, Lat = 78.65562)

coordinates(colony) <- ~Lon+Lat
proj4string(colony) <- proj.dec
colony.sp <- spTransform(colony, CRS(proj.utm))
colony.df <- as.data.frame(colony.sp)

## Project back to original coordinates
gpsDat.sp2 <- spTransform(gpsDat.sp, CRS(proj.dec))
gpsDat.df2 <- as.data.frame(gpsDat.sp2)


## BUILD PLOT 

## Base global plot
library(rnaturalearth); library(rnaturalearthdata)

world <- ne_countries(scale = "medium", returnclass = "sf")
world2 <- sf::st_transform(world, crs = proj.utm)

## ZOomed out 
png("Figures/full_map.png", width = 7, height = 7, units = "in", res = 700)
ggplot(data = world2) + 
  geom_sf(fill = "cadetblue", colour = "grey") +
  ggspatial::annotation_scale(location = "bl", width_hint = 0.25, style = "bar") +
  coord_sf(crs = proj.utm, xlim = c(-305377.9, 203524.07), ylim = c(-373192.5, 257386.9),
           label_axes = list(top = "E", left = "N", bottom = "E", right = "N")) +
  geom_path(aes(x = Longitude, y = Latitude, col = Ring), dat = gpsDat.df) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white"), 
        panel.grid.major = element_line(colour = "grey80"),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x.bottom = element_blank(), axis.title.x.bottom = element_blank(), 
        axis.text.y.right = element_blank(), axis.title.y.left = element_blank()) +
  annotate("point", shape = 17, colony.df$Lon, colony.df$Lat) +
  annotate("text", label = "Pyramiden", colony.df$Lon, (colony.df$Lat+25000))
dev.off()

## Zoomed in - need to add glacier
# Download high res Svalbard map
sval <- getData('GADM', country='SJM', level=1)
sval1 <- subset(sval, NAME_1 == "Svalbard")

# Set projection
sval2 <- sf::st_as_sf(sval1)
sval3 <- sf::st_transform(sval2, crs = proj.utm)

# Plot
png("Figures/zoomed_map.png", width = 7, height = 7, units = "in", res = 700)
ggplot(data = sval3) + 
  geom_sf(fill = "cadetblue", colour = "grey") + 
  coord_sf(crs = proj.utm, xlim = c(-15377.9, 13524.07), ylim = c(-33192.5, 17386.9),
           label_axes = list(top = "E", left = "N", bottom = "E", right = "N")) +
  geom_path(aes(x = Longitude, y = Latitude, col = Ring), dat = gpsDat.df) +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white"), 
        panel.grid.major = element_line(colour = "grey80"),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x.bottom = element_blank(), axis.title.x.bottom = element_blank(), 
        axis.text.y.right = element_blank(), axis.title.y.left = element_blank()) +
  annotate("point", shape = 17, size = 3, colony.df$Lon, colony.df$Lat) +
  annotate("text", label = "Pyramiden", colony.df$Lon, (colony.df$Lat+3000))
dev.off()




# DETAILED ZOOMED MAP -------------------------------------------------------

## Building detailed topography map to view glaciar use



## Get elevation
elevation <- elevatr::get_elev_raster(sval1, z = 9)
elevation_poly <- rasterToPolygons(elevation)

# Define colour ramp palette
colr <- colorRampPalette(ColorBrewer::brewer.pal(11, 'RdYlBu'))

## ggplot

# Transform data




## Plot elevation

elevation.pts <- rasterToPoints(elevation, spatial = T)
elevation.df <- data.frame(elevation.pts)

ggplot() + geom_raster(data = elevation.df, aes(x = x, y = y, fill =))