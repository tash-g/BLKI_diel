## ---------------------------
##
## Script name: 4_plot-actograms.R
##
## Purpose of script: Analyse nest visitation using attendance datasets
##
## Author: Dr. N. P. Huffeldt, adapted from S. Sjöberg (Lund University), adapted
## by N Gillies
##
## Created: 2024-04-04
## 
## ---------------------------


# Preamble =====================================================================

# Define the packages
packages <- c("tidyverse", "lubridate", "cowplot", "viridis", "stringr")

# Install packages not yet installed - change lib to library path
# installed_packages <- packages %in% rownames(installed.packages())
# 
#  if (any(installed_packages == FALSE)) {
#   install.packages(packages[!installed_packages])
#  }

# Load packages and functions
invisible(lapply(packages, library, character.only = TRUE))



# Load data --------------------------------------------------------------------

load("Data_inputs/BLKI_actogram_data.RData")

# Create acc_List
ACC_LIST <- split(gls.hourly, f = gls.hourly$ring)

# names each list with the ID
names(ACC_LIST) <- sapply(ACC_LIST, function(x) x$ring[1])
T_names <- names(ACC_LIST)


### Create dataframes ----------------------------------------------------------

Acto_list <- list()

for(i in 1:length(ACC_LIST)){
  startdt <- as.POSIXct(ACC_LIST[[i]]$date[1], tz="UTC") #define the first date to use data from
  df.a <- ACC_LIST[[i]] %>% 
    filter(date > startdt)
  df.a$Date2 <- as.Date(df.a$date)
  df.a$hour <- as.numeric(df.a$hour)

  df.b <- ACC_LIST[[i]] %>% 
    filter(date > startdt)
  df.b$Date2 <- as.Date(df.b$date) -1
  df.b$hour <- as.numeric(df.b$hour) + 24
  
  df.c <- rbind(df.a, df.b)
  df.c$mean.act <- df.c$total_imm/12
  df.c$Rday <- as.numeric(as.Date(df.c$Date2) - as.Date(startdt)-1 + (yday(startdt)+1))
  

  x <- list(df.c)
  DF_name <- names(ACC_LIST[i])
  names(x) <- DF_name
  Acto_list <- append(Acto_list, x)

}  



### Plot vertical Actograms-----------------------------------------------------

## Select ring numbers
selected_list <- readxl::read_xlsx("example_actograms.xlsx")
selected_birds.a <- c(selected_list$primary, selected_list$secondary) %>% 
  strsplit("_") %>% sapply("[", c(4)) 
selected_birds.b <- c(selected_list$primary, selected_list$secondary) %>% 
  strsplit("_") %>% sapply("[", c(5)) 
selected_birds <- paste(selected_birds.a, selected_birds.b, sep = "_")
selected_birds <- gsub("_NA", "", selected_birds)


#### Plot selection of birds ---------------------------------------------------

acto_name.list <- vector(mode = "list", length = length(selected_birds))

for(i in 1:length(selected_birds)){
  
  df.c <- Acto_list[grepl(selected_birds[i], names(Acto_list))]
  df.c <- do.call("rbind", df.c)
  rownames(df.c) <- NULL
  
  df.c$Rdayc <- as.factor(df.c$Rday)
  df.c$month <- format(df.c$date,"%m")
  
  label <- c(ifelse(unique(ceiling(df.c$Rday/30.5))<13, unique(ceiling(df.c$Rday/30.5)), unique(ceiling(df.c$Rday/30.5))-12))
  label <- c(month.abb[(ifelse(label<13, label, label-12))])
  
  breaks <- c(as.character(unique(ceiling(df.c$Rday/30.5)*30.5-15)))
  
  actogram <- ggplot()+
    geom_tile(aes(x = hour, y = Rdayc, fill = mean.act), data=df.c) + # colors tiles based on (mean) hourly activity. Must make column for hourly activity before running script
    scale_fill_gradientn(colours = c("white", "#009688", "#00796B", "#004d40", "black")) +
    geom_tile(data = df.c[is.na(df.c$Activity_NAs),], aes(x = hour, y = Rdayc), fill = "firebrick3", group = 1)  + #  colors NAs red
    coord_fixed(0.15) + # this adjusts the width of saved figure (smaller = wider)
    scale_x_continuous(name = "Time of day (Local)", expand = c(0, 0), limits=c(-0.001,47.001),
                       breaks = c( 6, 18, 30, 42), labels = c("06", "18", "06", "18")) +
    scale_y_discrete(name = NULL, expand=c(0,0), limits=rev, breaks = breaks, labels=label)+
    theme_bw(base_size = 12)+
    theme(panel.border = element_rect(fill = NA), legend.key.size = unit(4, "mm"), 
          plot.title = element_text(hjust = 0.5,), legend.position = "none",
          plot.subtitle = element_text(face = "italic", size = 8, hjust = 0.5)) +
    labs(title = paste0(str_to_title(df.c$colony[1]), "\n", df.c$lat, "°N"), 
         fill = "Mean\nimmersion", size = 10, family = "Arial",
         subtitle = df.c$ring[1])
  
  acto_name <- paste0("act_", df.c$ring[1], "_", df.c$lat[1], "_", df.c$colony[1])
  assign(acto_name, actogram)
  
  acto_name.list[[i]] <- acto_name
  
}

acto_name.df <- do.call("rbind", acto_name.list)
acto_name.df <- data.frame(acto_name.df)
acto_name.df$colony <- acto_name.df$acto_name.df %>% strsplit("_") %>% sapply("[", 5)
acto_name.df$lat <- acto_name.df$acto_name.df %>% strsplit("_") %>% sapply("[", 4)
acto_name.df$lat[acto_name.df$acto_name.df == "act_104401809_47.25914_witless bay"] <- 47.24
acto_name.df$colony[acto_name.df$acto_name.df == "act_104401809_47.25914_witless bay"] <- "Witless Bay"
acto_name.df <- acto_name.df[order(acto_name.df$lat),]



#### Big plots with them all 
multiPlot1 <- plot_grid(get(acto_name.df$acto_name.df[1])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,-1,0,0), "cm")),
          get(acto_name.df$acto_name.df[2])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,-1,0,0), "cm")),
          get(acto_name.df$acto_name.df[3])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,0,0,0), "cm")),
          get(acto_name.df$acto_name.df[4])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,0,0,0), "cm")),
          
          get(acto_name.df$acto_name.df[5])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,-1,0,0), "cm")),
          get(acto_name.df$acto_name.df[6])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,-1,0,0), "cm")),
          get(acto_name.df$acto_name.df[7])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,0,0,0), "cm")),
          get(acto_name.df$acto_name.df[8])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,0,0,0), "cm")),
                    
          get(acto_name.df$acto_name.df[9])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,-1,0,0), "cm")),
          get(acto_name.df$acto_name.df[10])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,-1,0,0), "cm")),
          get(acto_name.df$acto_name.df[11])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,0,0,0), "cm")),
          get(acto_name.df$acto_name.df[12])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,0,0,0), "cm")),
          
          get(acto_name.df$acto_name.df[13])+ theme(plot.margin = unit(c(0,-1,0,0), "cm")),
          get(acto_name.df$acto_name.df[14])+ theme(plot.margin = unit(c(0,-1,0,0), "cm")),
          get(acto_name.df$acto_name.df[15])+ theme(plot.margin = unit(c(0,0,0,0), "cm")),
          get(acto_name.df$acto_name.df[16])+ theme(plot.margin = unit(c(0,0,0,0), "cm")),
          
                  ncol = 4, nrow = 4,
                  align = c("hv"))


multiPlot2 <- plot_grid(get(acto_name.df$acto_name.df[17])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,-1,0,0), "cm")),
                        get(acto_name.df$acto_name.df[18])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,-1,0,0), "cm")),
                        get(acto_name.df$acto_name.df[19])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,0,0,0), "cm")),
                        get(acto_name.df$acto_name.df[20])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,0,0,0), "cm")),
                        
                        get(acto_name.df$acto_name.df[21])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,-1,0,0), "cm")),
                        get(acto_name.df$acto_name.df[22])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,-1,0,0), "cm")),
                        get(acto_name.df$acto_name.df[23])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,0,0,0), "cm")),
                        get(acto_name.df$acto_name.df[24])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,0,0,0), "cm")),
                        
                        get(acto_name.df$acto_name.df[25])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,-1,0,0), "cm")),
                        get(acto_name.df$acto_name.df[26])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,-1,0,0), "cm")),
                        get(acto_name.df$acto_name.df[27])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,0,0,0), "cm")),
                        get(acto_name.df$acto_name.df[28])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,0,0,0), "cm")),
                        
                        get(acto_name.df$acto_name.df[29])+ theme(plot.margin = unit(c(0,-1,0,0), "cm")),
                        get(acto_name.df$acto_name.df[30])+ theme(plot.margin = unit(c(0,-1,0,0), "cm")),
                        get(acto_name.df$acto_name.df[31])+ theme(plot.margin = unit(c(0,0,0,0), "cm")),
                        get(acto_name.df$acto_name.df[32])+ theme(plot.margin = unit(c(0,0,0,0), "cm")),
                        
                        ncol = 4, nrow = 4,
                        align = c("hv"))

multiPlot3 <- plot_grid(get(acto_name.df$acto_name.df[33])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,-1,0,0), "cm")),
                        get(acto_name.df$acto_name.df[34])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,-1,0,0), "cm")),
                        get(acto_name.df$acto_name.df[35])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,0,0,0), "cm")),
                        get(acto_name.df$acto_name.df[36])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,0,0,0), "cm")),
                        
                        get(acto_name.df$acto_name.df[37])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,-1,0,0), "cm")),
                        get(acto_name.df$acto_name.df[38])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,-1,0,0), "cm")),
                        get(acto_name.df$acto_name.df[39])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,0,0,0), "cm")),
                        get(acto_name.df$acto_name.df[40])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,0,0,0), "cm")),
                        
                        get(acto_name.df$acto_name.df[41])+ theme(plot.margin = unit(c(0,-1,0,0), "cm")),
                        get(acto_name.df$acto_name.df[42])+ theme(plot.margin = unit(c(0,-1,0,0), "cm")),
                        get(acto_name.df$acto_name.df[43])+ theme(plot.margin = unit(c(0,0,0,0), "cm")),
                        get(acto_name.df$acto_name.df[44])+ theme(plot.margin = unit(c(0,0,0,0), "cm")),
                        
                        ncol = 4, nrow = 4,
                        align = c("hv"))

png("Figures/Figure SX_actograms_1.png", width = 9, height = 14, units = "in", res = 600)
multiPlot1
dev.off()

png("Figures/Figure SX_actograms_2.png", width = 9, height = 14, units = "in", res = 600)
multiPlot2
dev.off()

png("Figures/Figure SX_actograms_3.png", width = 9, height = 14, units = "in", res = 600)
multiPlot3
dev.off()



###### PDF version -------------------------------------------------------------

multiPlot1.pdf <- plot_grid(get(acto_name.df$acto_name.df[1])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,-1,0,0), "cm")),
                        get(acto_name.df$acto_name.df[2])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,-1,0,0), "cm")),
                        get(acto_name.df$acto_name.df[3])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,0,0,0), "cm")),
                        get(acto_name.df$acto_name.df[4])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,0,0,0), "cm")),
                        
                        get(acto_name.df$acto_name.df[5])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,-1,0,0), "cm")),
                        get(acto_name.df$acto_name.df[6])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,-1,0,0), "cm")),
                        get(acto_name.df$acto_name.df[7])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,0,0,0), "cm")),
                        get(acto_name.df$acto_name.df[8])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,0,0,0), "cm")),
                        
                        get(acto_name.df$acto_name.df[9])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,-1,0,0), "cm")),
                        get(acto_name.df$acto_name.df[10])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,-1,0,0), "cm")),
                        get(acto_name.df$acto_name.df[11])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,0,0,0), "cm")),
                        get(acto_name.df$acto_name.df[12])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,0,0,0), "cm")),
                        
                        get(acto_name.df$acto_name.df[13])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,-1,0,0), "cm")),
                        get(acto_name.df$acto_name.df[14])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,-1,0,0), "cm")),
                        get(acto_name.df$acto_name.df[15])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,0,0,0), "cm")),
                        get(acto_name.df$acto_name.df[16])+ theme(axis.title.x = element_blank(), plot.margin = unit(c(0,0,0,0), "cm")),
                        
                        get(acto_name.df$acto_name.df[17])+ theme(plot.margin = unit(c(0,-1,0,0), "cm")),
                        get(acto_name.df$acto_name.df[18])+ theme(plot.margin = unit(c(0,-1,0,0), "cm")),
                        get(acto_name.df$acto_name.df[19])+ theme(plot.margin = unit(c(0,0,0,0), "cm")),
                        get(acto_name.df$acto_name.df[20])+ theme(plot.margin = unit(c(0,0,0,0), "cm")),
                        
                        
                        ncol = 4, nrow = 5,
                        align = c("hv"))


pdf("../Manuscript/SI_actogram_1.pdf")
multiPlot1.pdf
dev.off()

png("Figures/Figure SX_actograms_2.png", width = 9, height = 14, units = "in", res = 600)
multiPlot2
dev.off()

png("Figures/Figure SX_actograms_3.png", width = 9, height = 10.4, units = "in", res = 600)
multiPlot3
dev.off()


#### Plot all birds ------------------------------------------------------------

for(i in 1:length(Acto_list)){
  breaks <- c(as.character(unique(ceiling(Acto_list[[i]]$Rday/30.5)*30.5-15)))
  
  df.c <- Acto_list[[i]]
  df.c$Rdayc <- as.factor(df.c$Rday)
  df.c$month <- format(df.c$date,"%m")
  
  label <- c(ifelse(unique(ceiling(df.c$Rday/30.5))<13, unique(ceiling(df.c$Rday/30.5)), unique(ceiling(df.c$Rday/30.5))-12))
  label <- c(month.abb[(ifelse(label<13, label, label-12))])
  
  actogram <- ggplot()+
    geom_tile(aes(x = hour, y = Rdayc, fill = mean.act), data=df.c) + # colors tiles based on (mean) hourly activity. Must make column for hourly activity before running script
    scale_fill_gradientn(colours = c("white", "#009688", "#00796B", "#004d40", "black")) +
    geom_tile(data = df.c[is.na(df.c$Activity_NAs),], aes(x = hour, y = Rdayc), fill = "firebrick3", group = 1)  + #  colors NAs red
    coord_fixed(0.3) + # this adjusts the width of saved figure (smaller = wider)
    scale_x_continuous(name = "Time of day (UTC)", expand = c(0, 0), limits=c(-0.001,47.001),
                       breaks = c( 6, 18, 30, 42), labels = c("06", "18", "06", "18")) +
    scale_y_discrete(name = NULL, expand=c(0,0), limits=rev, breaks = breaks, labels=label)+
    theme_bw(base_size = 12)+
  theme(panel.border = element_rect(fill = NA), legend.key.size = unit(4, "mm"), 
          plot.title = element_text(hjust = 0.5,), legend.position = "right") +
    labs(title = Acto_list[[i]]$ID[1], fill = "Mean\nactivity", size = 12, family = "Arial")
  
  mypath <- paste0("Figures/Actograms/Act_", df.c$lat[1], "_",  df.c$colony[1], "_", names(Acto_list[i]), ".pdf")
  
  # png(filename = mypath, width = 7, height = 9, pointsize = 8, units = "in", res = 300)
  # actogram
  # dev.off()

  ggsave(
    filename = mypath,
    plot = actogram,
    device = NULL,
    path = NULL,
    scale = 1,
    width = 90,
    height = 84,
    units = c("mm"),
    dpi = 300, 
    bg = NULL)
}
