# -------------------------------------
# Script: 3b_statistical_analysis-randomised-GAM.R
# Author: Dr Natasha Gillies
# Purpose: Randomise data used in the GAM of immersion against sun elevation
# Notes:
# Date: 2023-09-20
# -------------------------------------

# Load functions & packages ----------------------------------------------------

# Functions
source("BLKI_functions.R")

# Packages
packages <- c("dplyr")

# Install packages not yet installed - change lib to library path
#installed_packages <- packages %in% rownames(installed.packages())

#if (any(installed_packages == FALSE)) {
#  install.packages(packages[!installed_packages], lib = "C:/Users/libraryPath")
#}

# Load packages
invisible(lapply(packages, library, character.only = TRUE))


# Repeatedly resample dataframe and fit model ----------------------------------

### This is a slow option using a for loop; in reality this has been run uisng
### the University of Liverpool Condor high-throughput system.

## Load the dataset
load("Data_inputs/BLKI_attendance_2013-2022.RData")
gam_data <- kits.att_all %>% select(ring, datetime, total_imm, sun_altitude, year)

# Set structure
gam_data$year <- as.factor(gam_data$year)
gam_data$ring <- as.factor(gam_data$ring)
gam_data <- na.omit(gam_data)


## Set loop parameters
n_iter <- 10
output_list <- vector(mode = "list", length = n_iter)

for (i in 1:n_iter) {
  
    tok <- Sys.time()
  
    ## Subset the data to random 72 hour samples to make manageable
    random_sample <- subsample_by_hours(gam_data, 72)
    
    random_sample[,c(1, 5)] <- lapply(random_sample[,c(1, 5)], as.factor)
   
    ## Fit the model
    gam_model <- mgcv::gam(total_imm ~ s(sun_altitude) +
                          s(ring, bs = 're'), 
                           data = random_sample)
    
    # Extract relevant parameters
    gam_summary <- summary(gam_model)
    sun_edf <- gam_summary$edf[1]
    sun_p <- gam_summary$s.table[1,4]
    
    summary.df <- data.frame(iter = i,
                             sun.edf = sun_edf, 
                             sun.p = sun_p)
    
    # Output to list
    output_list[[i]] <- summary.df
    
    tik <- Sys.time() 
    
    print(tik-tok)
    
}

## Examine results -------------------------------------------------------------

random_gam <- do.call("rbind", output_list)
save(random_gam, file = "Data_outputs/randomised_GAM_results.RData")



# Examine results from high-throughput system ==================================

# Point to results folder 
results_folder <- "./Condor/results/"
filenames <- list.files(path = results_folder, pattern = "*.RData")

results_list <- vector("list", length(filenames))

# Loop through each and extract results

for (i in 1:length(filenames)) {
  print(i)
  load(paste0(results_folder, filenames[i]))
  summary.df$iter <- i
  results_list[[i]] <- summary.df
  
}

resultsTable <- do.call("rbind", results_list)  


hist(resultsTable$sun.edf)
mean(resultsTable$sun.edf)
# [1] 0.06358375
sd(resultsTable$sun.edf)
# [1] 8.825313
