# -------------------------------------
# Script: 4_randomise-visits
# Author: Dr Natasha Gillies
# Purpose: Randomise starts and durations of foraging trips to see if random
# attendance matches observed
# Notes:
# Date: 2022-11-14
# -------------------------------------

### 0.0 Load functions & packages ----------------------------------------------

# 0.0.0 Functions

# 0.0.1 Define packages
packages <- c("dplyr", "lubridate", "lomb")

# 0.0.2 Install packages not yet installed - change lib to library path
#installed_packages <- packages %in% rownames(installed.packages())

#if (any(installed_packages == FALSE)) {
#  install.packages(packages[!installed_packages], lib = "C:/Users/libraryPath")
#}

# 0.0.3 Load packages and functions
invisible(lapply(packages, library, character.only = TRUE))
source("BLKI_functions.R")

# 0.04 Set for loop parameters
niter <- 1000 # set randomisation iterations
period_out <- vector(mode = "list", length = niter)

# 0.0.5 Load GLS data and split by ring
load("Data_inputs/BLKI_gls-trips.RData") 
glsTrips_split <- split(glsTrips, glsTrips$Ring)

### FOR LOOOP START 

### PROBLEM - PROPORTIONAL ATTENDANCE LOOKS WAY TOO LOW

for (i in 1:niter) {

  print(paste0("Commencing iteration ", i))
  
    ### 1.0 Load and randomise GLS data by ring --------------------------------------------
    
    print("Randomising data...")
  
    # 1.0.2 Use randomised GLS function to randomise data
    randomised_GLS <- lapply(glsTrips_split,
                      function(x) randomise_GLS(x$Ring[1], x$start, x$duration_mins) )
    
    # 1.0.3 Bind into single dataset
    r_glsTrips <- do.call("rbind", randomised_GLS)
  
    
    ### 2.0 Find number minutes attendance at colony by hour -----------------------
    
    print("Finding attendance...")
    
    # 2.0.0 Format datetime variables
    r_glsTrips <- rename(r_glsTrips, start_time = start, end_time = end)
    r_glsTrips$start_time <- as.POSIXct(r_glsTrips$start_time, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    r_glsTrips$end_time <- as.POSIXct(r_glsTrips$end_time, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    
    # 2.0.1 Find attendance per hour
    r_kit.att <- att_by_hour(r_glsTrips)
    
    
    ### 3.0 Run periodogram analysis -----------------------------------------------
    
    # 3.0.0 Split data by individuals 
    r_kits.ind <- r_kit.att[c("ring", "att_time", "time_ID")]
    r_kits.ind <- split(r_kits.ind, r_kits.ind$ring)
    
    print("Calculating period length...")
    
    # 3.0.1 Run the Lomb-Scargle function for each individual
    LS.ind <- lapply(r_kits.ind, function(x) lsp(x$att_time, times = x$time_ID,
                                               alpha = 0.01, type = c("period"), 
                                               ofac = 5, plot = F))

    
    ## 3.1 Summarize the LS output for each individual -----------------------------
    
    print("Summarising output...")
    
    LS.summary_list <- lapply(LS.ind,function(x)as.data.frame(t(summary.lsp(x))))
    LS.summary <- do.call("rbind", LS.summary_list)
    
    # 3.1.1 Get mean and standard error of individual period
    LSmean <- mean(as.numeric(LS.summary$`At  period`)) # 89.64485
    LSsd <- sd(as.numeric(LS.summary$`At  period`)) # 37.35065
    
    # 3.1.2 Output result
    period_out[[i]] <- cbind(i, LSmean, LSsd)

}

### FOR LOOP END


### 4.0 Output results ---------------------------------------------------------

r_results <- do.call("rbind", period_out)
save(r_results, file = "Data_outputs/randomised_periods.RData")

r_results <- data.frame(r_results)

r_LSmean <- mean(r_results$LSmean) # 106.6065

n <- nrow(r_results)
margin <- qt(0.975, df = n - 1)*LSsd/sqrt(n)

lower <- r_LSmean - margin # 104.2887
upper <- r_LSmean + margin # 108.9242

range(r_results$LSmean) # 71.8309 182.8184

# Real values : 
#  LSmean     LSsd
#  89.64485 37.35065

# Proportion period lengths greater than observed
nrow(subset(r_results, LSmean > 89.64485))/1000
nrow(subset(r_results, LSmean < 89.64485))/1000

# 4.0.1 Plot outputs

results_hist <- ggplot(aes(x = LSmean), data = r_results) + 
  geom_histogram(col = "black", fill = "grey") + 
  geom_vline(xintercept = 89.64485, 
             col = "red", linetype = "dashed", size = 1) +
  labs(y = "Count", x = "Mean period length (hours)") +
  BLKI_theme

png(paste0("Figures/randomised_periods.png"), 
    width = 9, height = 8, units = "in", res = 300)
results_hist
dev.off()
