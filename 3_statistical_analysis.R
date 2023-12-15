## ---------------------------
##
## Script name: 3_statistical_analysis.R
##
## Purpose of script: Analyse nest visitation using attendance datasets
##
## Author: Dr. Natasha Gillies
##
## Created: 2022-10-03
##
## Email: gilliesne@gmail.com
##
## Methods guided by those outlined in : 
## Huffeldt Nicholas Per and Merkel Flemming R. 2016 Sex-specific, inverted 
## rhythms of breeding-site attendance in an Arctic seabird. Biol. Lett. 122016028920160289
## 
## ---------------------------


# Preamble =====================================================================

# Define the packages
packages <- c("zoo", "stringr", "ggplot2", "binom", "grid", "mgcv", "lomb", 
              "lubridate", "purrr", "dplyr","ggpubr", "magrittr", "lme4",
              "reshape2", "sjPlot", "glmmTMB", "DHARMa", "gamlss", "broom.mixed",
              "ggeffects")

# Install packages not yet installed - change lib to library path
# installed_packages <- packages %in% rownames(installed.packages())
# 
#  if (any(installed_packages == FALSE)) {
#   install.packages(packages[!installed_packages])
#  }

# Load packages and functions
invisible(lapply(packages, library, character.only = TRUE))

source("BLKI_functions.R")

# Create data_outputs folder if doesn't already exist
out.path <- "./Data_outputs/"

# if(dir.exists(out.path) == FALSE) {
#   dir.create(out.path)
# }


# Load & process data ==========================================================
load("Data_inputs/BLKI_attendance_2013-2022.RData")

# Quick summary stats
kits.att_all %>% 
  summarise(n_colonies = n_distinct(colony),
            n_rings = n_distinct(ring),
            year_min = min(year), 
            year_max = max(year))

# n_colonies n_rings year_min year_max
# 23     962     2013     2022


### Extract mean hourly attendance, mean visits, mean immersion ----------------
# att_by_colony <- kits.att_all %>% 
#   group_by(year, colony, col_lat, col_lon, datetime) %>%
#   summarise(hourly_att = mean(att_time/60),
#             hourly_vis = sum(num_visits)/n_distinct(ring),
#             hourly_imm = mean(total_imm)) %>%
#   mutate(hour = hour(datetime))

#save(att_by_colony, file = "Data_outputs/attendance_by_colony.Rdata")
load("Data_outputs/attendance_by_colony.Rdata")

#### Summarise mean and CIs for each year, colony, hour -------------------------
# summary_colony_att <- att_by_colony %>%
#   group_by(year, colony, col_lat, col_lon, hour) %>%
#   summarise(mean_att = mean(hourly_att),
#             lower_att = confint(lm(hourly_att~1))[1],
#             upper_att = confint(lm(hourly_att~1))[2],
#             mean_vis = mean(hourly_vis),
#             lower_vis = confint(lm(hourly_vis~1))[1],
#             upper_vis = confint(lm(hourly_vis~1))[2],
#             mean_imm = mean(hourly_imm),
#             lower_imm = confint(lm(hourly_imm~1))[1],
#             upper_imm = confint(lm(hourly_imm~1))[2])

#save(summary_colony_att, file = "Data_outputs/summarised_colony_att.Rdata")
load("Data_outputs/summarised_colony_att.Rdata")



### Get immersion for all birds per time at each colony ------------------------
# kits.pop.imm <- kits.att_all %>% 
#   group_by(colony, col_lat, col_lon, datetime) %>%
#   summarise(tot.imm = sum(total_imm),
#             tot.kits = n_distinct(ring),
#             prop.imm = tot.imm/tot.kits) %>%
#   data.frame()

#save(kits.pop.imm, file = "Data_outputs/BLKI_pop_imm.RData")
load("Data_outputs/BLKI_pop_imm.RData")


### Load behavioural data ------------------------------------------------------
load("Data_inputs/BLKI_gls-daily-behaviour.RData")

# Isolate summer months
gls_daily %<>% mutate(month = as.numeric(lubridate::month(date))) %>% 
  filter(month >= 6 & month <= 8) %>% dplyr::select(-month)

# Get year
gls_daily %<>% mutate(year = format(as.Date(date, format = "%Y-%m-%d"),"%Y"))



### Colony attendance data -----------------------------------------------------

# Prepare data
gls_daily$col_mins <- gls_daily$col_att24*gls_daily$daylight.mins

col_behav <- gls_daily %>%
  dplyr::select(ring, date, col_att.day, n_recs.day, col_lat, year, col_mins) %>% 
  filter(!is.na(n_recs.day))


# Visualise population-level behaviour (attendance/visits/immersion) ===========

## ~ FIGURES ~ Plot attendance, visits, immersion for each year by colony ------------------------------

# Create plots by colony, and output to png
colonies <- unique(summary_colony_att$colony)

for (i in 1:length(colonies)) {
  
  mean_att.plot <- subset(att_by_colony, colony == colonies[i])
  mean_att.plot$date <- as.Date(mean_att.plot$datetime)
  summary_att.plot <- subset(summary_colony_att, colony == colonies[i])
  
  # Attendance plot
  att_plot <- ggplot() +
    geom_point(data = mean_att.plot,
               aes(x = hour, y = hourly_att, colour = as.factor(date))) +
    theme(legend.position = "none") +
    geom_smooth(data = summary_att.plot,
                aes(x = hour, y = mean_att, colour = year, group = year,
                    ymin = lower_att, ymax = upper_att),
                method = "gam", formula = y~s(x, k = 24), se = T,
                linetype = 2, col = "grey25", alpha = 0.25, lwd = 0.4) +
    labs(x = "Hour of day", y = "Mean attendance time (proportion of hour)") +
    ggtitle(paste0(colonies[i], ", ", mean_att.plot$col_lat[1])) +
    scale_x_continuous(breaks = seq(0, 24, 4), limits = c(0, 23)) +
    scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1)) +
    BLKI_theme +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "none")
  
  # Visits plot
  vis_plot <- ggplot() +
    geom_point(data = mean_att.plot,
               aes(x = hour, y = hourly_vis, colour = as.factor(date))) +
    theme(legend.position = "none") +
    geom_smooth(data = summary_att.plot,
                aes(x = hour, y = mean_vis, colour = year, group = year,
                    ymin = lower_vis, ymax = upper_vis),
                method = "gam", formula = y~s(x, k = 24), se = T,
                linetype = 2, col = "grey25", alpha = 0.25, lwd = 0.4) +
    labs(x = "Hour of day", y = "Mean visitation rate") +
    BLKI_theme +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "none")
  
  
  # Immersion plot
  imm_plot <- ggplot() +
    geom_point(data = mean_att.plot,
               aes(x = hour, y = hourly_imm, colour = as.factor(date))) +
    theme(legend.position = "none") +
    geom_smooth(data = summary_att.plot,
                aes(x = hour, y = mean_imm, colour = year, group = year,
                    ymin = lower_imm, ymax = upper_imm),
                method = "gam", formula = y~s(x, k = 24), se = T,
                linetype = 2, col = "grey25", alpha = 0.25, lwd = 0.4) +
    labs(x = "Hour of day", y = "Mean hourly immersion") +
    scale_x_continuous(breaks = seq(0, 24, 4), limits = c(0, 23)) +
    BLKI_theme +
    theme(legend.position = "none")
  
  
  # Output the plot
  final_plot <-  ggarrange(att_plot, vis_plot, imm_plot,
                           ncol = 1,
                           nrow = 3)
  
  png(paste0("Figures/Colony attendance/", colonies[i],"_attendance.png"), 
                    width = 12, height = 24, units = "in", res = 300)
  print(final_plot)
  dev.off()
 
}

beepr::beep(8)


## ~ FIGURES ~ Plot immersion against environment metrics ----------------------

colonies <- unique(summary_colony_att$colony)

for (i in 1:length(colonies)) {
  
  print(paste0("Processing colony ", i, " of ", length(colonies)))
  
  mean_att.plot <- subset(att_by_colony, colony == colonies[i])
  mean_att.plot$date <- as.Date(mean_att.plot$datetime)
  summary_att.plot <- subset(summary_colony_att, colony == colonies[i])

  # Calculate mean sun elevation per hour, temperature per hour and year
  env_summary <- kits.att_all %>%
                  filter(colony == colonies[i]) %>% 
                  distinct(datetime, .keep_all = T) %>%
                  group_by(hour) %>%
                  summarise(mean_alt = mean(sun_altitude),
                            mean_temp = mean(temp, na.rm = T),
                            sd_temp = sd(temp, na.rm = T),
                            time = format(datetime[1], "%H%:%M:%S")) %>%
                  mutate(hour = as.numeric(hour))
  
    ### (A) Plot sun elevation angle by hour -----------------------------------
    sun_alt_plot <- ggplot(data = env_summary, aes(x = hour, y = mean_alt)) +
      geom_point(size = 5, shape = "\u263C", col = "#FFE133") +
      geom_line(stat = "smooth", method = "gam", formula = y~s(x, k = 24), se = F,
                linetype = 2, col = "grey25", alpha = 0.5) +
      labs(y = expression("Sun elevation angle ("*~degree*")")) +
      scale_x_continuous(breaks = seq(0, 24, 4), limits = c(0, 23)) +
      scale_y_continuous(breaks = seq(-30, 70, 30), limits = c(-30, 70)) +
      geom_hline(yintercept = 0, linetype = 2) +
      ggtitle(paste0(colonies[i], ", ", mean_att.plot$col_lat[1])) +
      BLKI_theme +
      theme(legend.position = "none",
            axis.text.x = element_blank(),
            axis.title.x = element_blank())
  

    ## (B) Plot average hourly temperature fluctuation (across years ) ---------
    temp_plot <- ggplot(env_summary, aes(x = hour, y = mean_temp)) +
      geom_line(group = 1) +
      geom_ribbon(aes(ymin = mean_temp - sd_temp,
                      ymax = mean_temp + sd_temp,
                      x = hour),
                  fill = "grey25", alpha = 0.2) +
      labs(y = expression("Mean Temperature ("*~degree*C*")")) +
      scale_x_continuous(breaks = seq(0, 24, 4), limits = c(0, 23)) +
      BLKI_theme +
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank())


  ## (C) Plot overall average immersion ----------------------------------------
  
  imm_plot <- ggplot() +
    geom_point(data = mean_att.plot,
               aes(x = hour, y = hourly_imm, colour = as.factor(date))) +
    theme(legend.position = "none") +
    geom_smooth(data = summary_att.plot,
                aes(x = hour, y = mean_imm, colour = year, group = year,
                    ymin = lower_imm, ymax = upper_imm),
                method = "gam", formula = y~s(x, k = 24), se = T,
                linetype = 2, col = "grey25", alpha = 0.25, lwd = 0.4) +
    labs(x = "Hour of day", y = "Mean hourly immersion") +
    scale_x_continuous(breaks = seq(0, 24, 4), limits = c(0, 23)) +
    BLKI_theme +
    theme(legend.position = "none")
  
  ## (D) Output the combined plot per colony -------------------------------------
  
  combined_plot <- ggarrange(sun_alt_plot, temp_plot, imm_plot, 
                             ncol = 1, align = "v",
                             heights = c(1, 1, 2))
  
  png(paste0("Figures/Colony attendance/", colonies[i], "_env_imm.png"), 
      width = 12, height = 12, units = "in", res = 300)
  print(combined_plot)
  dev.off()
  
}



# Find population-level ACF ====================================================

## ~ FIGURE 2A ~ Plot colony prop immersion by latitude ------------------------
colonies <- unique(kits.pop.imm$colony[order(kits.pop.imm$col_lat)])
plot_list <- vector(mode = "list", length = length(colonies))

for (i in 1:length(colonies)) {
  
  print(paste0("Processing colony ", i, " of ", length(colonies)))
  
  pop.imm_plot <- kits.pop.imm %>% 
    filter(colony == colonies[[i]]) %>% 
    mutate(hour = hour(datetime), date = as.Date(datetime))
  
  pop.imm_figure <- ggplot() +
    geom_point(data = pop.imm_plot,
               aes(x = hour, y = prop.imm, colour = as.factor(date))) +
    geom_smooth(data = pop.imm_plot,
                aes(x = hour, y = prop.imm),
                method = "gam") +
    theme(legend.position = "none") +
    ggtitle(paste0(colonies[i], ", ", pop.imm_plot$col_lat[1])) +
    labs(x = "Hour of day", y = "Mean proportional immersion") +
    scale_x_continuous(breaks = seq(0, 24, 4), limits = c(0, 23)) +
    BLKI_theme +
    theme(legend.position = "none")
  
  plot_list[[i]] <- pop.imm_figure  
  
}

png("Figures/Figure 2A_imm_by_colony.png", 
    width = 30, height = 24, units = "in", res = 300)
gridExtra::grid.arrange(grobs = plot_list, nrow = 4)
dev.off()


## Run autocorrelatIon function on each colony ---------------------------------
colonies <- unique(kits.pop.imm$colony[order(kits.pop.imm$col_lat)])
plot_list <- vector(mode = "list", length = length(colonies))
acf_data_list <- vector(mode = "list", length = length(colonies))

for (i in 1:length(colonies)) {
  
  print(paste0("Processing colony ", i, " of ", length(colonies)))
  
  pop.imm <- subset(kits.pop.imm, colony == colonies[i])
  
  imm.ACF <- acf(pop.imm$prop.imm, na.action = na.pass, lag = 60)
 
  ## Extract the lag time and strength for the peak autocorrelation
  
  # Find the lag time corresponding to the maximum value after the first trough
  acf_values <- imm.ACF$acf
  trough_lag <- which(diff(sign(diff(acf_values))) == 2)[1] + 1 
  second_peak_lag <- trough_lag + which.max(acf_values[trough_lag:length(acf_values)])
  
  peak_strength <- imm.ACF$acf[second_peak_lag] - imm.ACF$acf[trough_lag]
  
  # Add to dataset
  acf_data_list[[i]] <- data.frame(colony = pop.imm$colony[1],
                       col_lat = pop.imm$col_lat[1],
                       col_lon = pop.imm$col_lon[1],
                       peak_lag = second_peak_lag,
                       peak_strength = peak_strength)
  
  ## Convert output to dataframe and process 
  imm.ACF_df <- as.data.frame(unlist(imm.ACF))
  imm.ACF_df <- cbind(Row.Names = rownames(imm.ACF_df), imm.ACF_df)
  
  # First 60 rows (equivalent to lag number) are the acf values
  imm.ACF_df <- as.data.frame(imm.ACF_df[1:61,])
  imm.ACF_df$Lag <- 0:60
  
  # Rename columns
  imm.ACF_df <- rename(imm.ACF_df, "ACF" = "unlist(imm.ACF)")
  imm.ACF_df$ACF <- as.numeric(imm.ACF_df$ACF)
  
  imm.ACF_df$Row.Names <- NULL

  # Calculate CIs
  n <- nrow(pop.imm)
  z <- qnorm((1 + 0.95) / 2)
  se_acf <- 1 / sqrt(n)
  
  imm.ACF_df$lower_ci <- imm.ACF_df$ACF - z * se_acf
  imm.ACF_df$upper_ci <- imm.ACF_df$ACF + z * se_acf
  
  ##  Plot attendance ACF
  imm.ACF_plot <- ggplot() + 
    geom_hline(yintercept = c(0)) +
    geom_hline(yintercept = c((1.96/sqrt(n)), -(1.96/sqrt(n))), linetype = 2) +
    geom_line(data = imm.ACF_df, aes(x = Lag, y = ACF)) +
    scale_y_continuous(limits = c(-1,1), expand = c(0,0)) +
    scale_x_continuous(breaks = c(0, 12, 24, 36, 48, 60), expand = c(0,0)) +
    BLKI_theme +
    ggtitle(paste0(colonies[i], ", ", pop.imm$col_lat[1])) +
    theme(plot.margin = margin(10, 10, 10, 10))
  
  plot_list[[i]] <- imm.ACF_plot  
  
}

acf_col <- do.call("rbind", acf_data_list)
#save(acf_col, file = "Data_outputs/acf_colony_data.RData")
load("Data_outputs/acf_colony_data.RData")


### ~ FIGURE 2B ~ ACF in immersion over all colonies ---------------------------
png("Figures/Figure 2B_colony_imm_ACF.png", 
    width = 30, height = 24, units = "in", res = 300)
gridExtra::grid.arrange(grobs = plot_list, nrow = 4)
dev.off()


## ACF peak strength and lag with colony latitude ------------------------------

peak_lag.lm <- lm(peak_lag ~ col_lat, data = acf_col)
summary(peak_lag.lm)
median(acf_col$peak_lag)
var(acf_col$peak_lag)

peak_strength.lm <- lm(peak_strength ~ col_lat, data = acf_col)
summary(peak_strength.lm)
mean(acf_col$peak_strength)

### ~ FIGURE 3 ~ Peak ACF by latitude ------------------------------------------

png("Figures/Figure 3_LM_peak_by_lat.png", 
    width = 7, height = 5, units = "in", res = 600)
ggplot(acf_col, aes(x = col_lat, y = peak_strength)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "brown2", lwd = 0.5, fill = "grey") +
  labs(x = expression("Colony latitude ("*~degree*" north)"), y = "ACF peak strength") +
  BLKI_theme
dev.off()


# Visualise individual-level rhythms ===========================================

# Find length of time series for each nest
# (must divide by 2 cause need at least 2 observations for a TS)
bird.ts <- tapply(kits.att_all$datetime, list(kits.att_all$ring), function(x)(length(x))/2)
mean(bird.ts); sd(bird.ts); qt(0.975,df=length(bird.ts)-1)*sd(bird.ts)/sqrt(length(bird.ts)) 

# [1] 2213.483
# [1] 1511.379
# [1] 95.62718


### Find autocorrelation function (ACF) for each individual --------------------

# Split the dataframe by ring so can apply function to each individual
kits.att_split <- data.frame(kits.att_all[,c("ring","total_imm")])
kits.att_split <- split(kits.att_split, kits.att_split$ring)

ind_ACF <- lapply(kits.att_split,
                  function(x) acf(x$total_imm, na.action = na.pass, lag = 60,
                                  main = x$ring[1])) 

# Create a dataframe from list output
ACF.df <- as.data.frame(unlist(ind_ACF))
ACF.df <- cbind(Row.Names = rownames(ACF.df), ACF.df)

## Extract rings 
ACF.df$ring <- substr(as.character(ACF.df$Row.Names),1,nchar(as.character(ACF.df$Row.Names))-5)
ACF.df$ring <- gsub('[[:punct:]]', '', ACF.df$ring)
ACF.df$ring <- as.factor(ACF.df$ring)

## Find lag ID
ACF.df$Lag <- substr(as.character(ACF.df$Row.Names),8,nchar(as.character(ACF.df$Row.Names)))
ACF.df$Lag <- gsub('[[:punct:]]', '', ACF.df$Lag)

## Isolate acf values
ACF.df$ACF_ID <- str_sub(ACF.df$Lag, 1,3)

ACF.df$Lag <- substring(ACF.df$Lag, 4)
ACF.df$Lag <- as.numeric(ACF.df$Lag)

ACF.df <- rename(ACF.df, ACF = "unlist(ind_ACF)")
ACF.df$ACF <- as.numeric(ACF.df$ACF)

ACF.df$Row.Names <- NULL

## Isolate acf values
ACF.df <- subset(ACF.df, ACF_ID == "acf")

summary(ACF.df)
length(unique(ACF.df$ring))

## Calculate CIs for plotting
n <- nrow(kits.att_all)
z <- qnorm((1 + 0.95) / 2)
se_acf <- 1 / sqrt(n)

ACF.df$lower_ci <- ACF.df$ACF - z * se_acf
ACF.df$upper_ci <- ACF.df$ACF + z * se_acf

# Merge in colony info
colonies <- kits.att_all %>% dplyr::select(ring, colony, col_lat) %>% distinct(ring, .keep_all = TRUE) 
ACF.df <- merge(ACF.df, colonies, by = "ring")


#### Plot ACF for each individual by colony ------------------------------------

colonies <- unique(ACF.df$colony[order(ACF.df$col_lat)])
plot_list <- vector(mode = "list", length = length(colonies))

for (i in 1:length(colonies)) {
  
  print(paste0("Processing colony ", i, " of ", length(colonies)))
  
  acf_plot <- subset(ACF.df, colony == colonies[i])
  
  ACF.ind_plot <- ggplot(acf_plot, aes(x = Lag, y = ACF, col = ring)) + 
    geom_line() +
    geom_hline(yintercept = c(0)) +
    geom_line(alpha = 0.5) +
    scale_y_continuous(breaks = seq(-1, 1, 0.5), limits = c(-1,1), expand = c(0,0)) +
    scale_x_continuous(breaks = seq(0, 60, 12), limits = c(0, 60), expand = c(0,0)) +
    BLKI_theme +
    ggtitle(paste0(colonies[i], ", ", acf_plot$col_lat[1])) +
    theme(legend.position = "none",
          plot.margin = margin(10, 10, 10, 10))
  
  plot_list[[i]] <- ACF.ind_plot  
  
}

### ~ FIGURE S4 ~ ACF in immersion by individuals over colonies ----------------
png("Figures/Figure S4_individual_imm_ACF.png", 
    width = 30, height = 24, units = "in", res = 300)
gridExtra::grid.arrange(grobs = plot_list, nrow = 4)
dev.off()


# GAM: immersion values ~ sun elevation ========================================

## Subset the data to random 72 hour samples to make manageable
set.seed(817)
random_sample <- subsample_by_hours(kits.att_all, 72)

random_sample[,c(3, 7)] <- lapply(random_sample[,c(3, 7)], as.factor)
random_sample$ring <- as.factor(random_sample$ring)

#save(random_sample, file = "Data_outputs/subsampled_kits.RData")
load("Data_outputs/subsampled_kits.RData")

random_sample$date <- as.Date(random_sample$datetime)


### Fit the model --------------------------------------------------------------

gam_model <- mgcv::gam(total_imm ~ s(sun_altitude) +
                   s(ring, bs = 're'), 
                 data = random_sample)

#save(gam_model, file = "Data_outputs/GAM_immersion_sun.RData")
load("Data_outputs/To delete/GAM_immersion_sun.RData")

gam.check(gam_model)
summary(gam_model)



### ~ FIGURE 4 ~ Immersion against sun elevation angle -------------------------

# Process data for plotting
p_obj <- plot(gam_model, residuals = TRUE)
p_obj <- p_obj[[1]] # just one smooth so select the first component
sm_df <- as.data.frame(p_obj[c("x", "se", "fit")])
sm_df$fit <- sm_df$fit + coef(gam_model)[1]
data_df <- as.data.frame(p_obj[c("raw", "p.resid")])

# Build the plot
gam_immersion.plot <- ggplot(sm_df, aes(x = x, y = fit)) +
  geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
              alpha = 0.25, col = "grey", linetype = 2) +
  geom_line(col = "brown2", lwd = 0.5) +
  labs(x = expression("Sun elevation angle ("*~degree*")"),
       y = "Predicted immersion") +
  BLKI_theme 


png("Figures/Figure 4_GAM_imm_by_sun.png", 
    width = 7, height = 5, units = "in", res = 600)
gam_immersion.plot
dev.off()

# Assess effects on behaviour ==================================================

### Look for autocorrelation in foraging effort --------------------------------

# Remove examples where < 14 days of data
gls_daily <- subset(gls_daily, !is.na(prop_forage.day))

acf_dat <- gls_daily %>% 
                  group_by(ring, year) %>% 
                  filter(n_distinct(date) >= 14) %>%
            mutate(ringYr = paste(ring, year, sep = "_")) %>%
            select(c(ringYr, prop_forage.day, date)) %>%
            data.frame()

acfs <- sapply(split(acf_dat$prop_forage.day, acf_dat$ringYr), 
               function (x) c(acf(x, lag.max = 10, plot = FALSE)$acf)[-1])
boxplot(t(acfs))

acfs <- data.frame(acfs)
id <- seq.int(1:nrow(acfs))
acfs <- cbind(id, acfs)

# Get mean and SD for each acf lag
acfs$row_std = apply(acfs[,-1], 1, sd)
acfs$row_mean = apply(acfs[,-1], 1, mean)

acfs_summary <- acfs[,c("id", "row_std", "row_mean"),]

# Plot autocorrelation
acf_forage_plot <- ggplot() + 
  geom_point(data = acfs_summary, aes(x = id, y = row_mean)) + 
  geom_errorbar(data =acfs_summary, aes(ymin = row_mean - row_std, 
                                  ymax = row_mean + row_std,
                                  x = id)) +
  ylim(c(-0.2,1)) + 
  scale_x_continuous(breaks = c(1:10)) +
  labs(y = "Mean ACF", x = "Lag") + 
  BLKI_theme


### ~ FIGURE S3 ~ Foraging autocorrelaiton -------------------------------------
png("Figures/Figure S3_acf_forage_plot.png", width = 9, height = 7, units = "in", res = 300)
acf_forage_plot
dev.off()


## Model effects on at-sea behaviour ===========================================

### Construct models -----------------------------------------------------------

# Remove full day colony visits
multi_behav <- gls_daily %>% filter(col_att.day < 0.95 & prop_flight.day < 0.95 & prop_rest.day < 0.95 & prop_forage.day < 0.95)

## Construct models for rest/flight/forage separately - proportion
forageProp_glmm <- glmmTMB(prop_forage.day ~ col_lat + (1|year) + (1|ring), 
                    family = beta_family(link = "logit"),
                    ziformula = ~.,
                    data = multi_behav)

restProp_glmm <- glmmTMB(prop_rest.day ~ col_lat + (1|year) + (1|ring), 
                       family = beta_family(link = "logit"),
                       ziformula = ~.,
                       data = multi_behav)

flightProp_glmm <- glmmTMB(prop_flight.day ~ col_lat + (1|year) + (1|ring), 
                       family = beta_family(link = "logit"),
                       ziformula = ~.,
                       data = multi_behav)

# save(forageProp_glmm, file = "Data_outputs/forageProp_glmm.RData")
# save(restProp_glmm, file = "Data_outputs/restProp_glmm.RData")
# save(flightProp_glmm, file = "Data_outputs/flightProp_glmm.RData")

load("Data_outputs/forageProp_glmm.RData")
load("Data_outputs/restProp_glmm.RData")
load("Data_outputs/flightProp_glmm.RData")

## Construct models for rest/flight/forage separately - mins

# Calculate mins
multi_behav$daylight.mins[multi_behav$daylight.mins == 0] <- 1440
multi_behav$forage_mins <- multi_behav$prop_forage.day * multi_behav$daylight.mins
multi_behav$rest_mins <- multi_behav$prop_rest.day * multi_behav$daylight.mins
multi_behav$flight_mins <- multi_behav$prop_flight.day * multi_behav$daylight.mins

# Fit models
forageMins_lmm <- lmer(forage_mins ~ col_lat + (1|ring) + (1|year),
                           data = multi_behav)
forageMins_lmm.null <-  lmer(forage_mins ~ 1 + (1|ring) + (1|year),
                             data = multi_behav)
anova(forageMins_lmm, forageMins_lmm.null)


restMins_lmm <- lmer(rest_mins ~ col_lat + (1|ring) + (1|year),
                       data = multi_behav)
restMins_lmm.null <-  lmer(rest_mins ~ 1 + (1|ring) + (1|year),
                             data = multi_behav)
anova(restMins_lmm, restMins_lmm.null)


flightMins_lmm <- lmer(flight_mins ~ col_lat + (1|ring) + (1|year),
                     data = multi_behav)
flightMins_lmm.null <-  lmer(flight_mins ~ 1 + (1|ring) + (1|year),
                           data = multi_behav)
anova(flightMins_lmm, flightMins_lmm.null)


# save(forageMins_lmm, file = "Data_outputs/forageMins_lmm.RData")
# save(restMins_lmm, file = "Data_outputs/restMins_lmm.RData")
# save(flightMins_lmm, file = "Data_outputs/flightMins_lmm.RData")

load("Data_outputs/forageMins_lmm.RData")
load("Data_outputs/restMins_lmm.RData")
load("Data_outputs/flightMins_lmm.RData")

#### Explore results -----------------------------------------------------------

latitudes <- seq(45, 80, by = 5)

## Proportional results
# Change coefficients to $cond = conditional; $zi = zero inflated coefficients

# Foraging #
cond_coefs.forProp <-  summary(forageProp_glmm)$coefficients$zi
confint(forageProp_glmm)

intercept <-  cond_coefs.forProp[1,1]
beta_collat <- cond_coefs.forProp[2,1]

means_cond_forage <- data.frame(latitude = latitudes)
means_cond_forage$prob <- logit2prob(intercept + beta_collat * means_cond_forage$latitude)
(means_cond_forage$prob[1]*100) - (means_cond_forage$prob[2]*100)


# Rest #
cond_coefs.restProp <-  summary(restProp_glmm)$coefficients$zi
confint(restProp_glmm)

intercept <-  cond_coefs.restProp[1,1]
beta_collat <- cond_coefs.restProp[2,1]

means_cond_rest <- data.frame(latitude = latitudes)
means_cond_rest$prob <- logit2prob(intercept + beta_collat * means_cond_rest$latitude)
(means_cond_rest$prob[1]*100) - (means_cond_rest$prob[2]*100)

# Flight #
cond_coefs.flightProp <-  summary(flightProp_glmm)$coefficients$zi
confint(flightProp_glmm)

intercept <-  cond_coefs.flightProp[1,1]
beta_collat <- cond_coefs.flightProp[2,1]

means_cond_flight <- data.frame(latitude = latitudes)
means_cond_flight$prob <- logit2prob(intercept + beta_collat * means_cond_flight$latitude)
(means_cond_flight$prob[1]*100) - (means_cond_flight$prob[2]*100)



## Minute effects

# Forage #
summary(forageMins_lmm)
confint(forageMins_lmm)
anova(forageMins_lmm, forageMins_lmm.null)

# Rest #
summary(restMins_lmm)
confint(restMins_lmm)
anova(restMins_lmm, restMins_lmm.null)

# Flight #
summary(flightMins_lmm)
confint(flightMins_lmm)
anova(flightMins_lmm, flightMins_lmm.null)



### Plotting ===================================================================

# Get raw data for plotting
plot_dat <- multi_behav %>% dplyr::select(c(daylight.mins, col_lat, ring,
                                            prop_forage.day, prop_rest.day, prop_flight.day,
                                            forage_mins, rest_mins, flight_mins))

## Set to long format
### Proportion
plot_dat.prop <- melt(plot_dat, id.vars = c("col_lat", "ring", "daylight.mins"),
     measure.vars = c("prop_forage.day", "prop_rest.day", "prop_flight.day"),
     variable.name = "behaviour", value.name = "proportion")

level_map <- c("prop_forage.day" = "forage",
               "prop_rest.day" = "rest",
               "prop_flight.day" = "flight")

plot_dat.prop$behaviour <- level_map[plot_dat.prop$behaviour]
plot_dat.prop$behaviour <- factor(plot_dat.prop$behaviour, levels = c("forage", "rest", "flight"))

### Minutes
plot_dat.mins <- melt(plot_dat, id.vars = c("col_lat", "ring", "daylight.mins"),
                      measure.vars = c("forage_mins", "rest_mins", "flight_mins"),
                      variable.name = "behaviour", value.name = "minutes")

level_map <- c("forage_mins" = "forage",
               "rest_mins" = "rest",
               "flight_mins" = "flight")

plot_dat.mins$behaviour <- level_map[plot_dat.mins$behaviour]
plot_dat.mins$behaviour <- factor(plot_dat.mins$behaviour, levels = c("forage", "rest", "flight"))

# Set plotting colours
forage_col <- "#D81B60"
rest_col <- "#1E88E5"
flight_col <- "#FFC107"
colony_col <- "#004D40"


#### Plot proportional effects -------------------------------------------------

## Get effect sizes for each behaviour
effects_for.prop <- data.frame(ggeffects::ggpredict(forageProp_glmm, terms = "col_lat")) %>% 
  rename(col_lat = x, behaviour = group) %>%
  mutate(behaviour = "forage")

effects_rest.prop <- data.frame(ggeffects::ggpredict(restProp_glmm, terms = "col_lat")) %>% 
  rename(col_lat = x, behaviour = group)  %>%
  mutate(behaviour = "rest") 

effects_flight.prop <- data.frame(ggeffects::ggpredict(flightProp_glmm, terms = "col_lat")) %>% 
  rename(col_lat = x, behaviour = group)  %>%
  mutate(behaviour = "flight")

# Bind together
proportional_effects <- rbind(effects_for.prop,
                              effects_rest.prop,
                              effects_flight.prop)


#save(proportional_effects, file = "Data_outputs/proportional_plotting_data.RData")
load("Data_outputs/proportional_plotting_data.RData")

### Build the plot ##

proportional_behaviour_plot <- ggplot() + 
  geom_point(aes(x = col_lat, y = proportion, group = behaviour, col = behaviour),
             position = position_jitter(w = 0.5),
             data = plot_dat.prop[sample(nrow(plot_dat.prop), round(nrow(plot_dat.prop)/100) ),]) +
  geom_line(aes(x = col_lat, y = predicted, group = behaviour, col = behaviour), 
            data = proportional_effects) +
  geom_ribbon(aes(x = col_lat, ymin = conf.low, ymax = conf.high, group = behaviour),
              data = proportional_effects,
              alpha = 0.5, fill = "grey") +
  scale_colour_viridis_d(labels = c("Forage", "Rest", "Flight"),
                      name = "Behaviour",
                      option = "plasma") +
  labs(x = "", 
       y = "Proportion daylight spent in behaviour") +
  BLKI_theme +
  theme(legend.position = c(0.07, 0.9),
        legend.background = element_rect(fill = "transparent"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())


#### Plot raw time effects --------------------------------------------------------

## Get effect sizes for each behaviour
effects_for.mins <- data.frame(ggeffects::ggpredict(forageMins_lmm, terms = "col_lat")) %>% 
  rename(col_lat = x, behaviour = group) %>%
  mutate(behaviour = "forage")

effects_rest.mins <- data.frame(ggeffects::ggpredict(restMins_lmm, terms = "col_lat")) %>% 
  rename(col_lat = x, behaviour = group) %>%
  mutate(behaviour = "rest")

effects_flight.mins <- data.frame(ggeffects::ggpredict(flightMins_lmm, terms = "col_lat")) %>% 
  rename(col_lat = x, behaviour = group) %>%
  mutate(behaviour = "flight")

# Bind together
mins_effects <- rbind(effects_for.mins,
                      effects_rest.mins,
                      effects_flight.mins)


#save(mins_effects, file = "Data_outputs/rawMins_plotting_data.RData")
load("Data_outputs/rawMins_plotting_data.RData")


## Build the plot ##

mins_behaviour_plot <- ggplot() + 
  geom_point(aes(x = col_lat, y = minutes, group = behaviour, col = behaviour),
             position = position_jitter(w = 0.5),
             data = plot_dat.mins[sample(nrow(plot_dat.mins), round(nrow(plot_dat.mins)/100) ),]) +
  geom_line(aes(x = col_lat, y = predicted, group = behaviour, col = behaviour), 
            data = mins_effects) +
  geom_ribbon(aes(x = col_lat, ymin = conf.low, ymax = conf.high, group = behaviour),
              data = mins_effects,
              alpha = 0.5, fill = "grey") +
  scale_colour_viridis_d(labels = c("Forage", "Rest", "Flight", "Colony"),
                      name = "Behaviour",
                      option = "plasma") +
  labs(x = expression("Colony latitude ("*~degree*" north)"), 
       y = "Minutes per day spent in behaviour") +
  BLKI_theme +
  theme(legend.position = "none")



### ~ FIGURE 5 ~ At-sea raw mins behaviours as a function of latitude ==========

png("Figures/Figure 5_behaviour_by_latitude.png", width = 9, height = 14, units = "in", res = 600)
ggarrange(proportional_behaviour_plot, mins_behaviour_plot, nrow = 2, align = "hv")
dev.off()


## Model effects on colony attendance ==========================================

### Construct models -----------------------------------------------------------

# Fit proportional model
colProp_BEINF <- gamlss(col_att.day ~ col_lat + 
                        re(random = ~1 | ring) + re(random = ~1 | year),
                  nu.formula = ~col_lat, 
                  tau.formula = ~col_lat, 
                  family = BEINF,
                  data = col_behav)

#save(colProp_BEINF, file = "Data_outputs/col_gamlss.RData")
load("Data_outputs/col_gamlss.RData")


# Fit raw minutes 
# These data are bounded but I can't see a way round this - no packages for Tobit
# models with random effects

## Set truncated parameters - this doesn't work due to error in the package
# trunc <- gamlss.tr::trun(par = c(0,1440),
#              type = "both",
#              family = TF)
# 
# ## Fit model
# colMins_TR <- gamlss(col_mins ~ col_lat + re(random = ~1 | ring),
#                      family = trunc,
#                      data = col_behav)
# 
# confint(colMins_TR)
# 
# #save(colMins_TR, file = "Data_outputs/col_truncated.RData")
# load("Data_outputs/col_truncated.RData")

colMins_lmm <- lmer(col_mins ~ col_lat + (1|ring),
                     data = col_behav)
colMins_lmm.null <- lmer(col_mins ~ 1 + (1|ring),
                     data = col_behav)

plot(colMins_lmm)
anova(colMins_lmm, colMins_lmm.null)


#### Explore results -----------------------------------------------------------

##### ~ ## Proportional ## -----------------------------------------------------

# Check diagnostic plots
plot(colProp_BEINF) 

# Model summary and confidence intervals (using broom.mixed) - quite slow
summary(colProp_BEINF)
colProp_BEINF.summary <- broom.mixed::tidy(colProp_BEINF, conf.int = TRUE)
colProp_BEINF.summary

# parameter term        estimate std.error statistic   p.value conf.low conf.high
# <chr>     <chr>          <dbl>     <dbl>     <dbl>     <dbl>    <dbl>     <dbl>
# 1 mu        (Intercept)  -0.671   0.0240      -27.9  2.93e-171  -0.718    -0.624 
# 2 mu        col_lat       0.0133  0.000342     39.0  0           0.0127    0.0140
# 3 sigma     (Intercept)  -0.0279  0.00246     -11.3  7.61e- 30  -0.0327   -0.0230
# 4 nu        (Intercept)   0.527   0.0740        7.12 1.08e- 12   0.381     0.672 
# 5 nu        col_lat      -0.0385  0.00108     -35.6  2.40e-276  -0.0407   -0.0364
# 6 tau       (Intercept)  -4.54    0.0809      -56.1  0          -4.70     -4.38  
# 7 tau       col_lat       0.0346  0.00113      30.7  6.95e-206   0.0323    0.0368


## Calculate probabilities for plotting

latitudes <- seq(45, 80, by = 5)

# Estimates and confint - mu #
intercept_mu <- colProp_BEINF.summary$estimate[1]
beta_collat_mu <- colProp_BEINF.summary$estimate[2]
std_error_intercept_mu <- colProp_BEINF.summary$std.error[1]
std_error_beta_collat_mu <- colProp_BEINF.summary$std.error[2]

# Create dateframe
mu_latitudes <- data.frame(latitude = latitudes)

# Calculate the critical value for the desired confidence level
z_critical <- qnorm(1 - 0.05 / 2)

# Calculate standard error 
mu_latitudes$std_error_mean <- sqrt((std_error_intercept_mu^2) + (mu_latitudes$latitude^2) * (std_error_beta_collat_mu^2))

# Calculate margin of error
mu_latitudes$margin_error <- z_critical * mu_latitudes$std_error_mean

# Calculate lower and upper bounds of confidence interval
mu_latitudes$lower_bound <- logit2prob(intercept_mu + beta_collat_mu * mu_latitudes$latitude - mu_latitudes$margin_error)
mu_latitudes$upper_bound <- logit2prob(intercept_mu + beta_collat_mu * mu_latitudes$latitude + mu_latitudes$margin_error)

# Convert the bounds from logit scale to probabilities
mu_latitudes$lower_bound_prob <- plogis(mu_latitudes$lower_bound)
mu_latitudes$upper_bound_prob <- plogis(mu_latitudes$upper_bound)

# Calculate probabilities
mu_latitudes$prob <- plogis(intercept_mu + beta_collat_mu * mu_latitudes$latitude)



# Estimates and confint - nu and tau #
intercept_nu <- colProp_BEINF.summary$estimate[4]
beta_collat_nu <- colProp_BEINF.summary$estimate[5]
std_error_intercept_nu <- colProp_BEINF.summary$std.error[4]
std_error_beta_collat_nu <- colProp_BEINF.summary$std.error[5]

intercept_tau <- colProp_BEINF.summary$estimate[6]
beta_collat_tau <- colProp_BEINF.summary$estimate[7]
std_error_intercept_tau <- colProp_BEINF.summary$std.error[6]
std_error_beta_collat_tau <- colProp_BEINF.summary$std.error[7]

# Create dataframe
parameters_latitudes <- data.frame(latitude = latitudes)

# Calculate the critical value for the desired confidence level
z_critical <- qnorm(1 - 0.05 / 2)

# Calculate standard errors
parameters_latitudes$std_error_nu <- sqrt((std_error_intercept_nu^2) + (parameters_latitudes$latitude^2) * (std_error_beta_collat_nu^2))
parameters_latitudes$std_error_tau <- sqrt((std_error_intercept_tau^2) + (parameters_latitudes$latitude^2) * (std_error_beta_collat_tau^2))

# Calculate margin of error for nu (shape) and tau (scale)
parameters_latitudes$margin_error_nu <- z_critical * parameters_latitudes$std_error_nu
parameters_latitudes$margin_error_tau <- z_critical * parameters_latitudes$std_error_tau

# Calculate the lower and upper bounds of confidence interval
parameters_latitudes$lower_bound_nu <- intercept_nu + beta_collat_nu * parameters_latitudes$latitude - parameters_latitudes$margin_error_nu
parameters_latitudes$upper_bound_nu <- intercept_nu + beta_collat_nu * parameters_latitudes$latitude + parameters_latitudes$margin_error_nu
parameters_latitudes$lower_bound_tau <- intercept_tau + beta_collat_tau * parameters_latitudes$latitude - parameters_latitudes$margin_error_tau
parameters_latitudes$upper_bound_tau <- intercept_tau + beta_collat_tau * parameters_latitudes$latitude + parameters_latitudes$margin_error_tau

# Exponentiate the bounds for nu and tau to obtain geometric mean and scale values
parameters_latitudes$lower_bound_nu_exp <- exp(parameters_latitudes$lower_bound_nu)
parameters_latitudes$upper_bound_nu_exp <- exp(parameters_latitudes$upper_bound_nu)
parameters_latitudes$lower_bound_tau_exp <- exp(parameters_latitudes$lower_bound_tau)
parameters_latitudes$upper_bound_tau_exp <- exp(parameters_latitudes$upper_bound_tau)

# Get exponentiates probabilities
parameters_latitudes$nu_prob <- exp(intercept_nu + beta_collat_nu * parameters_latitudes$latitude)
parameters_latitudes$tau_prob <- exp(intercept_tau + beta_collat_tau * parameters_latitudes$latitude)


##### ~ ## Minutes ## ----------------------------------------------------------

summary(colMins_lmm)
colMins_lmm.conf <- confint(colMins_lmm)

### Plotting ====================================================================

#### Plot proportional effects --------------------------------------------------

mu_plot.prob <- ggplot() +
  geom_ribbon(data = mu_latitudes, aes(x = latitude, ymin = lower_bound , ymax = upper_bound), alpha = 0.25, col = "grey",
              linetype = 2) +
  geom_line(data = mu_latitudes, aes(x = latitude, y = prob),
            col = "brown1", lwd = 0.4) +
  labs(y = "Proportion daylight hours at colony",
       x = expression("Colony latitude ("*~degree*" north)"),
       tag = "(A)") +
  geom_rug(data = col_behav, aes(x = col_lat, y = col_att.day), sides = "b") +
  ylim(c(0.45, 0.65)) +
  BLKI_theme



nu_plot.prob <- ggplot() +
  geom_ribbon(data = parameters_latitudes, aes(x = latitude, ymin = lower_bound_nu_exp , ymax = upper_bound_nu_exp), 
              alpha = 0.25, col = "grey",
              linetype = 2) +
  geom_line(data = parameters_latitudes, aes(x = latitude, y = nu_prob),
            col = "brown1", lwd = 0.4) +
  labs(y = "p | spend all daylight at the colony",
       x = expression("Colony latitude ("*~degree*" north)"),
       tag = "(B)") +
  geom_rug(data = col_behav, aes(x = col_lat, y = col_att.day), sides = "b") +
  ylim(c(0, 0.4)) +
  BLKI_theme



tau_plot.prob <- ggplot() +
  geom_ribbon(data = parameters_latitudes, aes(x = latitude, ymin = lower_bound_tau_exp , ymax = upper_bound_tau_exp), 
              alpha = 0.25, col = "grey",
              linetype = 2) +
  geom_line(data = parameters_latitudes, aes(x = latitude, y = tau_prob),
            col = "brown1", lwd = 0.4) +
  labs(y = "p | spend none of daylight at the colony",
       x = expression("Colony latitude ("*~degree*" north)"),
       tag = "(C)") +
  geom_rug(data = col_behav, aes(x = col_lat, y = col_att.day), sides = "b") +
  ylim(c(0, 0.25)) +
  BLKI_theme


#### ~ FIGURE 6 ~ Effects of colony latitude on proportional colony attendance ----

png("Figures/Figure 6_proportional_attendance.png", width = 12, height = 5, units = "in", res = 300)
ggarrange(mu_plot.prob, nu_plot.prob, tau_plot.prob,
          nrow = 1,
          align = "h")
dev.off()


#### Plot minutes effects ------------------------------------------------------

# Make predictions
summary(colMins_lmm)

colMins_lmm.df <- data.frame(ggpredict(colMins_lmm, terms = "col_lat")) %>% 
  rename(col_lat = x) 

# Build the plot
colMins_lmm.plot <- ggplot() +
  geom_point(data = col_behav %>% sample_frac(0.1), aes(x = col_lat, y = col_mins)) +
  geom_ribbon(data = colMins_lmm.df, aes(y = predicted, x = col_lat,
                                           ymin = conf.low, ymax = conf.high),
              alpha = 0.5, fill = "coral") +
  geom_line(data = colMins_lmm.df, aes(y = predicted, x = col_lat),
            linewidth = 1, col = "red") +
  labs(y = "Minutes spent at colony", x = "Colony latitude (°N)") +
  BLKI_theme


#### ~ FIGURE 7 ~ Effects of colony latitude on minutes colony attendance -------

png("Figures/Figure 7_minutes_attendance.png", width = 7, height = 5, units = "in", res = 300)
colMins_lmm.plot
dev.off()


# Model breeding success against latitude ======================================

### Load and process data -------------------------------------------------------
breeding <- read.csv("Data_inputs/BLKI_breeding-success.csv")
# Fix encoding issues
breeding$colony[breeding$colony == "r\xf8st"] <- "røst"
breeding$colony[breeding$colony == "bj\xf8rn\xf8ya"] <- "bjørnøya"
breeding$colony[breeding$colony == "skj\xe1lfandi"] <- "skjalfandi"

load("Data_inputs/BLKI_metadata.RData")
meta$colony <- tolower(meta$colony)
meta %<>% distinct(colony, .keep_all = TRUE) %>% dplyr::select(-ring)

breeding <- merge(breeding, meta, by = "colony")


## Get mean number of nests per colony and year
mean_nests <- breeding %>% group_by(colony) %>% summarise(mean_nests = mean(n_nests, na.rm = T))
mean(mean_nests$mean_nests, na.rm = TRUE)
range(mean_nests$mean_nests, na.rm = TRUE)

#### Visualise the breeding success relationship ---------------------------------

# Get aggregate and variance breeding success by colony
breeding_sum <- breeding %>% group_by(colony, col_lat, col_lon) %>% 
  summarise(mean_success = mean(large_chicks_per_nest, na.rm = TRUE),
            n_nests = sum(n_nests),
            var_success = var(large_chicks_per_nest, na.rm = TRUE)*(n_nests-1)/n_nests,
            sample_var = var(large_chicks_per_nest, na.rm = TRUE))

breeding_sum <- breeding_sum[order(breeding_sum$col_lat),]

# Build plots of mean and variance
## Mean plot
ggplot() +
  geom_point(data = breeding, aes(x = col_lat, y = large_chicks_per_nest)) +
  geom_smooth(method = "lm", data = breeding_sum, aes(x = col_lat, y = mean_success), 
              fill = "coral", col = "red") +
  labs(y = "Number of chicks produced per nest", x = "Colony latitude (°N)") +
  scale_y_continuous(breaks = seq(0, 1.5, 0.5)) +
  ylim(0, 1.5) +
  BLKI_theme
  

## Variance plot
breeding_var.plot <- ggplot() +
 geom_point(data = breeding_sum, aes(x = col_lat, y = var_success)) +
 geom_smooth(method = "lm", data = breeding_sum, aes(x = col_lat, y = sample_var), 
                fill = "coral", col = "red") +
 labs(y = "Variance in number of chicks produced per nest", x = "Colony latitude (°N)") +
 BLKI_theme


#### ~ FIGURE S5: Breeding success variation with latitude --------------------------------

png("Figures/Figure S5_breeding_variation.png", width = 7, height = 5, units = "in", res = 300)
breeding_var.plot
dev.off()

### Construct model -------------------------------------------------------------

# Remove NA values
breeding_mod <- subset(breeding, !is.na(large_chicks_per_nest))

# Specify the model
breeding_lmm <- lmerTest::lmer(large_chicks_per_nest ~ col_lat + (1|colony), 
                          data = breeding_mod)
breeding_lmm.null <- lmer(large_chicks_per_nest ~ 1 + (1|colony), 
                      data = breeding_mod)

anova(breeding_lmm, breeding_lmm.null)
summary(breeding_lmm)
confint(breeding_lmm)

# Make predictions and plot them
breeding_glmm.df <- data.frame(ggpredict(breeding_lmm, terms = "col_lat")) %>% 
  rename(col_lat = x) 


#### Build the final plot ------------------------------------------------------
breeding_glmm.plot <- ggplot() +
  geom_jitter(data = breeding, aes(x = col_lat, y = large_chicks_per_nest),
             width = 0.3) +
  # geom_ribbon(data = breeding_glmm.df, aes(y = predicted, x = col_lat,
  #                                          ymin = conf.low, ymax = conf.high),
  #             alpha = 0.5, fill = "coral") +
  # geom_line(data = breeding_glmm.df, aes(y = predicted, x = col_lat),
  #           linewidth = 1, col = "red", linetype = "dashed") +
  labs(y = "Number of chicks produced per nest", x = "Colony latitude (°N)") +
  BLKI_theme

#### ~ FIGURE 8: Breeding success with latitude --------------------------------

png("Figures/Figure 8_breeding_success.png", width = 7, height = 5, units = "in", res = 300)
breeding_glmm.plot
dev.off()


# APPENDIX ----------------------------------------------------------------


## Construct multivariate model to examine rest/flight/forage simultaneously

# Melt responses into a long format
behaviour_long <- melt(multi_behav, id.vars = c("col_lat", "ring", "n_recs24", "year", "daylight.mins"),
                       measure.vars = c("prop_forage.day", "prop_rest.day", "prop_flight.day"),
                       variable.name = "behaviour", value.name = "proportion")

# Construct the model
multiMod <- glmmTMB(proportion ~ daylight.mins * behaviour + (1|ring), 
                    family = beta_family(link = "logit"),
                    ziformula = ~.,
                    data = behaviour_long)



## Get interaction effects

## Intercept and overall values for each behaviour don't matter - these are 
## estimates of proportion when col_lat = 0. Only interested in effect of col lat 
## and interaction

# Critical value for confidence intervals
critical_value <- qnorm((1 + 0.95) / 2)

# Conditional model
cond_coefs <- multiMod_summary$coefficients$cond
cond_vcov <- vcov(multiMod)$cond


## Rest
calc_test_statistics(cond_coefs, cond_vcov, 2, 5)
#   beta           se   lower_conf   upper_conf   z_score p_value
#1 -0.009127003 0.0007502676 -0.0105975 -0.007656505 -12.165       0
## Flight
calc_test_statistics(cond_coefs, cond_vcov, 2, 6)
#    beta           se   lower_conf    upper_conf   z_score    p_value
# 1 -0.003233648 0.0007100692 -0.004625358 -0.001841938 -4.55399 5.263795e-06


# Estimated marginal means for foraging and flight
latitudes <- seq(45, 80, by = 5)
intercept <-  -0.183 
beta_collat <- -0.0134072

means_cond_forage <- data.frame(latitude = latitudes)
means_cond_forage$prob <- logit2prob(intercept + beta_collat * means_cond_forage$latitude)
(means_cond_forage$prob[1]*100) - (means_cond_forage$prob[2]*100)

# Same for flight
beta_collat <- -0.001315552

means_cond_flight <- data.frame(latitude = latitudes)
means_cond_flight$prob <- logit2prob(intercept + beta_collat * means_cond_flight$latitude)
(means_cond_flight$prob[1]*100) - (means_cond_flight$prob[2]*100)


effects_forage <- ggeffects::ggpredict(multiMod, terms = c("col_lat", "behaviour [prop_forage.day]"))
effects_flight <- ggeffects::ggpredict(multiMod, terms = c("col_lat", "behaviour [prop_flight.day]"))


# Zero inflated model
zi_coefs <- multiMod_summary$coefficients$zi
ze_vcov <- vcov(multiMod)$zi

# Rest
calc_test_statistics(zi_coefs, zi_vcov, 2, 5)
#    beta          se  lower_conf  upper_conf   z_score    p_value
# 1 -0.004679466 0.003267947 -0.01108452 0.001725593 -1.431928 0.1521643
# Flight
calc_test_statistics(zi_coefs, zi_vcov, 2, 6)
#   beta           se   lower_conf    upper_conf   z_score    p_value
# 1 -0.01501392 0.003258953 -0.02140135 -0.008626493 -4.606978 4.085636e-06


# Estimated marginal means for zero-inflation foraging
latitudes <- seq(45, 80, by = 5)
intercept <-  -2.18 
beta_collat <- -0.0346

means_zi_forage <- data.frame(latitude = latitudes)
means_zi_forage$prob <- logit2prob(intercept + beta_collat * means_zi_forage$latitude)

# Same for flight
latitudes <- seq(45, 80, by = 5)
intercept <-  -2.18 
beta_collat <- -0.02427242

means_zi_flight <- data.frame(latitude = latitudes)
means_zi_flight$prob <- logit2prob(intercept + beta_collat * means_zi_flight$latitude)


## Tidied coefficients table

broom.mixed::tidy(multiMod, conf.int = TRUE)

# effect   component group term                             estimate std.error statistic    p.value conf.low conf.high
# <chr>    <chr>     <chr> <chr>                               <dbl>     <dbl>     <dbl>      <dbl>    <dbl>     <dbl>
#  1 fixed    cond      NA    (Intercep… -0.102    0.0593       -1.73  8.40e-  2 -0.219     0.0138 
#  2 fixed    cond      NA    col_lat    -0.0163   0.000709    -22.9   1.72e-116 -0.0177   -0.0149 
#  3 fixed    cond      NA    behaviour… -1.19     0.0321      -37.1   2.13e-301 -1.25     -1.13   
#  4 fixed    cond      NA    behaviour… -1.07     0.0271      -39.6   0         -1.13     -1.02   
#  5 fixed    cond      NA    col_lat:b…  0.00714  0.000460     15.5   2.70e- 54  0.00624   0.00804
#  6 fixed    cond      NA    col_lat:b…  0.0130   0.000389     33.5   4.17e-246  0.0123    0.0138 
#  7 fixed    zi        NA    (Intercep… -4.29     0.244       -17.6   1.69e- 69 -4.77     -3.82   
#  8 fixed    zi        NA    col_lat     0.00874  0.00326       2.68  7.31e-  3  0.00235   0.0151 
#  9 fixed    zi        NA    behaviour…  3.89     0.147        26.5   1.04e-154  3.61      4.18   
# 10 fixed    zi        NA    behaviour… -0.656    0.424        -1.55  1.22e-  1 -1.49      0.175  
# 11 fixed    zi        NA    col_lat:b… -0.0134   0.00209      -6.43  1.28e- 10 -0.0175   -0.00933
# 12 fixed    zi        NA    col_lat:b… -0.0238   0.00612      -3.88  1.04e-  4 -0.0358   -0.0118 
# 13 ran_pars cond      ring  sd__(Inte…  0.149   NA            NA    NA          0.142     0.157  
# 14 ran_pars cond      year  sd__(Inte…  0.0972  NA            NA    NA          0.0623    0.152  
# 15 ran_pars zi        ring  sd__(Inte…  0.594   NA            NA    NA          0.564     0.626  
# 16 ran_pars zi        year  sd__(Inte…  0.235   NA            NA    NA          0.145     0.379 

