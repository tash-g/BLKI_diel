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
  filter(month >= 5 & month <= 7) %>% dplyr::select(-month)

# Get year
gls_daily %<>% mutate(year = format(as.Date(date, format = "%Y-%m-%d"),"%Y"))



### Colony attendance data -----------------------------------------------------

# Prepare data
gls_daily$col_hrs <- gls_daily$col_att24*gls_daily$daylight.mins

col_behav <- gls_daily %>%
  dplyr::select(ring, date, col_att.day, n_recs.day, col_lat, year, col_hrs, colony) %>% 
  filter(!is.na(n_recs.day)) %>%
  mutate(col_hrs = col_hrs/60) %>% select(-col_hrs)

# ==============================================================================
# Latitudinal variation in rhythmicity and period length =======================

## Main analyses in 5_individual-LombScargle_periodograms

# Visualise individual-level rhythms -------------------------------------------

# Find length of time series for each nest
bird.ts <- tapply(kits.att_all$datetime, list(kits.att_all$ring), function(x)(length(x))/2)
mean(bird.ts); sd(bird.ts); qt(0.975,df=length(bird.ts)-1)*sd(bird.ts)/sqrt(length(bird.ts)) 

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
colonies <- kits.att_all %>% dplyr::select(c(ring, colony, col_lat)) %>% distinct(ring, .keep_all = TRUE) 
ACF.df <- merge(ACF.df, colonies, by = "ring")


### Plot ACF for each individual by colony -------------------------------------

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
    ggtitle(paste0(colonies[i], ", ", round(acf_plot$col_lat[1], digits = 2))) +
    theme(legend.position = "none",
          plot.margin = unit(c(0, 0.25, 0, 0), "cm"),
          plot.title = element_text(size = 20, hjust = 1),
          axis.text.x = element_text(size = 18), 
          axis.text.y = element_text(size = 18), 
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          plot.title.position = 'plot')
  
  if (i %in% c(1:6)) {
    ACF.ind_plot <- ACF.ind_plot + theme(plot.margin = unit(c(0.5, 0.25, 0, 0), "cm"))
  }
  
  if (!(i %in% c(1, 7, 13, 19))) {
    ACF.ind_plot <- ACF.ind_plot + theme(axis.text.y = element_blank(), axis.title.y = element_blank())
  }
  
  if (!(i %in% c(17:22))) {
    ACF.ind_plot <- ACF.ind_plot + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
  }
  
  plot_list[[i]] <- ACF.ind_plot  
  
}

### ~ FIGURE S6 ~ ACF in immersion by individuals over colonies ----------------
png("Figures/FigureS6_individual_imm_ACF.png", 
    width = 25, height = 20, units = "in", res = 600)
cowplot::plot_grid(plotlist = plot_list, nrow = 4, align = "hv",
                   ncol = 6)
dev.off()



### ~ FIGURE S8 ~ Plot colony prop immersion by latitude -----------------------
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
    ggtitle(paste0(colonies[i], ", ", round(pop.imm_plot$col_lat[1], digits = 2))) +
    labs(x = "Hour of day", y = "Mean proportional immersion") +
    scale_x_continuous(breaks = seq(0, 24, 4), limits = c(0, 23)) +
    BLKI_theme +
    theme(legend.position = "none",
          plot.margin = unit(c(0, 0.25, 0, 0), "cm"),
          plot.title = element_text(size = 20, hjust = 1),
          axis.text.x = element_text(size = 18), 
          axis.text.y = element_text(size = 18), 
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          plot.title.position = 'plot')
  
  if (i %in% c(1:6)) {
    pop.imm_figure <- pop.imm_figure + theme(plot.margin = unit(c(0.5, 0.25, 0, 0), "cm"))
  }
  
  if (!(i %in% c(1, 7, 13, 19))) {
    pop.imm_figure <- pop.imm_figure + theme(axis.text.y = element_blank(), axis.title.y = element_blank())
  }
  
  if (!(i %in% c(17:22))) {
    pop.imm_figure <- pop.imm_figure + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
  }
  
  plot_list[[i]] <- pop.imm_figure  
  
}

png("Figures/FigureS8_imm_by_colony.png", 
    width = 25, height = 20, units = "in", res = 650)
cowplot::plot_grid(plotlist = plot_list, nrow = 4, align = "hv",
                   ncol = 6)
dev.off()


## Run autocorrelatIon function on each colony ----------------------------------
colonies <- unique(kits.pop.imm$colony[order(kits.pop.imm$col_lat)])
plot_list <- vector(mode = "list", length = length(colonies))
acf_data_list <- vector(mode = "list", length = length(colonies))

for (i in 1:length(colonies)) {
  
  print(paste0("Processing colony ", i, " of ", length(colonies)))
  
  pop.imm <- subset(kits.pop.imm, colony == colonies[i])
  
  imm.ACF <- acf(pop.imm$prop.imm, na.action = na.pass, lag = 60, plot = F)
 
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
    ggtitle(paste0(colonies[i], ", ", round(pop.imm$col_lat[1], digits = 2))) +
    theme(plot.margin = unit(c(0, 0.25, 0, 0), "cm"),
          plot.title = element_text(size = 20, hjust = 1),
          axis.text.x = element_text(size = 18), 
          axis.text.y = element_text(size = 18), 
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20))
  
  if (i %in% c(1:6)) {
    imm.ACF_plot <- imm.ACF_plot + theme(plot.margin = unit(c(0.5, 0.25, 0, 0), "cm"))
  }

  if (!(i %in% c(1, 7, 13, 19))) {
    imm.ACF_plot <- imm.ACF_plot + theme(axis.text.y = element_blank(), axis.title.y = element_blank())
  }

  if (!(i %in% c(17:22))) {
    imm.ACF_plot <- imm.ACF_plot + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
  }

  plot_list[[i]] <- imm.ACF_plot  
  
}

acf_col <- do.call("rbind", acf_data_list)
#save(acf_col, file = "Data_outputs/acf_colony_data.RData")
load("Data_outputs/acf_colony_data.RData")


### ~ FIGURE S9 ~ ACF in immersion over all colonies --------------------------
png("Figures/FigureS9_colony_imm_ACF.png", 
    width = 25, height = 20, units = "in", res = 600)
cowplot::plot_grid(plotlist = plot_list, nrow = 4, align = "hv",
                   ncol = 6)
dev.off()




# ACF peak strength and lag with colony latitude -------------------------------

acf_col$peak_strength <- abs(acf_col$peak_strength)

# Peak lag
peak_lag_lm <- lm(peak_lag ~ col_lat, data = acf_col)
summary(peak_lag_lm)
median(acf_col$peak_lag)
sd(acf_col$peak_lag)

# Peak strength
peak_strength_lm <- lm(log(peak_strength) ~ col_lat, data = acf_col)
summary(peak_strength_lm)
exp(summary(peak_strength_lm)$coefficients)
exp(coef(peak_strength_lm))
mean(acf_col$peak_strength)



### ~ FIGURE 4 ~ Peak ACF strength by latitude --------------------------------

peak_strength_lm.df <- data.frame(ggpredict(peak_strength_lm, terms = "col_lat", back_transform = F)) %>% 
  dplyr::rename(col_lat = x) 

# Build the plot
peak_strength_lm.plot <- ggplot() +
  geom_point(data = acf_col, aes(x = col_lat, y = log(peak_strength))) +
  geom_ribbon(data = peak_strength_lm.df, aes(y = predicted, x = col_lat,
                                              ymin = conf.low, ymax = conf.high),
              alpha = 0.5, fill = "grey") +
  geom_line(data = peak_strength_lm.df, aes(y = predicted, x = col_lat),
            linewidth = 1, col = "darkorange") +
  labs(x = expression("Colony latitude °N"), y = "log(ACF peak strength)") +
  scale_x_continuous(breaks = seq(45, 80, by = 5)) + 
  BLKI_theme 

png("Figures/Figure4_LM_peak_by_lat.png", 
    width = 9, height = 7, units = "in", res = 600)
peak_strength_lm.plot
dev.off()



# GAM: immersion values ~ sun elevation ----------------------------------------

## Subset the data to random 72 hour samples to make manageable
set.seed(817)
random_sample <- subsample_by_hours(kits.att_all, 72)

random_sample[,c(3, 7)] <- lapply(random_sample[,c(3, 7)], as.factor)
random_sample$ring <- as.factor(random_sample$ring)
random_sample$date <- as.Date(random_sample$datetime)

#save(random_sample, file = "Data_outputs/subsampled_kits.RData")
load("Data_outputs/subsampled_kits.RData")


## Do the same for Arctic colonies 
high_lat <- subset(kits.att_all, col_lat > 66.5)

# Subset the data
set.seed(718)
high_lat.random <- subsample_by_hours(high_lat, 72)

high_lat.random[,c(3, 7)] <- lapply(high_lat.random[,c(3, 7)], as.factor)
high_lat.random$ring <- as.factor(high_lat.random$ring)
high_lat.random$date <- as.Date(high_lat.random$datetime)

#save(high_lat.random, file = "Data_outputs/high_lat_subsampled_kits.RData")
load("Data_outputs/high_lat_subsampled_kits.RData")



### Fit the model --------------------------------------------------------------

gam_model <- mgcv::gam(total_imm ~ s(sun_altitude) +
                   s(ring, bs = 're'), 
                 data = random_sample)

#save(gam_model, file = "Data_outputs/GAM_immersion_sun.RData")
load("Data_outputs/GAM_immersion_sun.RData")

gam.check(gam_model)
summary(gam_model)

## Do the same for Arctic colonies 
high_lat.gam_model <- mgcv::gam(total_imm ~ s(sun_altitude) +
                         s(ring, bs = 're'), 
                       data = high_lat.random)

#save(high_lat.gam_model, file = "Data_outputs/high_lat_GAM_immersion_sun.RData")
load("Data_outputs/high_lat_GAM_immersion_sun.RData")

gam.check(high_lat.gam_model)
summary(high_lat.gam_model)
plot(high_lat.gam_model)

### ~ FIGURE S10 ~ Immersion against sun elevation angle -------------------------

# Process data for plotting
p_obj <- plot(gam_model, residuals = TRUE)
p_obj <- p_obj[[1]] # just one smooth so select the first component
sm_df <- as.data.frame(p_obj[c("x", "se", "fit")])
sm_df$fit <- sm_df$fit + coef(gam_model)[1]
data_df <- as.data.frame(p_obj[c("raw", "p.resid")])

# Build the plot
gam_immersion.plot <- ggplot(sm_df, aes(x = x, y = fit)) +
  geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
              alpha = 0.5, fill = "grey") +
  geom_line(col = "darkorange", lwd = 1) +
  labs(x = expression("Sun elevation angle ("*~degree*")"),
       y = "Predicted hourly immersion",
       tag = "(A)") +
  BLKI_theme 

### High latitude sun elevation 

# Process data for plotting
p_obj.high <- plot(high_lat.gam_model, residuals = TRUE)
p_obj.high <- p_obj.high[[1]] # just one smooth so select the first component
sm_df.high <- as.data.frame(p_obj.high[c("x", "se", "fit")])
sm_df.high$fit <- sm_df.high$fit + coef(high_lat.gam_model)[1]
data_df.high <- as.data.frame(p_obj.high[c("raw", "p.resid")])

# Build the plot
high_lat_gam_immersion.plot <- ggplot(sm_df.high, aes(x = x, y = fit)) +
  geom_ribbon(aes(ymin = fit - se, ymax = fit + se, y = NULL),
              alpha = 0.5, fill = "grey") +
  geom_line(col = "darkorange", lwd = 1) +
  labs(x = expression("Sun elevation angle ("*~degree*")"),
       y = "",
       tag = "(B)") +
  BLKI_theme


png("Figures/FigureS11_GAM_imm_by_sun.png", 
    width = 12, height = 7, units = "in", res = 600)
ggarrange(gam_immersion.plot, high_lat_gam_immersion.plot,
          ncol = 2,
          align = "h")
dev.off()



# ==============================================================================
# Latitudinal variation in behaviour ===========================================

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


### ~ FIGURE S4 ~ Foraging autocorrelaiton -------------------------------------
png("Figures/FigureS4_acf_forage_plot.png", width = 9, height = 7, units = "in", res = 300)
acf_forage_plot
dev.off()


## Latitudinal variation in behaviour - at-sea ---------------------------------

### Construct models -----------------------------------------------------------

# Remove full day colony visits
multi_behav <- gls_daily %>% filter(col_att.day < 0.95 & prop_flight.day < 0.95 & prop_rest.day < 0.95 & prop_forage.day < 0.95)

## Calculate hours
multi_behav$forage_hrs <- round(multi_behav$prop_forage24 * 24)
multi_behav$rest_hrs <- round(multi_behav$prop_rest24 * 24)
multi_behav$flight_hrs <- round(multi_behav$prop_flight24 * 24)

# Get daily means
mean(multi_behav$prop_forage.day); sd(multi_behav$prop_forage.day)
mean(multi_behav$prop_flight.day); sd(multi_behav$prop_flight.day)
mean(multi_behav$prop_rest.day); sd(multi_behav$prop_rest.day)

## Construct models for rest/flight/forage separately - proportion
load("Data_outputs/forageProp_glmm.RData")
load("Data_outputs/restProp_glmm.RData")
load("Data_outputs/flightProp_glmm.RData")

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


## Construct models for rest/flight/forage separately

## Remove NA values
multi_behav %<>% filter(!is.na(forage_hrs) & !is.na(rest_hrs) & !is.na(flight_hrs))

# Fit models
load("Data_outputs/forageHrs_glmm.RData")
load("Data_outputs/restHrs_glmm.RData")
load("Data_outputs/flightHrs_glmm.RData")

forageHrs_glmm <- glmmTMB(forage_hrs ~ col_lat + (1|ring) + (1|year),
                           data = multi_behav, family = nbinom1)
forageHrs_glmm.null <- glmmTMB(forage_hrs ~ 1 + (1|ring) + (1|year),
                           data = multi_behav, family = nbinom1)

anova(forageHrs_glmm, forageHrs_glmm.null)
summary(forageHrs_glmm)

performance::check_overdispersion(forageHrs_glmm)



restHrs_glmm <- glmmTMB(rest_hrs ~ col_lat + (1|ring) + (1|year),
                           data = multi_behav, family = nbinom1)
restHrs_glmm.null <- glmmTMB(rest_hrs ~ 1 + (1|ring) + (1|year),
                                data = multi_behav, family = nbinom1)

anova(restHrs_glmm, restHrs_glmm.null)
summary(restHrs_glmm)

performance::check_overdispersion(restHrs_glmm)



flightHrs_glmm <- glmmTMB(flight_hrs ~ col_lat + (1|ring) + (1|year),
                         data = multi_behav, family = nbinom1)
flightHrs_glmm.null <- glmmTMB(flight_hrs ~ 1 + (1|ring) + (1|year),
                              data = multi_behav, family = nbinom1)

anova(flightHrs_glmm, flightHrs_glmm.null)
summary(flightHrs_glmm)

performance::check_overdispersion(flightHrs_glmm)


# save(forageHrs_glmm, file = "Data_outputs/forageHrs_glmm.RData")
# save(restHrs_glmm, file = "Data_outputs/restHrs_glmm.RData")
# save(flightHrs_glmm, file = "Data_outputs/flightHrs_glmm.RData")


#### Explore results -----------------------------------------------------------

latitudes <- seq(45, 80, by = 5)

## Proportional results
# Change coefficients to $cond = conditional; $zi = zero inflated coefficients

# Foraging #
tab_model(forageProp_glmm)
cond_coefs.forProp <-  summary(forageProp_glmm)$coefficients$cond
confint(forageProp_glmm)

intercept <-  cond_coefs.forProp[1,1]
beta_collat <- cond_coefs.forProp[2,1]

means_cond_forage <- data.frame(latitude = latitudes)
means_cond_forage$prob <- logit2prob(intercept + beta_collat * means_cond_forage$latitude)
((means_cond_forage$prob[1] - means_cond_forage$prob[2]) / means_cond_forage$prob[1] ) * 100

means_cond_forage$prob[1]*100; means_cond_forage$prob[nrow(means_cond_forage)]*100

# Rest #
tab_model(restProp_glmm)
cond_coefs.restProp <-  summary(restProp_glmm)$coefficients$cond
confint(restProp_glmm)

intercept <-  cond_coefs.restProp[1,1]
beta_collat <- cond_coefs.restProp[2,1]

means_cond_rest <- data.frame(latitude = latitudes)
means_cond_rest$prob <- logit2prob(intercept + beta_collat * means_cond_rest$latitude)
(( means_cond_rest$prob[2] - means_cond_rest$prob[1]) / means_cond_rest$prob[1] ) * 100 

means_cond_rest$prob[1]*100; means_cond_rest$prob[nrow(means_cond_rest)]*100

# Flight #
tab_model(flightProp_glmm)
cond_coefs.flightProp <-  summary(flightProp_glmm)$coefficients$cond
confint(flightProp_glmm)

intercept <-  cond_coefs.flightProp[1,1]
beta_collat <- cond_coefs.flightProp[2,1]

means_cond_flight <- data.frame(latitude = latitudes)
means_cond_flight$prob <- logit2prob(intercept + beta_collat * means_cond_flight$latitude)
(( means_cond_flight$prob[1] - means_cond_flight$prob[2]) / means_cond_flight$prob[1] ) * 100 

means_cond_flight$prob[1]*100; means_cond_flight$prob[nrow(means_cond_flight)]*100

## Minute effects

# Forage - decrease #
summary(forageHrs_glmm)
tab_model(forageHrs_glmm)
exp(confint(forageHrs_glmm))
anova(forageHrs_glmm, forageHrs_glmm.null)

# Rest - decrease #
summary(restHrs_glmm)
tab_model(restHrs_glmm)
exp(confint(restHrs_glmm))
anova(restHrs_glmm, restHrs_glmm.null)

# Flight - increase #
summary(flightHrs_glmm)
tab_model(flightHrs_glmm)
exp(confint(flightHrs_glmm))
anova(flightHrs_glmm, flightHrs_glmm.null)


## Get effect sizes for each behaviour
effects_for.hrs <- data.frame(ggeffects::ggpredict(forageHrs_glmm, terms = "col_lat")) %>% 
  rename(col_lat = x, behaviour = group) %>%
  mutate(behaviour = "forage")

effects_rest.hrs <- data.frame(ggeffects::ggpredict(restHrs_glmm, terms = "col_lat")) %>% 
  rename(col_lat = x, behaviour = group) %>%
  mutate(behaviour = "rest")

effects_flight.hrs <- data.frame(ggeffects::ggpredict(flightHrs_glmm, terms = "col_lat")) %>% 
  rename(col_lat = x, behaviour = group) %>%
  mutate(behaviour = "flight")

### Plotting ===================================================================

# Get raw data for plotting
plot_dat <- multi_behav %>% dplyr::select(c(col_lat, ring,
                                            prop_forage.day, prop_rest.day, prop_flight.day,
                                            forage_hrs, rest_hrs, flight_hrs))

## Set to long format
### Proportion
plot_dat.prop <- melt(plot_dat, id.vars = c("col_lat", "ring"),
     measure.vars = c("prop_forage.day", "prop_rest.day", "prop_flight.day"),
     variable.name = "behaviour", value.name = "proportion")

level_map <- c("prop_forage.day" = "forage",
               "prop_rest.day" = "rest",
               "prop_flight.day" = "flight")

plot_dat.prop$behaviour <- level_map[plot_dat.prop$behaviour]
plot_dat.prop$behaviour <- factor(plot_dat.prop$behaviour, levels = c("forage", "rest", "flight"))

### Minutes
plot_dat.hrs <- melt(plot_dat, id.vars = c("col_lat", "ring"),
                      measure.vars = c("forage_hrs", "rest_hrs", "flight_hrs"),
                      variable.name = "behaviour", value.name = "hours")

level_map <- c("forage_hrs" = "forage",
               "rest_hrs" = "rest",
               "flight_hrs" = "flight")

plot_dat.hrs$behaviour <- level_map[plot_dat.hrs$behaviour]
plot_dat.hrs$behaviour <- factor(plot_dat.hrs$behaviour, levels = c("forage", "rest", "flight"))

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
       y = "Proportion daylight spent in behaviour",
       tag = "(A)") +
  BLKI_theme +
  theme(legend.position = c(0.07, 0.9),
        legend.background = element_rect(fill = "transparent"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())


#### Plot hour effects ---------------------------------------------------------

## Get effect sizes for each behaviour
effects_for.hrs <- data.frame(ggeffects::ggpredict(forageHrs_glmm, terms = "col_lat")) %>% 
  rename(col_lat = x, behaviour = group) %>%
  mutate(behaviour = "forage")

effects_rest.hrs <- data.frame(ggeffects::ggpredict(restHrs_glmm, terms = "col_lat")) %>% 
  rename(col_lat = x, behaviour = group) %>%
  mutate(behaviour = "rest")

effects_flight.hrs <- data.frame(ggeffects::ggpredict(flightHrs_glmm, terms = "col_lat")) %>% 
  rename(col_lat = x, behaviour = group) %>%
  mutate(behaviour = "flight")

# Bind together
hrs_effects <- rbind(effects_for.hrs,
                      effects_rest.hrs,
                      effects_flight.hrs)


#save(hrs_effects, file = "Data_outputs/rawHrs_plotting_data.RData")
load("Data_outputs/rawHrs_plotting_data.RData")


## Build the plot ##

hrs_behaviour_plot <- ggplot() + 
  geom_point(aes(x = col_lat, y = hours, group = behaviour, col = behaviour),
             position = position_jitter(w = 0.5),
             data = plot_dat.hrs[sample(nrow(plot_dat.hrs), round(nrow(plot_dat.hrs)/100) ),]) +
  geom_line(aes(x = col_lat, y = predicted, group = behaviour, col = behaviour), 
            data = hrs_effects) +
  geom_ribbon(aes(x = col_lat, ymin = conf.low, ymax = conf.high, group = behaviour),
              data = hrs_effects,
              alpha = 0.5, fill = "grey") +
  scale_colour_viridis_d(labels = c("Forage", "Rest", "Flight", "Colony"),
                      name = "Behaviour",
                      option = "plasma") +
  labs(x = expression("Colony latitude ("*~degree*" north)"), 
       y = "Hours per 24 hour-period spent in behaviour",
       tag = "(B)") +
  BLKI_theme +
  theme(legend.position = "none")



### ~ FIGURE 5 ~ At-sea raw hours behaviours as a function of latitude =========

png("Figures/Figure5_behaviour_by_latitude.png", width = 9, height = 14, units = "in", res = 600)
ggarrange(proportional_behaviour_plot, hrs_behaviour_plot, nrow = 2, align = "hv")
dev.off()


# Latitudinal variation in behaviour - colony attendance -----------------------

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


# Fit hourly response

## Scale colony latitude
col_behav$col_lat.sc <- scale(col_behav$col_lat)

colHrs_glmm <- glmmTMB(col_hrs ~ col_lat.sc + (1|ring),
                        data = col_behav,
                        family = nbinom1)
colHrs_glmm.null <- glmmTMB(col_hrs ~ 1 + (1|ring),
                             data = col_behav,
                             family = nbinom1)


#save(colHrs_glmm, file = "Data_outputs/col_glmm.RData")
load("Data_outputs/col_glmm.RData")

summary(colHrs_glmm)
overdisp_fun(colHrs_glmm)
performance::check_overdispersion(colHrs_glmm)




#### Explore results -----------------------------------------------------------

##### ~ ## Proportional ## -----------------------------------------------------

# Check diagnostic plots
plot(colProp_BEINF) 

# Model summary and confidence intervals (using broom.mixed) - quite slow
summary(colProp_BEINF)
tab_model(colProp_BEINF)
colProp_BEINF.summary <- broom.mixed::tidy(colProp_BEINF, conf.int = TRUE)

var(colProp_BEINF$mu.coefSmo[[1]]$coefficients$random$ring)

#save(colProp_BEINF.summary, file = "Data_outputs/col_gamlss.summary.RData")
load("Data_outputs/col_gamlss.summary.RData")

colProp_BEINF.summary

# # A tibble: 7 × 8
# parameter term        estimate std.error statistic   p.value conf.low conf.high
# <chr>     <chr>          <dbl>     <dbl>     <dbl>     <dbl>    <dbl>     <dbl>
#   1 mu        (Intercept)  0.587    0.0239       24.5  1.20e-132  0.540     0.634  
# 2 mu        col_lat     -0.00300  0.000342     -8.79 1.54e- 18 -0.00368  -0.00233
# 3 sigma     (Intercept) -0.0466   0.00244     -19.1  3.33e- 81 -0.0514   -0.0417 
# 4 nu        (Intercept) -4.88     0.132       -36.9  1.08e-296 -5.14     -4.62   
# 5 nu        col_lat      0.0236   0.00186      12.7  5.41e- 37  0.0200    0.0272 
# 6 tau       (Intercept) -3.06     0.0701      -43.7  0         -3.20     -2.92   
# 7 tau       col_lat      0.0169   0.000990     17.1  2.95e- 65  0.0149    0.0188

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

### Change in colony attendance for 5 degree increase in latitude
(mu_latitudes$prob[2] - mu_latitudes$prob[1])*100


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

# Get exponentiated probabilities
parameters_latitudes$nu_prob <- exp(intercept_nu + beta_collat_nu * parameters_latitudes$latitude)
parameters_latitudes$tau_prob <- exp(intercept_tau + beta_collat_tau * parameters_latitudes$latitude)


### Change in p no or constant colony attendance for 5 degree increase in latitude
(parameters_latitudes$nu_prob[2] - parameters_latitudes$nu_prob[1])*100
(parameters_latitudes$tau_prob[2] - parameters_latitudes$tau_prob[1])*100


##### ~ ## Minutes ## ----------------------------------------------------------

colHrs_glmm.summary <- summary(colHrs_glmm)
colHrs_glmm.conf <- confint(colHrs_glmm)
anova(colHrs_glmm, colHrs_glmmnull)

## Percentage increase for 5 degree change in latitude
intercept_glmm <- colHrs_glmm.summary$coefficients[1]
beta_collat_glmm <- colHrs_glmm.summary$coefficients[2]

(((intercept_glmm + beta_collat_glmm*50) - (intercept_glmm + beta_collat_glmm*45)) / 
  (intercept_glmm + beta_collat_glmm*50)) * 100




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
  labs(y = "p | spend none of daylight at the colony",
       x = expression("Colony latitude ("*~degree*" north)"),
       tag = "(B)") +
  geom_rug(data = col_behav, aes(x = col_lat, y = col_att.day), sides = "b") +
  ylim(c(0, 0.1)) +
  BLKI_theme



tau_plot.prob <- ggplot() +
  geom_ribbon(data = parameters_latitudes, aes(x = latitude, ymin = lower_bound_tau_exp , ymax = upper_bound_tau_exp), 
              alpha = 0.25, col = "grey",
              linetype = 2) +
  geom_line(data = parameters_latitudes, aes(x = latitude, y = tau_prob),
            col = "brown1", lwd = 0.4) +
  labs(y = "p | spend all of daylight at the colony",
       x = expression("Colony latitude ("*~degree*" north)"),
       tag = "(C)") +
  geom_rug(data = col_behav, aes(x = col_lat, y = col_att.day), sides = "b") +
  ylim(c(0, 0.3)) +
  BLKI_theme


#### ~ FIGURE 6 ~ Effects of colony latitude on proportional colony attendance ----

png("Figures/Figure6_proportional_attendance.png", width = 12, height = 5, units = "in", res = 300)
ggarrange(mu_plot.prob, nu_plot.prob, tau_plot.prob,
          nrow = 1,
          align = "h")
dev.off()


#### Plot hour effects ---------------------------------------------------------

# Make predictions
summary(colHrs_glmm)

colHrs_glmm.unsc <- glmmTMB(col_hrs ~ col_lat + (1|ring),
                        data = col_behav,
                        family = nbinom1)

#save(colHrs_glmm.unsc, file = "Data_outputs/col_glmm_unsc.RData")
load("Data_outputs/col_glmm_unsc.RData")

colHrs_glmm.df <- data.frame(ggeffects::ggpredict(colHrs_glmm.unsc, terms = "col_lat")) %>% 
  dplyr::rename(col_lat = x)

# Build the plot
colHrs_glmm.plot <- ggplot() +
  geom_boxplot(data = col_behav, aes(x = col_lat, y = col_hrs, group = colony),
               position = position_dodge(width = 0.5), width = 0.5) +
  geom_ribbon(data = colHrs_glmm.df, aes(y = predicted, x = col_lat,
                                           ymin = conf.low, ymax = conf.high),
              alpha = 0.5, fill = "grey") +
  geom_line(data = colHrs_glmm.df, aes(y = predicted, x = col_lat),
            linewidth = 1, col = "darkorange") +
  scale_y_continuous(breaks = c(0, 6, 12, 18, 24), limits = c(0, 24)) + 
  labs(y = "Hours spent at colony", x = "Colony latitude (°N)") +
  BLKI_theme


#### ~ FIGURE S11 ~ Effects of colony latitude on hours colony attendance ------

png("Figures/FigureS11_hours_attendance.png", width = 9, height = 7, units = "in", res = 600)
colHrs_glmm.plot
dev.off()




