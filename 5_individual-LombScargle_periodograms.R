## ---------------------------
##
## Script name: 5_individual-LombScargle_periodograms.R
##
## Purpose of script: Analyse nest visitation using attendance datasets
##
## Author: Dr. N. P. Huffeldt, adapted by N Gillies
##
## Created: 2024-04-04
## 
## ---------------------------


# Preamble =====================================================================

# Define the packages
packages <- c("tidyverse", "lomb", "plyr", "dplyr", "magrittr", "ggridges", "lme4",
              "lmerTest", "ggeffects")

# Install packages not yet installed - change lib to library path
# installed_packages <- packages %in% rownames(installed.packages())
# 
#  if (any(installed_packages == FALSE)) {
#   install.packages(packages[!installed_packages])
#  }

# Load packages and functions
invisible(lapply(packages, library, character.only = TRUE))

source("BLKI_functions.R")

# Load data --------------------------------------------------------------------

load("Data_inputs/BLKI_attendance_2013-2022.RData")
kits.att_all %<>% select(c(ring, datetime, colony, col_lat, total_imm, hour))

## Count samples per individual - only 2 with < 48 hours
kits.att_all %>% group_by(ring) %>% dplyr::summarise(nrecs = n()) %>% filter(nrecs < 48)

### Remove two short individuals
kits.att_all <- subset(kits.att_all, ring != "NOS_6188689" & ring != "NOS_6234065")

# Split dataframe by ID
acc_List <- split(kits.att_all, kits.att_all$ring)
names(acc_List) <- sapply(acc_List, function(x) x$ring[1])
T_names <- names(acc_List)

median_ci.func <- function (values) {
  
  n = length(values)
  h.re = replicate(10000, median(sample(values, n, rep=T)))
  return(quantile(h.re, c(.025,.975)))
  
}


# Periodogram analysis ---------------------------------------------------------
# Makes empty dataframe to be filled within loop
DF_lsp <- data.frame(ID = as.character(), 
			# month and year maybe not relevant for all analyses, but another variable to group results by if relevant
                     #Month_Year = as.character(),
                     #Month = as.character(),
                     Period = as.numeric(),
                     PNmax = as.numeric(),
                     P_value = as.numeric())

# Lomb-Scargle Periodogram analyses handles NAs, but the input data must be a two column dataframe with date and 
# variable of interest. To capture a rhythm at least double the period length is needed. For 24 h rhythm that 
# means a minimum of 48 h of sampling. Also, at minimum, the sampling frequency has to be half the rate or higher 
# (more frequent) than the period of the rhythm (i.e. for a 24 h rhythm, sampling has to be at least every 12 h).
# read the documentation for lsp and randlsp. randlsp can take a long time to run if many data are included 
# in a dataframe. lsp works well and is a fine alternative for large datasets. periodicities don't change much, 
# but p values are more precise with randlsp.
# In the below code I have "to = 48" because the input data have a sampling fequency of 1 h, and I'm interest 
# in perdiocities up to 48 h. If your data are collected at a higher or lower frequency this needs to be adjusted accordingly

for(i in 1:length(acc_List)){
  pb <- txtProgressBar(min = 0, max = length(acc_List), style = 3)
  setTxtProgressBar(pb, i)
	
  # Subset dataframe to include only two relevant columns, NAs OK
  TS <- acc_List[[i]][,c("datetime","total_imm")]
  res_lsp <- lsp(repeats = 1000, TS$total_imm, to = 48, type = "period", ofac = 5, plot = F)
  
  # run loop through results to obtain the frequency with highest peak/test statistic and the associated test statistic and p value

  DF_lsp_summary <- as.data.frame(summary(res_lsp))
    
  DF_lsp_indv <- data.frame(ID = as.character(names(acc_List[i])), 
		                          Period = as.numeric(DF_lsp_summary[10,]),
                              PNmax = as.numeric(DF_lsp_summary[9,]),  
                              P_value = as.numeric(DF_lsp_summary[11,]))
   DF_lsp <-  rbind(DF_lsp, DF_lsp_indv)
    
}

close(pb)

#save(DF_lsp, file = "Data_outputs/individual_periodograms.RData")
load("Data_outputs/individual_periodograms.RData")

# Summarise period lengths -----------------------------------------------------

### Integrate meta to get colonies ---------------------------------------------

load("Data_outputs/individual_periodograms.RData")
load("Data_inputs/BLKI_metadata_anon.RData")
load("rings_anonymised.RData")

meta$col_lat[meta$colony == "Witless Bay"] <- 47.20

DF_lsp %<>% dplyr::rename(ring = ID)
DF_lsp <- merge(DF_lsp, lookup_table, by = "ring", all.x = T) %>%
  select(-ring) %>%
  mutate(ring = anon_rings) %>%
  select(-anon_rings) %>%
  relocate(ring, .before = Period)

DF_lsp <- merge(DF_lsp, meta, by = "ring")


### Median/mean and standard error by colony -----------------------------------

LSP_Results_All <- DF_lsp %>%
  dplyr::group_by(colony, col_lat) %>%
  dplyr::summarise(Period_mean = mean(Period),
                   Period_median = median(Period),
                   Period_sd = sd(Period),
                   PNmax_mean = mean(PNmax),
                   PNmax_sd = sd(PNmax),
                   P_value_mean = mean(P_value),
                   P_value_sd = sd(P_value),
                   N = n(),
                   mad_med = mad(Period, constant = 1)) %>%
  mutate(Period_se = Period_sd / sqrt(N),
         Period_lower_ci = Period_mean - qt(1 - (0.05 / 2), N - 1) * Period_se,
         Period_upper_ci = Period_mean + qt(1 - (0.05 / 2), N - 1) * Period_se,
         Period_lower_mad = Period_median - mad_med / sqrt(N),
         Period_upper_mad = Period_median + mad_med / sqrt(N)) %>%
  arrange(col_lat)

# Significant periods only
LSP_Results_All.signif <- signif_lsp %>%
  dplyr::group_by(colony, col_lat) %>%
  dplyr::summarise(Period_mean = mean(Period),
                   Period_median = median(Period),
                   Period_sd = sd(Period),
                   PNmax_mean = mean(PNmax),
                   PNmax_sd = sd(PNmax),
                   P_value_mean = mean(P_value),
                   P_value_sd = sd(P_value),
                   N = n(),
                   mad_med = mad(Period, constant = 1)) %>%
  mutate(Period_se = Period_sd / sqrt(N),
         Period_lower_ci = Period_mean - qt(1 - (0.05 / 2), N - 1) * Period_se,
         Period_upper_ci = Period_mean + qt(1 - (0.05 / 2), N - 1) * Period_se,
         Period_lower_mad = Period_median - mad_med / sqrt(N),
         Period_upper_mad = Period_median + mad_med / sqrt(N)) %>%
  arrange(col_lat)


# Analysis ---------------------------------------------------------------------

### Proportion of rhythmic vs arrhythmic individuals -----------------------------

# Overall signif vs non-signif
nrow(subset(DF_lsp, P_value < 0.05))/nrow(DF_lsp)
nrow(subset(DF_lsp, P_value >= 0.05))/nrow(DF_lsp)

# 0.87 vs 0.13 

# Proportion significant for each colony
signif_period <- DF_lsp %>%
  dplyr::group_by(colony, col_lat) %>%
  dplyr::summarise(signif = sum(P_value < 0.05),
                   nonsignif = sum(P_value >= 0.05),
                   prop_signif = sum(P_value < 0.05)/n_distinct(ring)) %>%
  arrange(col_lat)

# Model prop significant against colony
signif.glm <- glm(cbind(signif, nonsignif) ~ col_lat,
                  data = signif_period, 
                  family = "binomial")

plot(signif.glm)
summary(signif.glm)
median(signif_period$prop_signif)

## Plot
effects.signif <- data.frame(ggeffects::ggpredict(signif.glm, terms = "col_lat")) %>% 
  dplyr::rename(col_lat = x)

plot.signif <- ggplot() +
  geom_point(data = signif_period, aes(x = col_lat, y = prop_signif)) +
  geom_line(data = effects.signif, aes(x = col_lat, y = predicted),
            linewidth = 1, col = "darkorange") +
  geom_ribbon(data = effects.signif, (aes(x = col_lat, ymin = conf.low, ymax = conf.high)),
              alpha = 0.5, fill = "grey") + 
  labs(y = "Proportion of birds with significant period lengths",
       x = "Colony latitude °N") +
  scale_x_continuous(breaks = seq(45, 80, by = 5)) + 
  BLKI_theme


### ~ FIGURE SX ~ Prop significant ~ col lat -----------------------------------

png(file = "Figures/FigureS5_prop_signif.png", width = 9, height = 7, units = "in", res = 600)
plot.signif
dev.off()


### Modelling period length with latitude --------------------------------------

period_lmm.signif <- lmer(Period ~ col_lat + (1|colony), data = signif_lsp)
summary(period_lmm)

period_lmm <- lmer(Period ~ col_lat + (1|colony), data = DF_lsp)
summary(period_lmm)

## Plot effects
period_lmm_signif.df <- data.frame(ggpredict(period_lmm.signif, terms = "col_lat")) %>%
  dplyr::rename(col_lat = x)

### Build the plot
period_lmm_signif.plot <- ggplot() +
  geom_point(data = signif_lsp, aes(x = col_lat, y = Period), 
             position = position_jitter(width = 0.25)) +
  geom_ribbon(data = period_lmm_signif.df, aes(y = predicted, x = col_lat,
                                               ymin = conf.low, ymax = conf.high),
              alpha = 0.5, fill = "grey") +
  geom_line(data = period_lmm_signif.df, aes(y = predicted, x = col_lat),
            linewidth = 1, col = "darkorange") +
  labs(x = expression("Colony latitude °N"), y = "Period length (hours)",
       caption = "") +
  scale_x_continuous(breaks = seq(45, 80, by = 5)) + 
  geom_hline(yintercept = 24, linetype = "dashed", col = "navy") +
  BLKI_theme


### ~ FIGURE 3 ~ Period length ~ col lat --------------------------------------

png(file = "Figures/Figure3_period_length.png", width = 9, height = 7, units = "in", res = 600)
period_lmm_signif.plot
dev.off()

# Visualization ----------------------------------------------------------------

## Reorder colony according to latitude
DF_lsp$colony <- as.factor(DF_lsp$colony)
DF_lsp$colony <- reorder(DF_lsp$colony, DF_lsp$col_lat)

### Ridgeplot and median periods for each colony -------------------------------

# Get only significant periods
signif_lsp <- subset(DF_lsp, P_value < 0.05)

# Build plot
period_ridges.signif <- ggplot(signif_lsp, aes(y = colony, x = Period)) + 
  geom_density_ridges(rel_min_height = 0.001, scale = 1.5, aes(fill = col_lat)) +
  geom_vline(xintercept = 24, linetype = "dashed", col = "navy") +
  labs(y = "", x = "Period length (hours)", fill = "Lat °N", tag = "(B)") +
  BLKI_theme +
  theme(legend.justification = "top")


# Get median periods
LSP_Results_All.signif %<>% arrange(col_lat)

median_plot <- ggplot(LSP_Results_All.signif, aes(x = col_lat, y = Period_median)) +
  geom_errorbar(aes(ymin = Period_median - mad_med, ymax = Period_median + mad_med)) +
  geom_point() + 
  labs(y = "Median period ± median absolute deviation (hours)", x = "Colony latitude °N", tag = "(A)") +
  geom_hline(yintercept = 24, linetype = "dashed", col = "navy") +
  geom_vline(xintercept = 66.5, linetype = "dotted", col = "darkorange", size = 1) +
  #geom_text(size = 5, hjust = 1.25) +
  scale_y_continuous(breaks = seq(25, 45, by = 5)) +
  BLKI_theme

### ~ FIGURE 2 ~ Ridgeplot - all periods ---------------------------------------
png(file = "Figures/Figure2_median_period.png", width = 21, height = 9, units = "in", res = 600)
ggpubr::ggarrange(median_plot, period_ridges.signif, nrow = 1,
                  widths = c(0.8, 1))
dev.off()



### Plot again for all periods -------------------------------------------------

period_ridges <- ggplot(DF_lsp, aes(y = colony, x = Period)) + 
  geom_density_ridges(rel_min_height = 0.001, scale = 1.5, aes(fill = col_lat)) +
  geom_vline(xintercept = 24, linetype = "dashed") +
  scale_x_continuous(limits = c(0, 48)) +
  labs(y = "", x = "Period length (hours)", fill = "Lat °N") +
  BLKI_theme +
  theme(legend.justification = "top")


### ~ FIGURE SX ~ Ridgeplot - all periods --------------------------------------

png(file = "Figures/FigureSX_ridgePlot_all.png", width = 9, height = 7, units = "in", res = 600)
period_ridges
dev.off()




