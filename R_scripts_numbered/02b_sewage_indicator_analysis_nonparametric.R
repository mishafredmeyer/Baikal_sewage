# This script is meant to perform univariate analysis of various sewage 
# indicators along Lake Baikal's southwestern shoreline. Each of the sewage
# indicators are regressed against log-transformed inverse distance weighted 
# population using a linear model. For each indicator, an individual plot is
# produced, and then all individual plots are merged into a single, 
# final figure.

# 1. Load packages --------------------------------------------------------

library(tidyverse)
library(viridis)
library(viridisLite)
library(ggpubr)
library(Kendall)

# 2. Load data ------------------------------------------------------------

# PPCP Data
ppcp <- read.csv(file = "../cleaned_data/ppcp.csv",
                 header = TRUE, stringsAsFactors = FALSE)

# Nutrient data
nutrients <- read.csv(file = "../cleaned_data/nutrients.csv",
                      header = TRUE, stringsAsFactors = FALSE)

# Stable isotope data
stable_isotopes <- read.csv(file = "../cleaned_data/stable_isotopes.csv",
                            header = TRUE, stringsAsFactors = FALSE)

# Chlorophyll a data
chlorophylla <- read.csv(file = "../cleaned_data/chlorophylla.csv",
                         header = TRUE, stringsAsFactors = FALSE)

# Microplastics data
microplastics <- read.csv(file = "../cleaned_data/microplastics.csv",
                          header = TRUE, stringsAsFactors = FALSE)

# Site metadata
metadata <- read.csv(file = "../cleaned_data/metadata.csv",
                     header = TRUE, stringsAsFactors = FALSE)

# Site distance data
distance <- read.csv(file = "../cleaned_data/distance_weighted_population_metrics.csv",
                     header = TRUE, stringsAsFactors = FALSE)

# Join site metadata with distance data
metadata_dist <- full_join(x = metadata, y = distance, by = "site")


# 3. PPCP analysis --------------------------------------------------------

# Join PPCP data with metadata/distance and create two custom metrics
ppcp_meta_dist <- full_join(x = ppcp, y = metadata_dist, by = "site")

# Analyze total PPCPs as a function of population intensity
ppcp_PI_model <- Kendall(y = (ppcp_meta_dist$ppcp_sum), 
                         x = (ppcp_meta_dist$distance_weighted_population))

# View model results
summary(ppcp_PI_model)

tau_repo <- rep(x = NA, 1000)

for(i in 1:1000){
  permuted_data <- ppcp_meta_dist %>%
    select(ppcp_sum) %>%
    rename("permuted_ppcp_sum" = "ppcp_sum") %>%
    sample_frac(size = 1) %>%
    cbind(., ppcp_meta_dist)
  
  permuted_model <- Kendall(y = (permuted_data$permuted_ppcp_sum), 
                           x = (permuted_data$distance_weighted_population))
  
  tau_repo[i] <- permuted_model$tau[[1]]
}

tau_repo <- data.frame(tau_repo)

permuted_ppcp_plot <- ggplot() +
  geom_histogram(data = tau_repo, aes(tau_repo), fill = "white",
                 color = "black") +
  geom_vline(aes(xintercept = ppcp_PI_model$tau[[1]])) +
  ggtitle("PPCP vs. IDW Population") +
  ylab("Frequency") +
  xlab("Kendall Tau") +
  theme_minimal()
  
# Export plot
ggsave(filename = "../figures/permuted_ppcp_PI_plot_permuted.png", plot = permuted_ppcp_plot,
       device = "png", width = 18, height = 12, units = "in")


# 4. Nutrient analysis ----------------------------------------------------

# Join nutrient data with metadata/distance and create custom metric
nutrients_meta_dist <- full_join(x = nutrients, y = metadata_dist, by = "site")


# 4.1 Phosphorus ----------------------------------------------------------

# Analyze total Phosphorus as a function of population intensity
phosphorus_PI_model <- Kendall(y = (nutrients_meta_dist$mean_tp_mg_dm3), 
                         x = (nutrients_meta_dist$distance_weighted_population))

# View model results
summary(phosphorus_PI_model)

tau_repo <- rep(x = NA, 1000)

for(i in 1:1000){
  permuted_data <- nutrients_meta_dist %>%
    select(mean_tp_mg_dm3) %>%
    rename("permuted_phosphorus" = "mean_tp_mg_dm3") %>%
    sample_frac(size = 1) %>%
    cbind(., nutrients_meta_dist)
  
  permuted_model <- Kendall(y = log10(permuted_data$permuted_phosphorus), 
                            x = log10(permuted_data$distance_weighted_population))
  
  tau_repo[i] <- permuted_model$tau[[1]]
}

tau_repo <- data.frame(tau_repo)

permuted_phosphorus_plot <- ggplot() +
  geom_histogram(data = tau_repo, aes(tau_repo), fill = "white",
                 color = "black") +
  geom_vline(aes(xintercept = phosphorus_PI_model$tau[[1]])) +
  ggtitle("Phosphorus vs. IDW Population") +
  ylab("Frequency") +
  xlab("Kendall Tau") +
  theme_minimal()

# Export plot
ggsave(filename = "../figures/permuted_phosphorus_PI_plot_permuted.png", 
       plot = permuted_phosphorus_plot,
       device = "png", width = 18, height = 12, units = "in")

# 4.2 Nitrate -------------------------------------------------------------

# Analyze nitrate as a function of population intensity
nitrate_PI_model <- Kendall(y = (nutrients_meta_dist$mean_no3_mg_dm3), 
                               x = (nutrients_meta_dist$distance_weighted_population))

# View model results
summary(nitrate_PI_model)

tau_repo <- rep(x = NA, 1000)

for(i in 1:1000){
  permuted_data <- nutrients_meta_dist %>%
    select(mean_no3_mg_dm3) %>%
    rename("permuted_nitrate" = "mean_no3_mg_dm3") %>%
    sample_frac(size = 1) %>%
    cbind(., nutrients_meta_dist)
  
  permuted_model <- Kendall(y = (permuted_data$permuted_nitrate), 
                            x = (permuted_data$distance_weighted_population))
  
  tau_repo[i] <- permuted_model$tau[[1]]
}

tau_repo <- data.frame(tau_repo)

permuted_nitrate_plot <- ggplot() +
  geom_histogram(data = tau_repo, aes(tau_repo), fill = "white",
                 color = "black") +
  geom_vline(aes(xintercept = nitrate_PI_model$tau[[1]])) +
  ggtitle("Nitrate vs. IDW Population") +
  ylab("Frequency") +
  xlab("Kendall Tau") +
  theme_minimal()

# Export plot
ggsave(filename = "../figures/permuted_nitrate_PI_plot_permuted.png", 
       plot = permuted_nitrate_plot,
       device = "png", width = 18, height = 12, units = "in")


# 4.3 Ammonium ------------------------------------------------------------

# Analyze ammonium as a function of population intensity
ammonium_PI_model <- Kendall(y = (nutrients_meta_dist$mean_nh4_mg_dm3), 
                            x = (nutrients_meta_dist$distance_weighted_population))

# View model results
summary(ammonium_PI_model)

tau_repo <- rep(x = NA, 1000)

for(i in 1:1000){
  permuted_data <- nutrients_meta_dist %>%
    select(mean_nh4_mg_dm3) %>%
    rename("permuted_ammonium" = "mean_nh4_mg_dm3") %>%
    sample_frac(size = 1) %>%
    cbind(., nutrients_meta_dist)
  
  permuted_model <- Kendall(y = (permuted_data$permuted_ammonium), 
                            x = (permuted_data$distance_weighted_population))
  
  tau_repo[i] <- permuted_model$tau[[1]]
}

tau_repo <- data.frame(tau_repo)

permuted_ammonium_plot <- ggplot() +
  geom_histogram(data = tau_repo, aes(tau_repo), fill = "white",
                 color = "black") +
  geom_vline(aes(xintercept = ammonium_PI_model$tau[[1]])) +
  ggtitle("Ammonium vs. IDW Population") +
  ylab("Frequency") +
  xlab("Kendall Tau") +
  theme_minimal()

# Export plot
ggsave(filename = "../figures/permuted_ammonium_PI_plot_permuted.png", 
       plot = permuted_ammonium_plot,
       device = "png", width = 18, height = 12, units = "in")


# 5. Stable isotopes analysis ---------------------------------------------

# Join stable isotope data with metadata/distance and create custom metric
stable_isotopes_meta_dist <- full_join(x = stable_isotopes, y = metadata_dist,
                                       by = "site")

# 5.1 N15 -----------------------------------------------------------------

# Analyze N15 as a function of population intensity
n15_PI_model <- Kendall(y = (stable_isotopes_meta_dist$N15), 
                             x = (stable_isotopes_meta_dist$distance_weighted_population))

# View model results
summary(n15_PI_model)

tau_repo <- rep(x = NA, 1000)

for(i in 1:1000){
  permuted_data <- stable_isotopes_meta_dist %>%
    select(N15) %>%
    rename("permuted_N15" = "N15") %>%
    sample_frac(size = 1) %>%
    cbind(., stable_isotopes_meta_dist)
  
  permuted_model <- Kendall(y = (permuted_data$permuted_N15), 
                            x = (permuted_data$distance_weighted_population))
  
  tau_repo[i] <- permuted_model$tau[[1]]
}

tau_repo <- data.frame(tau_repo)

permuted_n15_plot <- ggplot() +
  geom_histogram(data = tau_repo, aes(tau_repo), fill = "white",
                 color = "black") +
  geom_vline(aes(xintercept = n15_PI_model$tau[[1]])) +
  ggtitle("N15 vs. IDW Population") +
  ylab("Frequency") +
  xlab("Kendall Tau") +
  theme_minimal()

# Export plot
ggsave(filename = "../figures/permuted_N15_PI_plot_permuted.png", 
       plot = permuted_n15_plot,
       device = "png", width = 18, height = 12, units = "in")

# 5.2 C13 -----------------------------------------------------------------

# Analyze C13 as a function of population intensity
c13_PI_model <- Kendall(y = (stable_isotopes_meta_dist$C13), 
                        x = (stable_isotopes_meta_dist$distance_weighted_population))

# View model results
summary(c13_PI_model)

tau_repo <- rep(x = NA, 1000)

for(i in 1:1000){
  permuted_data <- stable_isotopes_meta_dist %>%
    select(C13) %>%
    rename("permuted_C13" = "C13") %>%
    sample_frac(size = 1) %>%
    cbind(., stable_isotopes_meta_dist)
  
  permuted_model <- Kendall(y = (permuted_data$permuted_C13), 
                            x = (permuted_data$distance_weighted_population))
  
  tau_repo[i] <- permuted_model$tau[[1]]
}

tau_repo <- data.frame(tau_repo)

permuted_c13_plot <- ggplot() +
  geom_histogram(data = tau_repo, aes(tau_repo), fill = "white",
                 color = "black") +
  geom_vline(aes(xintercept = c13_PI_model$tau[[1]])) +
  ggtitle("C13 vs. IDW Population") +
  ylab("Frequency") +
  xlab("Kendall Tau") +
  theme_minimal()

# Export plot
ggsave(filename = "../figures/permuted_C13_PI_plot_permuted.png", 
       plot = permuted_c13_plot,
       device = "png", width = 18, height = 12, units = "in")


# 6. Chlorophyll a analysis -----------------------------------------------

# Join Chl a data with metadata/distance and create custom metric
chlorophylla_meta_dist <- full_join(x = chlorophylla, y = metadata_dist,
                                    by = "site")

# Analyze chl a as a function of population intensity
chla_PI_model <- Kendall(y = (chlorophylla_meta_dist$mean_chlorophylla), 
                        x = (chlorophylla_meta_dist$distance_weighted_population))

# View model results
summary(c13_PI_model)

tau_repo <- rep(x = NA, 1000)

for(i in 1:1000){
  permuted_data <- chlorophylla_meta_dist %>%
    select(mean_chlorophylla) %>%
    rename("permuted_chla" = "mean_chlorophylla") %>%
    sample_frac(size = 1) %>%
    cbind(., chlorophylla_meta_dist)
  
  permuted_model <- Kendall(y = (permuted_data$permuted_chla), 
                            x = (permuted_data$distance_weighted_population))
  
  tau_repo[i] <- permuted_model$tau[[1]]
}

tau_repo <- data.frame(tau_repo)

permuted_chla_plot <- ggplot() +
  geom_histogram(data = tau_repo, aes(tau_repo), fill = "white",
                 color = "black") +
  geom_vline(aes(xintercept = chla_PI_model$tau[[1]])) +
  ggtitle("Chlorophyll a vs. IDW Population") +
  ylab("Frequency") +
  xlab("Kendall Tau") +
  theme_minimal()

# Export plot
ggsave(filename = "../figures/permuted_chla_PI_plot_permuted.png", 
       plot = permuted_chla_plot,
       device = "png", width = 18, height = 12, units = "in")

# 7. Microplastics analysis -----------------------------------------------

# Format microplastics data before join
microplastics <- microplastics %>%
  group_by(site) %>%
  summarize(mean_total = mean(x = total_microplastics, na.rm = TRUE),
            mean_density = mean(x = density, na.rm = TRUE),
            mean_fragment_density = mean(x = fragment_density, na.rm = TRUE),
            mean_fiber_density = mean(x = fiber_density, na.rm = TRUE),
            mean_bead_density = mean(x = bead_density, na.rm = TRUE))

# Join microplastics data with metadata/distance and create custom metric
microplastics_meta_dist <- full_join(x = microplastics, y = metadata_dist,
                                     by = "site")


# 7.1 Mean total microplastics --------------------------------------------

# Analyze mean total microplastics as a function of population intensity
microplastics_mean_PI_model <- Kendall(y = (microplastics_meta_dist$mean_total), 
                                       x = (microplastics_meta_dist$distance_weighted_population))

# View model results
summary(microplastics_mean_PI_model)

tau_repo <- rep(x = NA, 1000)

for(i in 1:1000){
  permuted_data <- microplastics_meta_dist %>%
    select(mean_total) %>%
    rename("permuted_mp_total" = "mean_total") %>%
    sample_frac(size = 1) %>%
    cbind(., microplastics_meta_dist)
  
  permuted_model <- Kendall(y = (permuted_data$permuted_mp_total), 
                            x = (permuted_data$distance_weighted_population))
  
  tau_repo[i] <- permuted_model$tau[[1]]
}

tau_repo <- data.frame(tau_repo)

permuted_mp_total_plot <- ggplot() +
  geom_histogram(data = tau_repo, aes(tau_repo), fill = "white",
                 color = "black") +
  geom_vline(aes(xintercept = microplastics_mean_PI_model$tau[[1]])) +
  ggtitle("Mean Total Microplastics vs. IDW Population") +
  ylab("Frequency") +
  xlab("Kendall Tau") +
  theme_minimal()

# Export plot
ggsave(filename = "../figures/permuted_mp_total_PI_plot_permuted.png", 
       plot = permuted_mp_total_plot,
       device = "png", width = 18, height = 12, units = "in")


# 7.2 Mean microplastic density -------------------------------------------

# Analyze mean microplastic density as a function of population intensity
microplastics_density_PI_model <- Kendall(y = (microplastics_meta_dist$mean_density), 
                                       x = (microplastics_meta_dist$distance_weighted_population))

# View model results
summary(microplastics_density_PI_model)

tau_repo <- rep(x = NA, 1000)

for(i in 1:1000){
  permuted_data <- microplastics_meta_dist %>%
    select(mean_density) %>%
    rename("permuted_mp_density" = "mean_density") %>%
    sample_frac(size = 1) %>%
    cbind(., microplastics_meta_dist)
  
  permuted_model <- Kendall(y = (permuted_data$permuted_mp_density), 
                            x = (permuted_data$distance_weighted_population))
  
  tau_repo[i] <- permuted_model$tau[[1]]
}

tau_repo <- data.frame(tau_repo)

permuted_mp_density_plot <- ggplot() +
  geom_histogram(data = tau_repo, aes(tau_repo), fill = "white",
                 color = "black") +
  geom_vline(aes(xintercept = microplastics_density_PI_model$tau[[1]])) +
  ggtitle("Mean Microplastics Density vs. IDW Population") +
  ylab("Frequency") +
  xlab("Kendall Tau") +
  theme_minimal()

# Export plot
ggsave(filename = "../figures/permuted_mp_density_PI_plot_permuted.png", 
       plot = permuted_mp_density_plot,
       device = "png", width = 18, height = 12, units = "in")


# 8. Combine plots --------------------------------------------------------

ggarrange(permuted_ppcp_plot, permuted_n15_plot, permuted_phosphorus_plot, permuted_chla_plot,
          permuted_nitrate_plot, permuted_mp_total_plot, permuted_ammonium_plot,
          permuted_mp_density_plot, ncol = 2, nrow = 4, labels = "AUTO") %>%
  ggexport(filename = "../figures/combined_permuted_plot.png",
           height = 1900, width = 1200, res = 120)
