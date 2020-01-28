library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(viridisLite)
library(ggpubr)

# 1. Load data ------------------------------------------------------------

ppcp <- read.csv("../cleaned_data/ppcp_20190320.csv", header = TRUE)
nutrients <- read.csv("../cleaned_data/nutrients_20190320.csv", header = TRUE)
stable_isotopes <- read.csv("../cleaned_data/stable_isotopes_20190320.csv", header = TRUE)
chlorophylla <- read.csv("../cleaned_data/chlorophylla_20190320.csv", header = TRUE)
microplastics <- read.csv("../cleaned_data/microplastics_20190320.csv", header = TRUE)

metadata <- read.csv("../cleaned_data/metadata_20190320.csv", header = TRUE)
distance <- read.csv("../cleaned_data/distance_20190320.csv", header = TRUE)

metadata_dist <- full_join(metadata, distance)


# 2. PPCP analysis --------------------------------------------------------

ppcp_meta_dist <- full_join(ppcp, metadata_dist) %>%
  mutate(POPULATION_INTENSITY = (SOUTH_SHORE_DIST*(POP_SOUTH_DEV/SOUTH_SHORE_AREA))/SOUTH_DEV_DIST,
         TOTAL_PPCP = Acetaminophen+Paraxanthine+Cotinine+Caffeine)

ppcp_PI_model <- lm(log10(TOTAL_PPCP) ~ log10(POPULATION_INTENSITY), data = ppcp_meta_dist)
summary(ppcp_PI_model)

ppcp_PI_plot <- ggplot(ppcp_meta_dist, aes(log10(POPULATION_INTENSITY), log10(TOTAL_PPCP))) +
 geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  ylab("log10([Total PPCP])") +
  xlab(expression(paste("log10(", over("Population density at nearest development", "Distance to nearest development"), ")"))) +
  ggtitle("PPCP vs. Distance-weighted Population") +
  # annotate("label",x = 1.75, y = -1.33,
  #          label = deparse(bquote(~p-value: ~  .(round(summary(ppcp_PI_model)$coefficients[2,4], 3)) ~ R^2: ~ .(round(summary(ppcp_PI_model)$r.squared,3)))),
  #          parse = T)
  annotate("label",x = 1.75, y = -1.33,  
           label = paste("p-value: ", round(summary(ppcp_PI_model)$coefficients[2,4], 3), 
                          "\nR-squared: ", round(summary(ppcp_PI_model)$r.squared,3)), parse = TRUE) +
  theme_minimal() #+ 
  # theme(title = element_text(size = 30), 
  #       axis.text = element_text(size = 20),
  #       axis.title = element_text(size = 30),
  #       axis.title.y=element_text(margin=margin(0,20,0,0)), 
  #       axis.title.x=element_text(margin=margin(20,0,0,0)))

ggsave("ppcp_PI_plot.png", ppcp_PI_plot, device = "png",
       path = "../figures", width = 18, height = 12, units = "in")


# Nutrient Analysis -----------------------------------------------------

##First combine the data 
nutrients_meta_dist <- full_join(nutrients, metadata_dist) %>%
  mutate(POPULATION_INTENSITY = (SOUTH_SHORE_DIST*(POP_SOUTH_DEV/SOUTH_SHORE_AREA))/SOUTH_DEV_DIST)

##Analyze phosphorus as a function of population intensity
phosphorus_PI_model <- lm(log10(mean_TP_mg_dm3) ~ log10(POPULATION_INTENSITY), data = nutrients_meta_dist)
summary(phosphorus_PI_model)

phosphorus_PI_plot <- ggplot(nutrients_meta_dist, aes(log10(POPULATION_INTENSITY), log10(mean_TP_mg_dm3))) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  ylab("log10([Total Phosphorus])") +
  xlab(expression(paste("log10(", over("Population density at nearest development", "Distance to nearest development"), ")"))) +
  ggtitle("Phosphorus vs. Distance-weighted Population") +
  annotate("label",x = 1.75, y = -1.45,  
           label = paste0("p-value: ", round(summary(phosphorus_PI_model)$coefficients[2,4], 3),
                          "\nr-squared: ", round(summary(phosphorus_PI_model)$r.squared,3))) +
  theme_minimal() #+ 
  # theme(title = element_text(size = 30),
  #       axis.text = element_text(size = 20),
  #       axis.title = element_text(size = 30),
  #       axis.title.y=element_text(margin=margin(0,20,0,0)), 
  #       axis.title.x=element_text(margin=margin(20,0,0,0)))

ggsave("phosphorus_PI_plot.png", phosphorus_PI_plot, device = "png",
       path = "../figures", width = 18, height = 12, units = "in")

##Analyze nitrate
nitrate_PI_model <- lm(log10(mean_NO3_mg_dm3) ~ log10(POPULATION_INTENSITY), data = nutrients_meta_dist)
summary(nitrate_PI_model)

nitrate_PI_plot <- ggplot(nutrients_meta_dist, aes(log10(POPULATION_INTENSITY), log10(mean_NO3_mg_dm3))) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  ylab("log10([Total Nitrate])") +
  xlab(expression(paste("log10(", over("Population density at nearest development", "Distance to nearest development"), ")"))) +
  ggtitle("Nitrate vs. Distance-weighted Population") +
  annotate("label",x = 1.75, y = -0.8,  
           label = paste0("p-value: ", round(summary(nitrate_PI_model)$coefficients[2,4], 3),
                          "\nr-squared: ", round(summary(nitrate_PI_model)$r.squared,3))) +
  theme_minimal() #+ 
  # theme(title = element_text(size = 30),
  #       axis.text = element_text(size = 20),
  #       axis.title = element_text(size = 30),
  #       axis.title.y=element_text(margin=margin(0,20,0,0)), 
  #       axis.title.x=element_text(margin=margin(20,0,0,0)))

ggsave("nitrate_PI_plot.png", nitrate_PI_plot, device = "png",
       path = "../figures", width = 18, height = 12, units = "in")

##Analyze ammonium
ammonium_PI_model <- lm(log10(mean_NH4_mg_dm3) ~ log10(POPULATION_INTENSITY), data = nutrients_meta_dist)
summary(ammonium_PI_model)

ammonium_PI_plot <- ggplot(nutrients_meta_dist, aes(log10(POPULATION_INTENSITY), log10(mean_NH4_mg_dm3))) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  ylab("log10([Total Ammonium])") +
  xlab(expression(paste("log10(", over("Population density at nearest development", "Distance to nearest development"), ")"))) +
  ggtitle("Ammonium vs. Distance-weighted Population") +
  annotate("label",x = 1.75, y = -1.25,  
           label = paste0("p-value: ", round(summary(ammonium_PI_model)$coefficients[2,4], 3),
                          "\nr-squared: ", round(summary(ammonium_PI_model)$r.squared,3))) +
  theme_minimal() #+ 
  # theme(title = element_text(size = 30),
  #       axis.text = element_text(size = 20),
  #       axis.title = element_text(size = 30),
  #       axis.title.y=element_text(margin=margin(0,20,0,0)), 
  #       axis.title.x=element_text(margin=margin(20,0,0,0)))

ggsave("ammonium_PI_plot.png", ammonium_PI_plot, device = "png",
       path = "../figures", width = 18, height = 12, units = "in")


# Stable Isotopes Analysis ------------------------------------------------

stable_isotopes_meta_dist <- full_join(stable_isotopes, metadata_dist) %>%
  mutate(POPULATION_INTENSITY = (SOUTH_SHORE_DIST*(POP_SOUTH_DEV/SOUTH_SHORE_AREA))/SOUTH_DEV_DIST)

##Analyze phosphorus as a function of population intensity
n15_PI_model <- lm(log10(N15) ~ log10(POPULATION_INTENSITY), data = stable_isotopes_meta_dist[stable_isotopes$Genus != "Sp.", ])
summary(n15_PI_model)

n15_PI_plot <- ggplot(stable_isotopes_meta_dist[stable_isotopes$Genus != "Sp.", ], 
                      aes(log10(POPULATION_INTENSITY), log10(N15))) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  ylab("log10(N15)") +
  xlab(expression(paste("log10(", over("Population density at nearest development", "Distance to nearest development"), ")"))) +
  ggtitle("N15 vs. Distance-weighted Population") +
  annotate("label",x = 1.75, y = 0.89,  
           label = paste0("p-value: ", round(summary(n15_PI_model)$coefficients[2,4], 3),
                          "\nr-squared: ", round(summary(n15_PI_model)$r.squared,3))) +
  theme_minimal() #+ 
  # theme(title = element_text(size = 30),
  #       axis.text = element_text(size = 20),
  #       axis.title = element_text(size = 30),
  #       axis.title.y=element_text(margin=margin(0,20,0,0)), 
  #       axis.title.x=element_text(margin=margin(20,0,0,0)))

ggsave("n15_PI_plot.png", n15_PI_plot, device = "png",
       path = "../figures", width = 18, height = 12, units = "in")

##Analyze phosphorus as a function of population intensity
c13_PI_model <- lm((C13) ~ log10(POPULATION_INTENSITY), data = stable_isotopes_meta_dist[stable_isotopes$Genus != "Sp.", ])
summary(c13_PI_model)

c13_PI_plot <- ggplot(stable_isotopes_meta_dist[stable_isotopes$Genus != "Sp.", ], 
                      aes(log10(POPULATION_INTENSITY), (C13))) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  ylab("C13%") +
  xlab(expression(paste("log10(", over("Population density at nearest development", "Distance to nearest development"), ")"))) +
  ggtitle("C13% vs. Distance-weighted Population") +
  annotate("label",x = 1.75, y = -17,  
           label = paste0("p-value: ", round(summary(c13_PI_model)$coefficients[2,4], 3),
                          "\nr-squared: ", round(summary(c13_PI_model)$r.squared,3))) +
  theme_minimal() #+ 
  # theme(title = element_text(size = 30),
  #       axis.text = element_text(size = 20),
  #       axis.title = element_text(size = 30),
  #       axis.title.y=element_text(margin=margin(0,20,0,0)), 
  #       axis.title.x=element_text(margin=margin(20,0,0,0)))

ggsave("c13_PI_plot.png", c13_PI_plot, device = "png",
       path = "../figures", width = 18, height = 12, units = "in")


# Chlorophyll a analysis --------------------------------------------------

chlorophylla_meta_dist <- full_join(chlorophylla, metadata_dist) %>%
  mutate(POPULATION_INTENSITY = (SOUTH_SHORE_DIST*(POP_SOUTH_DEV/SOUTH_SHORE_AREA))/SOUTH_DEV_DIST)

chlorophylla_PI_model <- lm((mean_chlorophylla) ~ log10(POPULATION_INTENSITY), data = chlorophylla_meta_dist)
summary(chlorophylla_PI_model)

chlorophylla_PI_plot <- ggplot(chlorophylla_meta_dist,aes(log10(POPULATION_INTENSITY), (mean_chlorophylla))) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  ylab("Chlorophyll a") +
  xlab(expression(paste("log10(", over("Population density at nearest development", "Distance to nearest development"), ")"))) +
  ggtitle("Chlorophyll a vs. Distance-weighted Population") +
  annotate("label",x = 1.75, y = 1.3,  
           label = paste0("p-value: ", round(summary(chlorophylla_PI_model)$coefficients[2,4], 3),
                          "\nr-squared: ", round(summary(chlorophylla_PI_model)$r.squared,3))) +
   theme_minimal() #+
  # theme(title = element_text(size = 30),
  #       axis.text = element_text(size = 20),
  #       axis.title = element_text(size = 30),
  #       axis.title.y=element_text(margin=margin(0,20,0,0)),
  #       axis.title.x=element_text(margin=margin(20,0,0,0)))

ggsave("chlorophylla_PI_plot.png", chlorophylla_PI_plot, device = "png",
       path = "../figures", width = 18, height = 12, units = "in")


# Microplastics Analysis --------------------------------------------------

microplastics_meta_dist <- full_join(microplastics, metadata_dist) %>%
  mutate(POPULATION_INTENSITY = (SOUTH_SHORE_DIST*(POP_SOUTH_DEV/SOUTH_SHORE_AREA))/SOUTH_DEV_DIST)

##Mean Total Microplastics
microplastics_total_PI_model <- lm((mean_total_microplastics) ~ log10(POPULATION_INTENSITY), data = microplastics_meta_dist)
summary(microplastics_total_PI_model)

microplastics_total_PI_plot <- ggplot(microplastics_meta_dist,
                                      aes(log10(POPULATION_INTENSITY), (mean_total_microplastics))) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  ylab("Mean Total Microplastics") +
  xlab(expression(paste("log10(", over("Population density at nearest development", "Distance to nearest development"), ")"))) +
  ggtitle("Total Microplastics vs. Distance-weighted Population") +
  annotate("label",x = 1.75, y = 9,  
           label = paste0("p-value: ", round(summary(microplastics_total_PI_model)$coefficients[2,4], 3),
                          "\nr-squared: ", round(summary(microplastics_total_PI_model)$r.squared,3))) +
  theme_minimal() #+ 
  # theme(title = element_text(size = 30),
  #       axis.text = element_text(size = 20),
  #       axis.title = element_text(size = 30),
  #       axis.title.y=element_text(margin=margin(0,20,0,0)), 
  #       axis.title.x=element_text(margin=margin(20,0,0,0)))

ggsave("microplastics_total_PI_plot.png", microplastics_total_PI_plot, device = "png",
       path = "../figures", width = 18, height = 12, units = "in")

##Mean Microplastic Density
microplastics_density_PI_model <- lm((mean_microplastic_density) ~ log10(POPULATION_INTENSITY), data = microplastics_meta_dist)
summary(microplastics_density_PI_model)

microplastics_density_PI_plot <- ggplot(microplastics_meta_dist,
                                      aes(log10(POPULATION_INTENSITY), (mean_microplastic_density))) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  ylab("Mean Microplastic Density") +
  xlab(expression(paste("log10(", over("Population density at nearest development", "Distance to nearest development"), ")"))) +
  ggtitle("Microplastics Density vs. Distance-weighted Population") +
  annotate("label",x = 1.75, y = 7,
           label = paste0("p-value: ", round(summary(microplastics_density_PI_model)$coefficients[2,4], 3),
                          "\nr-squared: ", round(summary(microplastics_density_PI_model)$r.squared,3))) +
  theme_minimal() #+ 
  # theme(title = element_text(size = 30),
  #       axis.text = element_text(size = 20),
  #       axis.title = element_text(size = 30),
  #       axis.title.y=element_text(margin=margin(0,20,0,0)), 
  #       axis.title.x=element_text(margin=margin(20,0,0,0)))

ggsave("microplastics_density_PI_plot.png", microplastics_density_PI_plot, device = "png",
       path = "../figures", width = 18, height = 12, units = "in")


# Combine Plots -----------------------------------------------------------

ggarrange(ppcp_PI_plot, n15_PI_plot, phosphorus_PI_plot, chlorophylla_PI_plot,
                            nitrate_PI_plot, microplastics_total_PI_plot, ammonium_PI_plot,
                            microplastics_density_PI_plot, ncol = 2, nrow = 4, labels = "AUTO") %>%
  ggexport(filename = "../figures/combined_plot.png", height = 1900, width = 1200, res = 120)
