library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(viridisLite)
library(vegan)
library(factoextra)
library(cluster)

# 1. Load the data --------------------------------------------------------

fatty_acid <- read.csv("../cleaned_data/fatty_acid_20190320.csv", header = TRUE)

metadata <- read.csv("../cleaned_data/metadata_20190320.csv", header = TRUE)
distance <- read.csv("../cleaned_data/distance_20190320.csv", header = TRUE)
nutrients <- read.csv("../cleaned_data/nutrients_20190320.csv", header = TRUE)

metadata_dist <- full_join(metadata, distance)

ppcp <- read.csv("../cleaned_data/ppcp_20190320.csv", header = TRUE)

ppcp_meta_dist <- full_join(ppcp, metadata_dist) %>% 
  mutate(POPULATION_INTENSITY = (SOUTH_SHORE_DIST*(POP_SOUTH_DEV/SOUTH_SHORE_AREA))/SOUTH_DEV_DIST) 

fatty_acid_ppcp_meta_dist <- inner_join(fatty_acid, ppcp_meta_dist) %>%
  unite(Taxon, c("Genus", "Species"), remove = FALSE)


# 2. Overall Fatty Acid Analysis ------------------------------------------

## Entire Fatty Acid 

fatty_acid_whole_wide <- fatty_acid %>%
  filter(Genus != "Hyalella") %>%
  select(Site:C24.0) %>%
  gather(Fatty_Acid, Concentration, C12.0:C24.0) %>%
  group_by(Site, Genus, Species) %>%
  mutate(Total_Fatty_Acid = sum(Concentration),
         Prop_Fatty_Acid = Concentration/Total_Fatty_Acid) %>%
  select(-Concentration, - Total_Fatty_Acid) %>%
  spread(Fatty_Acid, Prop_Fatty_Acid) %>%
  unite(Taxon, c("Genus", "Species"), sep = "_", remove = FALSE)

mean <- as.vector(sapply(fatty_acid_whole_wide[,5:63], mean))
var <- as.vector(sapply(fatty_acid_whole_wide[,5:63], var))
mean_var <- data.frame(cbind(mean[1:59], var[1:59]))
colnames(mean_var)[colnames(mean_var) == "X1"] <- "Mean"
colnames(mean_var)[colnames(mean_var) == "X2"] <- "Variance"
mean_var <- dplyr::mutate(mean_var, Var_Mean_RATIO = Variance/Mean)
row.names(mean_var) <- colnames(fatty_acid_whole_wide[,5:63])
mean_var[order(mean_var$Var_Mean_RATIO, decreasing = TRUE), ]

whole_fatty_acid_metaMDS <- metaMDS(fatty_acid_whole_wide[ , 5:63], distance = "bray",
                                         k = 2, try = 100)
whole_fatty_acid_metaMDS

data_scores <- as.data.frame(scores(whole_fatty_acid_metaMDS))
data_scores$Site <- fatty_acid_whole_wide$Site
data_scores$Taxon <- fatty_acid_whole_wide$Taxon
data_scores$Taxon <- gsub("_", " ", data_scores$Taxon)
data_scores$Taxon <- gsub("NA", "", data_scores$Taxon)

species_scores <- as.data.frame(scores(whole_fatty_acid_metaMDS, "species")) 
species_scores$species <- rownames(species_scores)

nmds <- ggplot() + 
  #geom_text(data=species_scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5, size = 5) +  # add the species labels
  geom_point(data=data_scores,aes(x=NMDS1,y=NMDS2, fill = Taxon), size = 10, shape = 21,
             color = "grey60", alpha = .75, stroke = 2) +  # add the site labels
  scale_fill_manual(values = inferno(49)[c(1, 4, 7, 15, 20, 27, 33, 38, 45)]) +
  #scale_size_continuous(range = c(5,20)) +
  ggtitle("NMDS with Entire FA Spectrum") +
  #coord_equal() +
  annotate("label", x = 0.4, y = 0.45, label = paste("Stress: ", round(whole_fatty_acid_metaMDS$stress, digits = 3)),
           size = 10) +
  theme_minimal() +
  theme(legend.position = "bottom",
        title = element_text(size = 20),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18), 
        legend.text = element_text(size = 12))
nmds

ggsave("All_species_all_FA.png", nmds, device = "png",
       path = "../figures/", width = 12, height = 10, units = "in", dpi = 300)

## Essential Fatty Acids 

fatty_acid_essential_wide <- fatty_acid %>%
  filter(Genus != "Hyalella") %>%
  select(Site, Genus, Species, C18.3w3, C18.4w3, C20.5w3, C22.5w3, C22.6w3, 
         C18.2w6, C18.2w6t, C20.4w6) %>%
  gather(Fatty_Acid, Concentration, C18.3w3:C20.4w6) %>%
  group_by(Site, Genus, Species) %>%
  mutate(Total_Fatty_Acid = sum(Concentration),
         Prop_Fatty_Acid = Concentration/Total_Fatty_Acid) %>%
  select(-Concentration, - Total_Fatty_Acid) %>%
  spread(Fatty_Acid, Prop_Fatty_Acid) %>%
  unite(Taxon, c("Genus", "Species"), sep = "_", remove = FALSE)

mean <- as.vector(sapply(fatty_acid_essential_wide[,5:12], mean))
var <- as.vector(sapply(fatty_acid_essential_wide[,5:12], var))
mean_var <- data.frame(cbind(mean[1:8], var[1:8]))
colnames(mean_var)[colnames(mean_var) == "X1"] <- "Mean"
colnames(mean_var)[colnames(mean_var) == "X2"] <- "Variance"
mean_var <- dplyr::mutate(mean_var, Var_Mean_RATIO = Variance/Mean)
row.names(mean_var) <- colnames(fatty_acid_essential_wide[,5:12])
mean_var[order(mean_var$Var_Mean_RATIO, decreasing = TRUE), ]

essential_fatty_acid_metaMDS <- metaMDS(fatty_acid_essential_wide[,5:12], distance = "bray",
                                    k = 2, try = 100)
essential_fatty_acid_metaMDS

data_scores <- as.data.frame(scores(essential_fatty_acid_metaMDS))
data_scores$Site <- fatty_acid_essential_wide$Site
data_scores$Taxon <- fatty_acid_essential_wide$Taxon
data_scores$Taxon <- gsub("_", " ", data_scores$Taxon)
data_scores$Taxon <- gsub("NA", "", data_scores$Taxon)

species_scores <- as.data.frame(scores(essential_fatty_acid_metaMDS, "species")) 
species_scores$species <- rownames(species_scores)

nmds <- ggplot() + 
  #geom_text(data=species_scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5, size = 7) +  # add the species labels
  geom_point(data=data_scores,aes(x=NMDS1,y=NMDS2, fill = Taxon), size = 10, shape = 21,
             color = "grey60", alpha = .75, stroke = 2) +  # add the site labels
  scale_fill_manual(values = inferno(49)[c(1, 4, 7, 15, 20, 27, 33, 38, 45)]) +
  #scale_size_continuous(range = c(5,20)) +
  ggtitle("NMDS with Essential FA Spectrum") +
  #coord_equal() +
  annotate("label", x = -0.5, y = 0.35, label = paste("Stress: ", round(essential_fatty_acid_metaMDS$stress, digits = 3)), size = 10) +
  theme_minimal() +
  theme(legend.position = "bottom",
        title = element_text(size = 20),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18), 
        legend.text = element_text(size = 12))
nmds

ggsave("All_species_Essential_FA.png", nmds, device = "png",
       path = "../figures/", width = 12, height = 10, units = "in",
       dpi = 300)

# 3. Periphyton Fatty Acid Analysis ---------------------------------------

## All Fatty Acids Considered

fatty_acid_wide <- fatty_acid %>%
  select(Site:C24.0) %>%
  filter(Genus %in% c("Drapa", "Splash")) %>%
  unite(Taxon, c("Genus", "Species"), remove = FALSE) %>%
  gather(Fatty_Acid, Concentration, C12.0:C24.0) %>%
  group_by(Site, Genus, Species) %>%
  mutate(Total_Fatty_Acid = sum(Concentration),
         Prop_Fatty_Acid = Concentration/Total_Fatty_Acid) %>%
  select(-Concentration, - Total_Fatty_Acid) %>%
  spread(Fatty_Acid, Prop_Fatty_Acid)
  
mean <- as.vector(sapply(fatty_acid_wide[,5:63], mean))
# Creates a vector of interspecific variances 
var <- as.vector(sapply(fatty_acid_wide[,5:63], var))
#Combines those mean and variance into one dataframe
mean_var <- data.frame(cbind(mean[1:59], var[1:59]))
# Renames columns as mean and variance
colnames(mean_var)[colnames(mean_var) == "X1"] <- "Mean"
colnames(mean_var)[colnames(mean_var) == "X2"] <- "Variance"
# Calculates the variance:mean ratio
mean_var <- dplyr::mutate(mean_var, Var_Mean_RATIO = Variance/Mean)
# Adds rownames, such that row names are protein families
row.names(mean_var) <- colnames(fatty_acid[,5:63])
mean_var[order(mean_var$Var_Mean_RATIO, decreasing = TRUE), ]

periphyton_fatty_acid_metaMDS <- metaMDS(fatty_acid_wide[ , 5:63], distance = "bray",
                                         k = 2, try = 100)

data_scores <- as.data.frame(scores(periphyton_fatty_acid_metaMDS))
data_scores$Site <- fatty_acid_wide$Site
data_scores$Taxon <- fatty_acid_wide$Taxon

species_scores <- as.data.frame(scores(periphyton_fatty_acid_metaMDS, "species")) 
species_scores$species <- rownames(species_scores)

nmds <- ggplot() + 
  geom_text(data=species_scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_text(data=data_scores,aes(x=NMDS1,y=NMDS2,label=Site, color = Taxon), size = 10) +  # add the site labels
  scale_color_manual(values = inferno(36)[c(1, 27, 33)]) +
  #scale_size_continuous(range = c(5,20)) +
  ggtitle("Periphyton NMDS with Entire FA Spectrum") +
  #coord_equal() +
  annotate("label", x = -0.6, y = 0.6, label = paste("Stress: ", round(periphyton_fatty_acid_metaMDS$stress, digits = 3))) +
  theme_minimal() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18), 
        legend.text = element_text(size = 18))
nmds

ggsave("Periphyton_All_FA_20190404.png", nmds, device = "png",
       path = "../figures/", width = 18, height = 18, units = "in")

## Essential Fatty Acids 

fatty_acid_essential_wide <- fatty_acid %>%
  filter(Genus %in% c("Drapa", "Splash")) %>%
  select(Site, Genus, Species, C18.3w3, C18.4w3, C20.5w3, C22.5w3, C22.6w3, 
         C18.2w6, C18.2w6t, C20.4w6) %>%
  gather(Fatty_Acid, Concentration, C18.3w3:C20.4w6) %>%
  group_by(Site, Genus, Species) %>%
  mutate(Total_Fatty_Acid = sum(Concentration),
         Prop_Fatty_Acid = Concentration/Total_Fatty_Acid) %>%
  select(-Concentration, - Total_Fatty_Acid) %>%
  spread(Fatty_Acid, Prop_Fatty_Acid) %>%
  unite(Taxon, c("Genus", "Species"), sep = "_", remove = FALSE)

mean <- as.vector(sapply(fatty_acid_essential_wide[,5:12], mean))
var <- as.vector(sapply(fatty_acid_essential_wide[,5:12], var))
mean_var <- data.frame(cbind(mean[1:8], var[1:8]))
colnames(mean_var)[colnames(mean_var) == "X1"] <- "Mean"
colnames(mean_var)[colnames(mean_var) == "X2"] <- "Variance"
mean_var <- dplyr::mutate(mean_var, Var_Mean_RATIO = Variance/Mean)
row.names(mean_var) <- colnames(fatty_acid_essential_wide[,5:12])
mean_var[order(mean_var$Var_Mean_RATIO, decreasing = TRUE), ]

essential_fatty_acid_metaMDS <- metaMDS(fatty_acid_essential_wide[,5:12], distance = "bray",
                                        k = 2, try = 100)
essential_fatty_acid_metaMDS

data_scores <- as.data.frame(scores(essential_fatty_acid_metaMDS))
data_scores$Site <- fatty_acid_essential_wide$Site
data_scores$Taxon <- fatty_acid_essential_wide$Taxon

species_scores <- as.data.frame(scores(essential_fatty_acid_metaMDS, "species")) 
species_scores$species <- rownames(species_scores)

nmds <- ggplot() + 
  geom_text(data=species_scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_text(data=data_scores,aes(x=NMDS1,y=NMDS2,label=Site, color = Taxon), size = 10) +  # add the site labels
  scale_color_manual(values = inferno(36)[c(1, 27, 33)]) +
  #scale_size_continuous(range = c(5,20)) +
  ggtitle("Periphyton NMDS with Essential FA Spectrum") +
  #coord_equal() +
  annotate("label", x = -0.5, y = -0.002, label = paste("Stress: ", round(essential_fatty_acid_metaMDS$stress, digits = 3))) +
  theme_minimal() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18), 
        legend.text = element_text(size = 18))
nmds

ggsave("Periphyton_Essential_FA_20190404.png", nmds, device = "png",
       path = "../figures/", width = 18, height = 18, units = "in")


# 4. Invertebrate Fatty Acids ---------------------------------------------

## All Fatty Acids Considered

fatty_acid_wide <- fatty_acid %>%
  select(Site:C24.0) %>%
  filter(!(Genus %in% c("Drapa", "Splash", "Hyalella", "Snail")))%>%
  gather(Fatty_Acid, Concentration, C12.0:C24.0) %>%
  group_by(Site, Genus, Species) %>%
  mutate(Total_Fatty_Acid = sum(Concentration),
         Prop_Fatty_Acid = Concentration/Total_Fatty_Acid) %>%
  select(-Concentration, - Total_Fatty_Acid) %>%
  spread(Fatty_Acid, Prop_Fatty_Acid) %>%
  unite(Taxon, c("Genus", "Species"), remove = FALSE)

mean <- as.vector(sapply(fatty_acid_wide[,4:62], mean))
# Creates a vector of interspecific variances 
var <- as.vector(sapply(fatty_acid_wide[,4:62], var))
#Combines those mean and variance into one dataframe
mean_var <- data.frame(cbind(mean[1:59], var[1:59]))
# Renames columns as mean and variance
colnames(mean_var)[colnames(mean_var) == "X1"] <- "Mean"
colnames(mean_var)[colnames(mean_var) == "X2"] <- "Variance"
# Calculates the variance:mean ratio
mean_var <- dplyr::mutate(mean_var, Var_Mean_RATIO = Variance/Mean)
# Adds rownames, such that row names are protein families
row.names(mean_var) <- colnames(fatty_acid[,4:62])
mean_var[order(mean_var$Var_Mean_RATIO, decreasing = TRUE), ]

invert_fatty_acid_metaMDS <- metaMDS(fatty_acid_wide[ , 5:63], distance = "bray",
                                         k = 2, try = 100)
invert_fatty_acid_metaMDS

data_scores <- as.data.frame(scores(invert_fatty_acid_metaMDS))
data_scores$Site <- fatty_acid_wide$Site
data_scores$Taxon <- fatty_acid_wide$Taxon

species_scores <- as.data.frame(scores(invert_fatty_acid_metaMDS, "species")) 
species_scores$species <- rownames(species_scores)

nmds <- ggplot() + 
  geom_text(data=species_scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_text(data=data_scores,aes(x=NMDS1,y=NMDS2,label=Site, color = Taxon), size = 10) +  # add the site labels
  scale_color_manual(values = inferno(36)[c(1, 15, 22, 33)]) +
  #scale_size_continuous(range = c(5,20)) +
  ggtitle("NMDS for Amphipods with Entire FA Spectrum") +
  #coord_equal() +
  annotate("label", x = 0.4, y = 0.35, label = paste("Stress: ", round(invert_fatty_acid_metaMDS$stress, digits = 3))) +
  theme_minimal() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18), 
        legend.text = element_text(size = 18))
nmds
ggsave("Amphipod_All_FA_20190404.png", nmds, device = "png",
       path = "../figures/", width = 18, height = 18, units = "in")

## Essential Fatty Acids 

fatty_acid_essential_wide <- fatty_acid %>%
  filter(!(Genus %in% c("Drapa", "Splash", "Hyalella", "Snail")))%>%
  select(Site, Genus, Species, C18.3w3, C18.4w3, C20.5w3, C22.5w3, C22.6w3, 
         C18.2w6, C18.2w6t, C20.4w6) %>%
  gather(Fatty_Acid, Concentration, C18.3w3:C20.4w6) %>%
  group_by(Site, Genus, Species) %>%
  mutate(Total_Fatty_Acid = sum(Concentration),
         Prop_Fatty_Acid = Concentration/Total_Fatty_Acid) %>%
  select(-Concentration, - Total_Fatty_Acid) %>%
  spread(Fatty_Acid, Prop_Fatty_Acid) %>%
  unite(Taxon, c("Genus", "Species"), sep = "_", remove = FALSE)

mean <- as.vector(sapply(fatty_acid_essential_wide[,5:12], mean))
var <- as.vector(sapply(fatty_acid_essential_wide[,5:12], var))
mean_var <- data.frame(cbind(mean[1:8], var[1:8]))
colnames(mean_var)[colnames(mean_var) == "X1"] <- "Mean"
colnames(mean_var)[colnames(mean_var) == "X2"] <- "Variance"
mean_var <- dplyr::mutate(mean_var, Var_Mean_RATIO = Variance/Mean)
row.names(mean_var) <- colnames(fatty_acid_essential_wide[,5:12])
mean_var[order(mean_var$Var_Mean_RATIO, decreasing = TRUE), ]

essential_fatty_acid_metaMDS <- metaMDS(fatty_acid_essential_wide[,5:12], distance = "bray",
                                        k = 2, try = 100)
essential_fatty_acid_metaMDS

data_scores <- as.data.frame(scores(essential_fatty_acid_metaMDS))
data_scores$Site <- fatty_acid_essential_wide$Site
data_scores$Taxon <- fatty_acid_essential_wide$Taxon

species_scores <- as.data.frame(scores(essential_fatty_acid_metaMDS, "species")) 
species_scores$species <- rownames(species_scores)

nmds <- ggplot() + 
  geom_text(data=species_scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_text(data=data_scores,aes(x=NMDS1,y=NMDS2,label=Site, color = Taxon), size = 10) +  # add the site labels
  scale_color_manual(values = inferno(36)[c(1,  15, 22,  33)]) +
  #scale_size_continuous(range = c(5,20)) +
  ggtitle("NMDS with Essential FA Spectrum") +
  #coord_equal() +
  annotate("label", x = 0.0, y = 0.1, label = paste("Stress: ", round(essential_fatty_acid_metaMDS$stress, digits = 3))) +
  theme_minimal() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18), 
        legend.text = element_text(size = 18))
nmds

ggsave("Amphipod_Essential_FA_20190404.png", nmds, device = "png",
       path = "../figures/", width = 18, height = 18, units = "in")

# 5. Correlating fatty acids with sewage ----------------------------------

fatty_acid__prop_ppcp_meta_dist <- fatty_acid_ppcp_meta_dist %>%
  select(Site, Genus, Species, Taxon, C18.3w3, C18.4w3, C20.5w3, C22.5w3, C22.6w3, 
         C18.2w6, C18.2w6t, C20.4w6, Caffeine:PPCP.SUM, POPULATION_INTENSITY) %>%
  gather(Fatty_acid, Concentration, C18.3w3:C20.4w6) %>%
  group_by(Site, Genus, Species) %>%
  mutate(Total_Fatty_Acid = sum(Concentration),
         Prop_Fatty_Acid = Concentration/Total_Fatty_Acid) %>%
  select(-Concentration, - Total_Fatty_Acid) %>%
  spread(Fatty_acid, Prop_Fatty_Acid) 

ppcp_fa_plot <- fatty_acid__prop_ppcp_meta_dist %>%
  filter(Genus %in% c("Eulimnogammarus", "Splash"),
         Species != "cyaneus") %>%
  mutate(Taxon = ifelse(Taxon == "Eulimnogammarus_verrucosus", "Eulimnogammarus verrucosus", Taxon),
         Taxon = ifelse(Taxon == "Eulimnogammarus_vitatus", "Eulimnogammarus vitatus", Taxon),
         Taxon = ifelse(Taxon == "Splash_zone", "Periphyton", Taxon)) %>%
ggplot(aes(log10(PPCP.SUM), ((C18.3w3+C18.4w3)/(C20.5w3+C22.5w3)))) +
  geom_point(size = 3) +
  facet_wrap(~ Taxon) +
  geom_smooth(method = "lm") +
  xlab("log10([Total PPCP])") +
  ylab(expression(frac(18:3~omega~3 + 18:4~omega~3, 20:5~omega~3 + 22:6~omega~3))) +
  theme_bw() +
  theme(text = element_text(size = 20))

ggsave("../figures/ppcp_fa_plot.png", ppcp_fa_plot, device = "png",
       path = "../figures/", width = 12, height = 6, units = "in")

peri_ppcp_lm <- lm(((C18.3w3+C18.4w3)/(C20.5w3+C22.5w3+C22.6w3)) ~ log10(PPCP.SUM),
               data = fatty_acid__prop_ppcp_meta_dist[fatty_acid__prop_ppcp_meta_dist$Genus == "Splash" , ])
summary(peri_ppcp_lm)

fa_prop_nutrients_meta_dist <- inner_join(fatty_acid__prop_ppcp_meta_dist, nutrients)

peri_phos_lm <- lm(((C18.3w3+C18.4w3)/(C20.5w3+C22.5w3+C22.6w3)) ~ log10(mean_TPO43_mg_dm3),
              data = fa_prop_nutrients_meta_dist[fa_prop_nutrients_meta_dist$Genus == "Splash", ])
summary(peri_phos_lm)

invert_ppcp_lm <- lm(((C18.3w3+C18.4w3)/(C20.5w3+C22.5w3+C22.6w3)) ~ log10(PPCP.SUM),
                   data = fatty_acid__prop_ppcp_meta_dist[fatty_acid__prop_ppcp_meta_dist$Genus != "Splash", ])
summary(invert_ppcp_lm)
