# This script aggregate the fatty acid data and performs
# several univariate and multivariate analyses aimed at
# relating fatty acid community compositions with sewage
# indicators in Lake Baikal.

library(tidyverse)
library(viridis)
library(viridisLite)
library(vegan)


# 1. Load the data --------------------------------------------------------

fatty_acid <- read.csv(file = "../cleaned_data/fatty_acid.csv",
                       header = TRUE, stringsAsFactors = FALSE)

metadata <- read.csv(file = "../cleaned_data/metadata.csv",
                     header = TRUE, stringsAsFactors = FALSE)

distance <- read.csv(file = "../cleaned_data/distance_weighted_population_metrics.csv",
                     header = TRUE, stringsAsFactors = FALSE)

nutrients <- read.csv(file = "../cleaned_data/nutrients.csv",
                      header = TRUE, stringsAsFactors = FALSE)

metadata_dist <- full_join(x = metadata, y = distance)

ppcp <- read.csv(file = "../cleaned_data/ppcp.csv",
                 header = TRUE, stringsAsFactors = FALSE)

ppcp_meta_dist <- full_join(x = ppcp, y = metadata_dist, by = c("site"))

fatty_acid_ppcp_meta_dist <- inner_join(x = fatty_acid, y = ppcp_meta_dist,
                                        by = c("site")) %>%
  unite(taxon, c("Genus", "Species"), remove = FALSE)


# 2. Overall Fatty Acid Analysis ------------------------------------------


# 2.1 Whole fatty acid spectrum -------------------------------------------

# Create a dataframe of the entire fatty acid spectrum
fatty_acid_whole_wide <- fatty_acid %>%
  filter(Genus != "Hyalella") %>%
  select(site:C24.0) %>%
  gather(key = fatty_acid, value = concentration, C12.0:C24.0) %>%
  group_by(site, Genus, Species) %>%
  mutate(total_fatty_acid = sum(concentration),
         prop_fatty_acid = concentration / total_fatty_acid) %>%
  select(-concentration, - total_fatty_acid) %>%
  spread(key = fatty_acid, value = prop_fatty_acid) %>%
  unite(taxon, c("Genus", "Species"), sep = "_", remove = FALSE)

# Identify mean, variance, and coefficient of variation across all
# sites for each taxonomic grouping
mean <- as.vector(sapply(fatty_acid_whole_wide[, 5:63], mean))
var <- as.vector(sapply(fatty_acid_whole_wide[, 5:63], var))
mean_var <- data.frame(cbind(Mean = mean[1:59], Variance = var[1:59]))
mean_var <- dplyr::mutate(mean_var, Var_Mean_Ratio = Variance / Mean)
row.names(mean_var) <- colnames(fatty_acid_whole_wide[, 5:63])
mean_var[order(mean_var$Var_Mean_Ratio, decreasing = TRUE), ]

# Perform NMDS
whole_fatty_acid_metaMDS <- metaMDS(comm = fatty_acid_whole_wide[, 5:63],
                                    distance = "bray", k = 2, try = 100)
whole_fatty_acid_metaMDS

# Extract scores for each site and add additional data for analysis and plotting
data_scores <- as.data.frame(scores(whole_fatty_acid_metaMDS)) %>%
  mutate(site = fatty_acid_whole_wide$site,
         taxon = fatty_acid_whole_wide$taxon,
         taxon = gsub(pattern = "_", replacement = " ", x = taxon),
         taxon = gsub(pattern = "NA", replacement = "", x = taxon),
         taxon = gsub(pattern = "Drapa", replacement = "Drapa spp.", x = taxon))

# Create plot
# This figure is associated with figure S1 in the associated manuscript
nmds <- ggplot() +
  geom_point(data = data_scores, aes(x = NMDS1, y = NMDS2, fill = taxon),
             size = 10, shape = 21, color = "grey60", alpha = .75, stroke = 2) +
  scale_fill_manual(values = inferno(69)[c(1, 13, 18, 20, 27, 33, 38, 45, 55)]) +
  ggtitle("NMDS with Entire FA Spectrum") +
  annotate("label", x = 0.4, y = 0.45,
           label = paste("Stress: ", round(whole_fatty_acid_metaMDS$stress, digits = 3)),
           size = 10) +
  theme_minimal() +
  theme(legend.position = "bottom",
        title = element_text(size = 20),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 12))
nmds

ggsave(filename = "all_species_all_FA.png", plot = nmds, device = "png",
       path = "../figures/", width = 12, height = 10, units = "in", dpi = 300)


# 2.2 Redo analysis with only essential fatty acids -----------------------

# Create a dataframe of only the essential fatty acids
fatty_acid_essential_wide <- fatty_acid %>%
  filter(Genus != "Hyalella") %>%
  select(site, Genus, Species, C18.3w3, C18.4w3, C20.5w3, C22.5w3, C22.6w3,
         C18.2w6, C18.2w6t, C20.4w6) %>%
  gather(key = fatty_acid, value = concentration, C18.3w3:C20.4w6) %>%
  group_by(site, Genus, Species) %>%
  mutate(total_fatty_acid = sum(concentration),
         prop_fatty_acid = concentration / total_fatty_acid) %>%
  select(-concentration, - total_fatty_acid) %>%
  spread(key = fatty_acid, value = prop_fatty_acid) %>%
  unite(col = taxon, c("Genus", "Species"), sep = "_", remove = FALSE)

# Identify mean, variance, and coefficient of variation across all
# sites for each taxonomic grouping
mean <- as.vector(sapply(fatty_acid_essential_wide[, 5:12], mean))
var <- as.vector(sapply(fatty_acid_essential_wide[, 5:12], var))
mean_var <- data.frame(cbind(Mean = mean[1:8], Variance = var[1:8]))
mean_var <- dplyr::mutate(mean_var, Var_Mean_Ratio = Variance / Mean)
row.names(mean_var) <- colnames(fatty_acid_essential_wide[, 5:12])
mean_var[order(mean_var$Var_Mean_Ratio, decreasing = TRUE), ]

# Perform NMDS
essential_fatty_acid_metaMDS <- metaMDS(comm = fatty_acid_essential_wide[, 5:12],
                                        distance = "bray", k = 2, try = 100)
essential_fatty_acid_metaMDS

# Aggregate data scores for plotting and analysis of the NMDS
data_scores <- as.data.frame(scores(essential_fatty_acid_metaMDS)) %>%
  mutate(site = fatty_acid_essential_wide$site,
         taxon = fatty_acid_essential_wide$taxon,
         taxon = gsub(pattern = "_", replacement = " ", x = taxon),
         taxon = gsub(pattern = "NA", replacement = "", x = taxon),
         taxon = gsub(pattern = "Drapa", replacement = "Drapa spp.", x = taxon))

# Plot NMDS for all species but only Essential Fatty Acids
# This plot is figure S2 in the associated ms.
nmds <- ggplot() +
  geom_point(data = data_scores, aes(x = NMDS1, y = NMDS2, fill = taxon),
             size = 10, shape = 21, color = "grey60", alpha = .75, stroke = 2) +
  scale_fill_manual(values = inferno(69)[c(1, 13, 18, 20, 27, 33, 38, 45, 55)]) +
  ggtitle("NMDS with Essential FA Spectrum") +
  annotate("label", x = -0.5, y = 0.35,
           label = paste("Stress: ",
                         round(essential_fatty_acid_metaMDS$stress, digits = 3)),
           size = 10) +
  theme_minimal() +
  theme(legend.position = "bottom",
        title = element_text(size = 20),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 12))
nmds

ggsave(filename = "all_species_essential_FA.png", plot = nmds, device = "png",
       path = "../figures/", width = 12, height = 10, units = "in",
       dpi = 300)


# 3. Correlating fatty acids with sewage ----------------------------------

fatty_acid_prop_ppcp_meta_dist <- fatty_acid_ppcp_meta_dist %>%
  select(site, taxon, Genus, Species, C18.3w3, C18.4w3, C20.5w3, C22.5w3, C22.6w3,
         C18.2w6, C18.2w6t, C20.4w6, caffeine:ppcp_sum, distance_weighted_population) %>%
  gather(key = fatty_acid, value = concentration, C18.3w3:C20.4w6) %>%
  group_by(site, Genus, Species) %>%
  mutate(total_fatty_acid = sum(concentration),
         prop_fatty_acid = concentration / total_fatty_acid) %>%
  select(-concentration, -total_fatty_acid) %>%
  spread(key = fatty_acid, value = prop_fatty_acid)

# This figure is Figure 7 within the body of the associated manuscript.
ppcp_fa_plot <- fatty_acid_prop_ppcp_meta_dist %>%
  filter(taxon %in% c("Eulimnogammarus_verrucosus", "Eulimnogammarus_vitatus",
                      "Periphyton_NA")) %>%
  mutate(taxon = ifelse(test = taxon == "Eulimnogammarus_verrucosus",
                        yes = "Eulimnogammarus verrucosus", no = taxon),
         taxon = ifelse(test = taxon == "Eulimnogammarus_vitatus",
                        yes = "Eulimnogammarus vitatus", no = taxon),
         taxon = ifelse(test = taxon == "Periphyton_NA",
                        yes = "Periphyton", no = taxon)) %>%
ggplot(aes(x = log10(ppcp_sum), y = ((C18.3w3 + C18.4w3) / (C20.5w3 + C22.5w3)))) +
  geom_point(size = 3) +
  facet_wrap(~ taxon) +
  geom_smooth(method = "lm") +
  xlab(label = "log10([Total PPCP])") +
  ylab(label = expression(frac(18:3~omega~3 + 18:4~omega~3, 20:5~omega~3 + 22:6~omega~3))) +
  theme_bw() +
  theme(text = element_text(size = 20))

ggsave(filename = "ppcp_fa_plot.png", plot = ppcp_fa_plot, device = "png",
       path = "../figures/", width = 12, height = 6, units = "in")

# Create linear models to explore which sewage indicators relate with fatty acid profiles
# Our analysis considers C18 fatty acids in comparison to C20,22 fatty acids.
# These two fatty acid groups roughly reflect Green:Diatom algal signature,
# which should inrease with an increasing sewage signal.

# First compare essential fatty acid ratios in periphyton with total PPCP concentration
peri_ppcp_lm <- lm(formula = ((C18.3w3 + C18.4w3) / (C20.5w3 + C22.5w3)) ~ log10(ppcp_sum),
               data = filter(fatty_acid_prop_ppcp_meta_dist, Genus == "Periphyton"))

summary(peri_ppcp_lm)

# Second compare essential fatty acid ratios in invertebrates with total PPCP conentration
invert_ppcp_lm <- lm(((C18.3w3 + C18.4w3) / (C20.5w3 + C22.5w3 + C22.6w3)) ~ log10(ppcp_sum),
                   data = filter(fatty_acid_prop_ppcp_meta_dist, Genus != "Periphyton"))

summary(invert_ppcp_lm)
