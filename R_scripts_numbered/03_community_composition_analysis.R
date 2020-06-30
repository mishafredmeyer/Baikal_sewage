# This script is meant to perform multivariate and some univariate analysis of
# macroinvertebrates and periphyton community composition along Lake Baikal's
# southwestern shoreline. These community compositions are then related to
# sewage indicator concentrations and inverse distance weight population
# along the shoreline.

# 1. Load packages --------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(viridisLite)
library(vegan)
library(factoextra)
library(cluster)
library(ggpubr)

# Pull in function for plotting correlations
source("panel_cor_function.R")


# 2. Load the data --------------------------------------------------------

# Invertebrate data
invertebrates <- read.csv("../cleaned_data/invertebrates.csv",
                          header = TRUE, stringsAsFactors = FALSE)

# Periphyton data
periphyton <- read.csv("../cleaned_data/periphyton.csv",
                       header = TRUE, stringsAsFactors = FALSE)

# Site metadata
metadata <- read.csv("../cleaned_data/metadata.csv",
                     header = TRUE, stringsAsFactors = FALSE)

# Site distance data
distance <- read.csv("../cleaned_data/distance_weighted_population_metrics.csv",
                     header = TRUE, stringsAsFactors = FALSE)

# Join site metadata with distance data
metadata_dist <- full_join(x = metadata, y = distance, by = "Site")

# PPCP data
ppcp <- read.csv("../cleaned_data/ppcp.csv",
                 header = TRUE, stringsAsFactors = FALSE)

# Join PPCP data with metadata/distance and select sites of interest
ppcp_meta_dist <- full_join(x = ppcp, y = metadata_dist, by = "Site") %>%
  filter(!(Site %in% c("OS-1", "OS-2", "OS-3")))


# Define low/moderate/high sites
# These groupings are derived from the inverse distance weighted population metric,
# where sites with "low" IDW population are recorded in the "low" vector, moderate in "mod",
# and high in the "high" vector. These same delineations are used in both the periphyton
# and invertebrate analysis.

low <- c("BGO-1", "BGO-2", "BGO-3", "KD-1", "KD-2", "MS-1")
mod <- c("BK-2", "BK-3", "SM-1")
high <- c("BK-1", "EM-1", "LI-3", "LI-2", "LI-1")


# 3. Periphyton analysis ---------------------------------------------------

# Join periphyton data with metadata/distance and create custom metric
periphyton_meta_dist <- full_join(x = periphyton, y = ppcp_meta_dist,
                                  by = "Site")


# 3.1 Univariate analysis -------------------------------------------------

# Convert periphyton data to long format and add total count
periphyton_meta_dist_long <- periphyton_meta_dist %>%
  gather(key = Taxon, value = Count, desmidales:ulothrix) %>%
  filter(!(Taxon %in% c("desmidales", "pediastrum"))) %>%
  group_by(Site) %>%
  mutate(Total_count = sum(Count))

# Rework Site column as a factor
periphyton_meta_dist_long$Site <- factor(x = periphyton_meta_dist_long$Site,
                                         levels = c("BGO-2",  "BGO-3", "KD-1",
                                                    "KD-2", "BGO-1", "MS-1",
                                                    "BK-3", "BK-2", "SM-1",
                                                    "BK-1", "EM-1", "LI-1",
                                                    "LI-2", "LI-3"))

# Remove data with Site = NA
periphyton_long_clean <- periphyton_meta_dist_long %>% filter(!is.na(Site))

# Plot periphyton counts as a function of site and taxa
periphyton_meta_dist_plot <- ggplot(data = periphyton_long_clean) +
  geom_bar(aes(x = Site, y = Total_count), fill = "grey80", stat = "identity") +
  geom_bar(aes(x = Site, y = Count, fill = Taxon), stat = "identity") +
  scale_fill_viridis(discrete = TRUE, option = "inferno") +
  facet_wrap(~ Taxon) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ylab("Number of cells") +
  xlab("Site (Arranged by increasing distance-weighted population)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(size = 20),
        strip.text.x = element_text(size = 27),
        strip.background = element_rect(fill = "white"),
        panel.background = element_rect(color = "black"),
        axis.title = element_text(size = 30),
        axis.text.y = element_text(size = 24),
        axis.title.y = element_text(margin = margin(0, 20, 0, 0)),
        axis.title.x = element_text(size = 24, margin = margin(20, 0, 0, 0)),
        legend.text = element_text(size = 16),
        axis.ticks.x = element_blank())

# This plot is not included in the associated manuscript but is included to
# offer insight for multivariate analyeses.

ggsave(filename = "periphyton_univariate.png", plot = periphyton_meta_dist_plot,
       device = "png", path = "../figures/", width = 11, height = 8.5,
       units = "in")


# 3.2 Multivariate analysis -----------------------------------------------

# Clean dataset and add IDW_pop_group
periphyton_meta_dist_wide <- periphyton_meta_dist %>%
  filter(!(Site %in% c("OS-1", "OS-2", "OS-3"))) %>%
  mutate(IDW_pop_group = ifelse(test = Site %in% high,
                                yes = "High", no = "NULL"),
         IDW_pop_group = ifelse(test = (Site %in% mod),
                                yes = "Mod", no = IDW_pop_group),
         IDW_pop_group = ifelse(test = (Site %in% low),
                                yes = "Low", no = IDW_pop_group)) %>%
  as.data.frame()

# Define community for NMDS
peri_community <- periphyton_meta_dist_wide %>%
  select(diatom:ulothrix)

# Run NMDS
periphyton_nmds <- metaMDS(comm = peri_community,  try = 100)
periphyton_nmds

# Pull scores from NMDS and add site data
data_scores <- as.data.frame(scores(x = periphyton_nmds))
data_scores$Site <- periphyton_meta_dist_wide %>%
  pull(Site)

# Join scores with PPCP data and code into PI groups
data_scores <- inner_join(x = data_scores, y = ppcp_meta_dist, by = "Site") %>%
  mutate(IDW_pop_group = ifelse(test = Site %in% high,
                                yes = "High", no = "NULL"),
         IDW_pop_group = ifelse(test = Site %in% mod,
                                yes = "Mod", no = IDW_pop_group),
         IDW_pop_group = ifelse(test = Site %in% low,
                                yes = "Low", no = IDW_pop_group))

# Rework IDW_pop_group column as a factor
data_scores$IDW_pop_group <- factor(x = data_scores$IDW_pop_group,
                                    levels = c("High", "Mod", "Low"))

# Pull species scores from NMDS
species_scores <- as.data.frame(scores(x = periphyton_nmds, display = "species"))
species_scores$species <- rownames(species_scores)

# Plot the NMDS
periphyton_IDW_pop_group_plot <- ggplot() +
  geom_point(data = data_scores,
             aes(x = NMDS1, y = NMDS2, size = log10(distance_weighted_population),
                 color = IDW_pop_group)) +
  scale_size_continuous(range = c(5, 20), guide = FALSE) +
  scale_color_manual(values = inferno(15)[c(3, 8, 11)],
                     name = "IDW Population Grouping") +
  guides(colour = guide_legend(override.aes = list(size = 10))) +
  annotate("label", x = 0, y = -0.35, size = 10,
           label = paste("Stress: ", round(periphyton_nmds$stress, digits = 3))) +
  theme_minimal() +
  theme(legend.position = "right",
        text = element_text(size = 24))

periphyton_IDW_pop_group_plot

# Export plot
ggsave(filename = "../figures/periphyton_IDW_pop_group_plot.png",
       plot = periphyton_IDW_pop_group_plot, device = "png", height = 10,
       width = 20, dpi = 300)

# Vizualize optimum cluster number
peri_cluster <- fviz_nbclust(x = peri_community, FUNcluster = kmeans,
                             method = "wss")

peri_cluster

# Export plot
# This plot correcspond with figure S1a in the associated manuscript.
ggsave(filename = "../figures/periphyton_kmeans_analysis.png",
       plot = peri_cluster, device = "png", height = 10, width = 20,
       dpi = 300)

# Run PERMANOVA
adonis(formula = periphyton_meta_dist_wide[, 2:7]
       ~ periphyton_meta_dist_wide[, 22],
       data = periphyton_meta_dist_wide,
       method = "bray", permutations = 999)

simper_results <- simper(comm = periphyton_meta_dist_wide[, 2:7],
                         group = periphyton_meta_dist_wide[, 22],
                         permutations = 999)

summary(simper_results)


# 4. Invertebrate analysis ---------------------------------------------------


# 4.1 Univariate analysis -------------------------------------------------

# Convert invertebrate data to long format and add total count
invertebrate_meta_dist <- full_join(x = invertebrates, y = ppcp_meta_dist,
                                    by = "Site")

# Define amphipod genera
amphipods <- c("Eulimnogammarus", "Poekilogammarus", "Pallasea", "Hyallela",
               "Cryptoropus", "Brandtia")
# Define mollusc genera
molluscs <- c("Acroloxidae", "Baicaliidae",  "Benedictidate", "Planorbidae",
              "Valvatidae")

# Make invertebrate data long format, sum counts, and group genera
invertebrates_long <- invertebrate_meta_dist %>%
  select(Site:Valvatidae) %>%
  gather(key = Taxon, value = Count, Acroloxidae:Valvatidae) %>%
  filter(!grepl("juvenile", Taxon)) %>%
  separate(col = Taxon, into = c("Genus", "Species")) %>%
  group_by(Site, Genus) %>%
  summarize(Total_Genus = sum(Count)) %>%
  ungroup() %>%
  group_by(Site) %>%
  mutate(Total_Site = sum(Total_Genus)) %>%
  unique() %>%
  mutate(GROUPING = ifelse(test = Genus %in% amphipods,
                           yes = "Amphipod", no = NA),
         GROUPING = ifelse(test = Genus %in% molluscs,
                           yes = "Mollusc", no = GROUPING),
         GROUPING = ifelse(test = Genus == "Asellidae",
                           yes = "Isopod", no = GROUPING),
         GROUPING = ifelse(test = Genus == "caddisflies",
                           yes = "Caddisflies", no = GROUPING),
         GROUPING = ifelse(test = Genus == "flatworms",
                           yes = "Planaria", no = GROUPING),
         GROUPING = ifelse(test = Genus == "Leeches",
                           yes = "Hirudinea", no = GROUPING)) %>%
  filter(!is.na(GROUPING))

# Rework Site column as a factor
invertebrates_long$Site <- factor(x = invertebrates_long$Site,
                                  levels = c("BGO-2",  "BGO-3", "KD-1",
                                             "KD-2", "BGO-1", "MS-1",
                                             "BK-3", "BK-2", "SM-1",
                                             "BK-1", "EM-1", "LI-1",
                                             "LI-2", "LI-3"))

# Plot invert counts by site and genus
ggplot(invertebrates_long) +
  geom_bar(aes(x = Site, y = Total_Site), alpha = 0.5, stat = "identity") +
  geom_bar(aes(x = Site, y = Total_Genus, fill = Genus), stat = "identity") +
  facet_wrap(~ Genus) +
  xlab("Locations (left-to-right reads South to North)") +
  ylab("Number of Individuals") +
  theme_classic() +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        strip.background = element_rect(fill = "white"),
        panel.background = element_rect(color = "black"),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(margin = margin(0, 20, 0, 0)),
        axis.title.x = element_text(margin = margin(20, 0, 0, 0)),
        legend.text = element_text(size = 16))

# Create a cleaned up wide format invertebrate genus dataset
invertebrates_condensed_wide <- invertebrate_meta_dist %>%
  select(Site:Valvatidae) %>%
  gather(key = Taxon, value = Count, Acroloxidae:Valvatidae) %>%
  filter(!grepl("juvenile", Taxon)) %>%
  separate(col = Taxon, into = c("Genus", "Species")) %>%
  group_by(Site, Genus) %>%
  summarize(Total_Genus = sum(Count)) %>%
  ungroup() %>%
  group_by(Site) %>%
  mutate(Total_Site = sum(Total_Genus)) %>%
  unique() %>%
  spread(key = Genus, value = Total_Genus)

# Check
head(invertebrates_condensed_wide)

# Check for correlations
pairs(invertebrates_condensed_wide[, 3:18], upper.panel = panel.cor)

# Choronomids and Hyallela were not preserved well, and therefore were not well
# identified in the dataset. This step removes them from further analysis.

poorly_preserved <- c("Choronomids", "Hyallela")

# Remove poorly preserved genera from dataset
invertebrates_well_preserved_long <- invertebrate_meta_dist %>%
  select(Site:Valvatidae) %>%
  gather(key = Taxon, value = Count, Acroloxidae:Valvatidae) %>%
  filter(!grepl("juvenile", Taxon)) %>%
  separate(col = Taxon, into = c("Genus", "Species")) %>%
  filter(!(Genus %in% poorly_preserved)) %>%
  group_by(Site, Genus) %>%
  summarize(Total_Genus = sum(Count)) %>%
  ungroup() %>%
  group_by(Site) %>%
  mutate(Total_Site = sum(Total_Genus)) %>%
  unique() %>%
  ungroup()

# Rework Genus and Site columns as factors
invertebrates_well_preserved_long_2 <- invertebrates_well_preserved_long %>%
  mutate(Genus = factor(x = Genus,
                        levels = c("Brandtia", "Cryptoropus", "Eulimnogammarus",
                                   "Pallasea", "Poekilogammarus",
                                   "Acroloxidae", "Baicaliidae", "Benedictidate",
                                   "Maackia", "Planorbidae", "Valvatidae",
                                   "Asellidae", "Caddisflies", "Flatworms",
                                   "Leeches")),
         Site = factor(x = Site,
                       levels = c("BGO-2",  "BGO-3", "KD-1",
                                  "KD-2", "BGO-1", "MS-1",
                                  "BK-3", "BK-2", "SM-1",
                                  "BK-1", "EM-1", "LI-1",
                                  "LI-2", "LI-3")))

# Prep invert data for plotting and then do so
invertebrate_well_preserved_plot <- invertebrates_well_preserved_long %>%
  mutate(Group = ifelse(test = Genus %in% amphipods,
                        yes = "Amphipod", no = Genus),
         Group = ifelse(test = Genus %in% molluscs,
                        yes = "Mollusc", no = Group),
         Genus = factor(x = Genus,
                        levels = c("Brandtia", "Cryptoropus", "Eulimnogammarus",
                                   "Pallasea", "Poekilogammarus",
                                   "Acroloxidae", "Baicaliidae", "Benedictidate",
                                   "Maackia", "Planorbidae", "Valvatidae",
                                   "Asellidae", "Caddisflies", "Flatworms",
                                   "Leeches"))) %>%
  ggplot() +
  geom_bar(aes(x = Site, y = Total_Site), alpha = 0.5, stat = "identity") +
  geom_bar(aes(x = Site, y = Total_Genus, fill = as.factor(Group)),
           stat = "identity") +
  scale_fill_viridis(discrete = TRUE, option = "inferno") +
  facet_wrap(~ Genus) +
  xlab("Site (Arranged by increasing population intensity)") +
  ylab("Number of Individuals") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        strip.background = element_rect(fill = "white"),
        panel.background = element_rect(color = "black"),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(margin = margin(0, 20, 0, 0)),
        axis.title.x = element_text(margin = margin(20, 0, 0, 0)),
        legend.text = element_text(size = 16))

# Export plot
ggsave(filename = "../figures/invertebrate_well_preserved_plot.png",
       plot = invertebrate_well_preserved_plot, device = "png",
       width = 18, height = 12, dpi = 300)


# 4.2 Multivariate analysis -----------------------------------------------

# Create a cleaned up wide format invertebrate genus dataset with well preserved
# genera. In this step we also remove rare species (i.e., taxa representing less
# than 1% of the intercommunity abundance).
invertebrates_well_preserved_wide <- invertebrate_meta_dist %>%
  select(Site:Valvatidae) %>%
  gather(key = Taxon, value = Count, Acroloxidae:Valvatidae) %>%
  filter(!grepl("juvenile", Taxon)) %>%
  separate(col = Taxon, into = c("Genus", "Species")) %>%
  filter(!(Genus %in% poorly_preserved)) %>%
  group_by(Site, Genus) %>%
  summarize(Total_Genus = sum(Count)) %>%
  ungroup() %>%
  group_by(Site) %>%
  mutate(Total_Site = sum(Total_Genus),
         prop_genus = Total_Genus / Total_Site) %>%
  ungroup() %>%
  group_by(Genus) %>%
  mutate(mean_prop = mean(prop_genus, na.rm = TRUE)) %>%
  filter(mean_prop >= 0.01) %>%
  select(-mean_prop, -prop_genus) %>%
  spread(key = Genus, value = Total_Genus)

# Square-root transform
inver_comm <- sqrt(invertebrates_well_preserved_wide[, 3:12])

# Vizualize optimum cluster number for invert community
invert_cluster <- fviz_nbclust(x = inver_comm,
                               FUNcluster = kmeans, method = "wss")
invert_cluster

# Export plot
ggsave(filename = "../figures/invertebrate_kmeans_analysis.png",
       plot = invert_cluster, device = "png",
       width = 18, height = 12, dpi = 300)

# Run invert NMDS
invertebrates_metaMDS <- metaMDS(comm = inver_comm, distance = "bray", try = 100)
invertebrates_metaMDS

# Pull scores from NMDS and add site data
data_scores <- as.data.frame(scores(x = invertebrates_metaMDS))
data_scores$Site <- invertebrates_well_preserved_wide$Site

# Join scores with PPCP data and code into IDW_pop_group
data_scores <- full_join(x = data_scores, y = ppcp_meta_dist, by = "Site") %>%
  mutate(IDW_pop_group = ifelse(test = Site %in% c(low, mod),
                                yes = "Low/Mod", no = NA),
         IDW_pop_group = ifelse(test = Site %in% high,
                                yes = "High", no = IDW_pop_group),
         IDW_pop_group = factor(x = IDW_pop_group,
                                levels = c("High", "Low/Mod")))

# Plot NMDS
inverts_well_preserved_nmds <- ggplot() +
  geom_point(data = data_scores,
             aes(x = NMDS1, y = NMDS2, size = log10(distance_weighted_population),
                 color = IDW_pop_group)) +
  scale_size_continuous(range = c(8, 20), guide = FALSE) +
  scale_color_manual(values = inferno(15)[c(3, 8, 11, 14)],
                     name = "IDW Population Grouping") +
  guides(colour = guide_legend(override.aes = list(size = 10))) +
  coord_equal() +
  ylim(c(-.5, .5)) +
  xlim(c(-.5, .5)) +
  annotate("label", x = -0.15, y = 0.4, size = 10,
           label = paste("Stress: ",
                         round(invertebrates_metaMDS$stress, digits = 3))) +
  theme(legend.position = "right",
        strip.text.x = element_text(size = 20, color = "grey80"),
        text = element_text(size = 24),
        axis.title.y = element_text(margin = margin(0, 20, 0, 0)),
        axis.title.x = element_text(margin = margin(20, 0, 0, 0)),
        panel.background = element_rect("white"),
        panel.grid.major = element_line(colour = "grey80"),
        panel.grid.minor = element_line(colour = "grey80"),
        axis.ticks = element_line(color = "grey80"))

inverts_well_preserved_nmds

# Export plot. This plot is Figure 4 in the associated manuscript.
ggsave(filename = "../figures/inverts_well_preserved_nmds.png",
       plot = inverts_well_preserved_nmds, device = "png",
       height = 10, width = 20, dpi = 300)

# Re-join the metadata and distance data back in with invert data
inverts_well_preserved_meta_dist_wide <- full_join(x = invertebrates_well_preserved_wide,
                                                   y = ppcp_meta_dist,
                                                   by = "Site") %>%
  mutate(IDW_pop_group = ifelse(test = Site %in% c(low, mod),
                                yes = "Low/Mod", no = NA),
         IDW_pop_group = ifelse(test = Site %in% high,
                                yes = "High", no = IDW_pop_group),
         distance_weighted_population = log10(distance_weighted_population)) %>%
  data.frame()

# Run PERMANOVA
adonis(formula = sqrt(inverts_well_preserved_meta_dist_wide[, 3:12]) ~
         inverts_well_preserved_meta_dist_wide[, 27],
       data = inverts_well_preserved_meta_dist_wide,
       method = "bray")

# Run SIMPER
simper_results <- simper(comm = sqrt(inverts_well_preserved_meta_dist_wide[, 3:12]),
                         group = inverts_well_preserved_meta_dist_wide[, 27],
                         permutations = 999)

summary(simper_results)

# Combine k-means WSS plots into one.
# This plot is associated with Figure S1 in the associated manuscript.

ggarrange(peri_cluster, invert_cluster, ncol = 2, nrow = 1, labels = "AUTO") %>%
  ggexport(filename = "../figures/keams_combined_plot.png",
           hieght = 600, width = 1200, res = 120)

ggarrange(periphyton_IDW_pop_group_plot, inverts_well_preserved_nmds,
          ncol = 2, nrow = 1, labels = "AUTO") %>%
  ggexport(filename = "../figures/nmds_combined_plot.png",
           hieght = 2000, width = 3000, res = 120)
