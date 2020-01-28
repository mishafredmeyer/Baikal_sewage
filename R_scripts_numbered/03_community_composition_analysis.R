library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(viridisLite)
library(vegan)
library(factoextra)
library(cluster)

# 1. Load the data --------------------------------------------------------

invertebrates <- read.csv("../cleaned_data/invertebrates_20190320.csv", header = TRUE)
periphyton <- read.csv("../cleaned_data/periphyton_20190320.csv", header = TRUE)

metadata <- read.csv("../cleaned_data/metadata_20190320.csv", header = TRUE)
distance <- read.csv("../cleaned_data/distance_20190320.csv", header = TRUE)

metadata_dist <- full_join(metadata, distance)

ppcp <- read.csv("../cleaned_data/ppcp_20190320.csv", header = TRUE)

ppcp_meta_dist <- full_join(ppcp, metadata_dist) %>%
  filter(!(Site %in% c("OS-1", "OS-2", "OS-3")))

# 2. Periphyton Analysis ---------------------------------------------------

low <- c("BGO-1", "BGO-2", "KD-1", "KD-2")
mod <- c("BGO-3", "BK-2", "BK-3", "MS-1")
high <- c("BK-1", "SM-1", "EM-1", "LI-3", "LI-2")

periphyton_meta_dist <- full_join(periphyton, ppcp_meta_dist) %>%
  mutate(POPULATION_INTENSITY = (SOUTH_SHORE_DIST*(POP_SOUTH_DEV/SOUTH_SHORE_AREA))/SOUTH_DEV_DIST)

##Univariate Analysis
periphyton_meta_dist_long <- periphyton_meta_dist %>%
  gather(Taxon, Count, desmidales:ulothrix) %>%
  filter(!(Taxon %in% c("desmidales", "pediastrum"))) %>%
  group_by(Site) %>%
  mutate(Total_count = sum(Count))

periphyton_meta_dist_long$Site <- factor(periphyton_meta_dist_long$Site, 
                                    levels = c("BGO-3", "BGO-1", "BGO-2", "KD-1", "KD-2",
                                               "MS-1", "BK-3", "BK-2", "BK-1", "SM-1", "EM-1",
                                               "LI-3", "LI-2", "LI-1"))

periphyton_meta_dist_plot <- ggplot(periphyton_meta_dist_long[!is.na(periphyton_meta_dist_long$Site), ]) +
  geom_bar(aes(Site, Total_count), fill = 'grey80', stat = "identity") +
  geom_bar(aes(Site, Count, fill = Taxon), stat = "identity") +
  scale_fill_viridis(discrete = TRUE, option = "inferno") +
  facet_wrap(~ Taxon) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ylab("Number of cells") +
  xlab("Site (Arranged by increasing population intensity)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x=element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(size=20),
        strip.text.x = element_text(size=27),
        strip.background = element_rect(fill = "white"),
        panel.background = element_rect(color = "black"),
        axis.title = element_text(size = 30),
        #axis.text.x = element_text("none"),
        axis.text.y = element_text(size = 24),
        axis.title.y=element_text(margin=margin(0,20,0,0)), 
        axis.title.x=element_text(size = 24, margin=margin(20,0,0,0)),
        legend.text=element_text(size=16), 
        #axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave("periphyton_univariate.png", periphyton_meta_dist_plot, "png",
       path = "../figures/", width = 11, height = 8.5, units = "in")


##Multivariate Analysis

periphyton_meta_dist_wide <- periphyton_meta_dist %>%
  filter(!(Site %in% c("OS-1", "OS-2", "OS-3"))) %>% 
  mutate(PI_group = ifelse(Site %in% high, "High", "NULL"),
         PI_group = ifelse((Site %in% mod) | (Site %in% low), "Mod/Low", PI_group)) %>%
  select(-desmidales, -pediastrum) %>%
  as.data.frame()

periphyton_nmds <- metaMDS(periphyton_meta_dist_wide[periphyton_meta_dist_wide$Site != "LI-1", 2:5],  try = 100)
periphyton_nmds

data_scores <- as.data.frame(scores(periphyton_nmds))
data_scores$Site <- periphyton_meta_dist_wide[periphyton_meta_dist_wide$Site != "LI-1",1]
data_scores <- inner_join(data_scores, ppcp_meta_dist) %>%
  mutate(PI_group = ifelse(Site %in% high, "High", "NULL"),
         PI_group = ifelse(Site %in% mod, "Mod", PI_group),
         PI_group = ifelse(Site %in% low, "Low", PI_group), 
         POPULATION_INTENSITY = (SOUTH_SHORE_DIST*(POP_SOUTH_DEV/SOUTH_SHORE_AREA))/SOUTH_DEV_DIST)

data_scores$PI_group <- factor(data_scores$PI_group,
                               levels = c("High", "Mod", "Low"))

species_scores <- as.data.frame(scores(periphyton_nmds, "species")) 
species_scores$species <- rownames(species_scores)


periphyton_PI_group_plot <- ggplot() + 
  #geom_text(data=species_scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5, size = 10) +  # add the species labels
  geom_point(data=data_scores[data_scores$Site != "LI-1", ],aes(x=NMDS1,y=NMDS2, size = log10(POPULATION_INTENSITY + 1), color = PI_group)) +  # add the site labels
  scale_size_continuous(range = c(5,20), guide = FALSE) +
  #ggtitle("NMDS with Periphyton Proportions") +
  xlim(c(-0.5, 0.5)) +
  scale_color_manual(values = inferno(15)[c(3,8,11)], name = "Distance-weighted Population Grouping") +
  guides(colour = guide_legend(override.aes = list(size=10))) +
  #coord_equal() +
  annotate("label", x = 0, y = -0.35, size = 10, 
           label = paste("Stress: ", round(periphyton_nmds$stress, digits = 3))) +
  theme_minimal() +
  theme(legend.position = "right",
        text = element_text(size = 24))
periphyton_PI_group_plot

ggsave("periphyton_PI_group_plot.png", periphyton_PI_group_plot, device = "png", 
       path = "../figures/",  height = 10, width = 20, dpi = 300)

per_cluster <- fviz_nbclust(periphyton_meta_dist_wide[periphyton_meta_dist_wide$Site != "LI-1", 2:5], kmeans, method = "wss")

per_cluster

adonis(periphyton_meta_dist_wide[, 2:5] 
       ~ periphyton_meta_dist_wide[, 28], 
       data = periphyton_meta_dist_wide[periphyton_meta_dist_wide$Site != "LI-1", ],
       method = "bray", permutations = 999)



# 3. Invertebrate Analysis ---------------------------------------------------

invertebrate_meta_dist <- full_join(invertebrates, ppcp_meta_dist) %>%
  mutate(POPULATION_INTENSITY = (SOUTH_SHORE_DIST*(POP_SOUTH_DEV/SOUTH_SHORE_AREA))/SOUTH_DEV_DIST)

amphipods <- c("Eulimnogammarus", "Poekilogammarus", "Pallasea", "Hyallela", "Cryptoropus", "Brandtia")
molluscs <- c("Acroloxidae", "Baicaliidae",  "Benedictidate", "Planorbidae", "Valvatidae")

panel.cor <- function(x, y, digits = 2, cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  # correlation coefficient
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste("r= ", txt, sep = "")
  text(0.5, 0.6, txt)
  
  # p-value calculation
  p <- cor.test(x, y)$p.value
  txt2 <- format(c(p, 0.123456789), digits = digits)[1]
  txt2 <- paste("p= ", txt2, sep = "")
  if(p<0.01) txt2 <- paste("p= ", "<0.01", sep = "")
  text(0.5, 0.4, txt2)
}

invertebrates_long <- invertebrate_meta_dist %>%
  select(Site:Valvatidae) %>%
  gather(Taxon, Count, Acroloxidae:Valvatidae) %>%
  filter(!grepl("juvenile", Taxon)) %>%
  separate(Taxon, c("Genus", "Species")) %>%
  group_by(Site, Genus) %>%
  summarize(Total_Genus = sum(Count)) %>%
  ungroup() %>%
  group_by(Site) %>%
  mutate(Total_Site = sum(Total_Genus)) %>%
  unique() %>%
  dplyr::mutate(GROUPING = ifelse(Genus %in% amphipods, "Amphipod", NA),
                GROUPING = ifelse(Genus %in% molluscs, "Mollusc", GROUPING),
                GROUPING = ifelse(Genus == "Asellidae", "Isopod", GROUPING),
                GROUPING = ifelse(Genus == "caddisflies", "Caddisflies", GROUPING),
                GROUPING = ifelse(Genus == "flatworms", "Planaria", GROUPING),
                GROUPING = ifelse(Genus == "Leeches", "Hirudinea", GROUPING)) %>%
  filter(!is.na(GROUPING))

invertebrates_long$Site <- factor(invertebrates_long$Site, levels = c("LI-1", "LI-3", "LI-2", "EM-1", "SM-1", "BK-2", "BK-1", "BK-3", "MS-1", "KD-2", "KD-1", "BGO-2", "BGO-1", "BGO-3"))

ggplot(invertebrates_long) +
  geom_bar(aes(x = Site, y = Total_Site), alpha =0.5, stat = "identity") +
  geom_bar(aes(x= Site, y = Total_Genus, fill = Genus), stat = "identity") +
  facet_wrap(~ Genus) +
  #scale_y_continuous(limits = c(0,1300), expand = c(0, 0)) +
  #scale_fill_manual(values = get_palette("Paired", length(unique(invertebrates_long$Genus)))) +
  xlab("Locations (left-to-right reads South to North)") +
  ylab("Number of Individuals") +
  theme_classic() +
  theme(legend.position = "bottom",
        axis.text.x=element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(size=20),
        strip.text.x = element_text(size=20),
        strip.background = element_rect(fill = "white"),
        panel.background = element_rect(color = "black"),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y=element_text(margin=margin(0,20,0,0)), 
        axis.title.x=element_text(margin=margin(20,0,0,0)),
        legend.text=element_text(size=16))

invertebrates_condensed_wide <- invertebrate_meta_dist %>%
  select(Site:Valvatidae) %>%
  gather(Taxon, Count, Acroloxidae:Valvatidae) %>%
  filter(!grepl("juvenile", Taxon)) %>%
  separate(Taxon, c("Genus", "Species")) %>%
  group_by(Site, Genus) %>%
  summarize(Total_Genus = sum(Count)) %>%
  ungroup() %>%
  group_by(Site) %>%
  mutate(Total_Site = sum(Total_Genus)) %>%
  unique() %>%
  spread(Genus, Total_Genus)
  
head(invertebrates_condensed_wide)

pairs(invertebrates_condensed_wide[,3:18], upper.panel=panel.cor)

correlated <- c("Acroloxidae", "Asellidae", "Baicaliidae", "Benedictidate",
                "Brandtia", "Hyallela", "Maackia", "Choronomids")

invertebrates_without_corr_long <- invertebrate_meta_dist %>%
  select(Site:Valvatidae) %>%
  gather(Taxon, Count, Acroloxidae:Valvatidae) %>%
  filter(!grepl("juvenile", Taxon)) %>%
  separate(Taxon, c("Genus", "Species")) %>%
  filter(!(Genus %in% correlated)) %>%
  group_by(Site, Genus) %>%
  summarize(Total_Genus = sum(Count)) %>% 
  ungroup() %>%
  group_by(Site) %>%
  mutate(Total_Site = sum(Total_Genus)) %>%
  unique() 

invertebrates_without_corr_long$Genus <- factor(invertebrates_without_corr_long$Genus, 
                                                   levels = c("Cryptoropus", "Eulimnogammarus", 
                                                              "Pallasea", "Poekilogammarus", 
                                                              "Planorbidae", "Valvatidae", 
                                                              "Caddisflies", "Flatworms", "Leeches"))
invertebrates_without_corr_long$Site <- factor(invertebrates_without_corr_long$Site, 
                                                levels = c("BGO-3", "BGO-1", "BGO-2", "KD-1", "KD-2",
                                                           "MS-1", "BK-3", "BK-2", "BK-1", "SM-1", "EM-1",
                                                           "LI-3", "LI-2", "LI-1"))

invertebrate_wo_corr_plot <- invertebrates_without_corr_long %>%
  mutate(Group = ifelse(Genus %in% amphipods, "Amphipod", Genus),
         Group = ifelse(Genus %in% molluscs, "Mollusc", Group),
         Genus = factor(Genus, levels = c("Cryptoropus", "Eulimnogammarus", "Pallasea", "Poekilogammarus", 
                                          "Planorbidae", "Valvatidae", "Caddisflies", "Flatworms", 
                                          "Leeches"))) %>%
  ggplot() +
  geom_bar(aes(x = Site, y = Total_Site), alpha =0.5, stat = "identity") +
  geom_bar(aes(x= Site, y = Total_Genus, fill = as.factor(Group)), stat = "identity") +
  scale_fill_viridis(discrete = TRUE, option = "inferno") +
  facet_wrap(~ Genus) +
  #scale_y_continuous(limits = c(0,1300), expand = c(0, 0)) +
  xlab("Site (Arranged by increasing population intensity)") +
  ylab("Number of Individuals") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x=element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(size=20),
        strip.text.x = element_text(size=20),
        strip.background = element_rect(fill = "white"),
        panel.background = element_rect(color = "black"),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y=element_text(margin=margin(0,20,0,0)), 
        axis.title.x=element_text(margin=margin(20,0,0,0)),
        legend.text=element_text(size=16))
ggsave("../figures/invertebrate_wo_corr_plot.png", invertebrate_wo_corr_plot, device = "png",
       width = 18, height = 12, dpi = 300)

invertebrates_without_corr_wide <- invertebrate_meta_dist %>%
  select(Site:Valvatidae) %>%
  gather(Taxon, Count, Acroloxidae:Valvatidae) %>%
  filter(!grepl("juvenile", Taxon)) %>%
  separate(Taxon, c("Genus", "Species")) %>%
  filter(!(Genus %in% correlated)) %>%
  group_by(Site, Genus) %>%
  summarize(Total_Genus = sum(Count)) %>%
  ungroup() %>%
  group_by(Site) %>%
  mutate(Total_Site = sum(Total_Genus)) %>%
  unique() %>%
  spread(Genus, Total_Genus)


invert_cluster <- fviz_nbclust(invertebrates_without_corr_wide[ , 3:11], kmeans, method = "wss")
invert_cluster

invertebrates_metaMDS <- metaMDS(invertebrates_without_corr_wide[ , 3:11], distance = "bray", try = 100)
invertebrates_metaMDS

data_scores <- as.data.frame(scores(invertebrates_metaMDS))
data_scores$Site <- invertebrates_without_corr_wide$Site
data_scores <- full_join(data_scores, ppcp_meta_dist) %>%
  mutate(POPULATION_INTENSITY = (SOUTH_SHORE_DIST*(POP_SOUTH_DEV/SOUTH_SHORE_AREA))/SOUTH_DEV_DIST, 
         POP_GROUP = ifelse(Site %in% low, "Low", NA),
         POP_GROUP = ifelse(Site %in% mod, "Mod", POP_GROUP),
         POP_GROUP = ifelse(Site %in% high, "High", POP_GROUP))
data_scores$POP_GROUP <- factor(data_scores$POP_GROUP, 
                                levels = c("High", "Mod", "Low"))

species_scores <- as.data.frame(scores(invertebrates_metaMDS, "species")) 
species_scores$species <- rownames(species_scores)

inverts_without_corr_nmds <- ggplot() + 
  #geom_polygon(data=data.scores,aes(x=NMDS1,y=NMDS2,fill=SITE.GROUP,group=SITE.GROUP),alpha=0.50) +
  #geom_text(data=species_scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5, size = 10, color = "grey20") +  # add the species labels
  #geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2),size=3) + # add the point markers
  geom_point(data= drop_na(data_scores),aes(x=NMDS1,y=NMDS2, size = log10(PPCP.SUM), 
                                 color = POP_GROUP)) +  # add the site labels
  scale_size_continuous(range = c(8,20), guide = FALSE) +
  scale_color_manual(values = inferno(15)[c(3,8,11, 14)], 
                     name = "Distance-weighted Population Grouping") + 
  guides(colour = guide_legend(override.aes = list(size=10))) +
  coord_equal() +
  ylim(c(-.5, .5)) +
  xlim(c(-.5, .5)) +
  #geom_hline(yintercept = -0.1, linetype = "dashed", color = "red") +
  #theme_grey() +
  annotate("label", x = -0.15, y = -0.4, size = 10, 
           label = paste("Stress: ", round(invertebrates_metaMDS$stress, digits = 3))) +
  theme(legend.position = "right",
        strip.text.x = element_text(size=20, color = "grey80"),
        text = element_text(size = 24), 
        axis.title.y=element_text(margin=margin(0,20,0,0)), 
        axis.title.x=element_text(margin=margin(20,0,0,0)),
        panel.background = element_rect("white"),
        panel.grid.major = element_line(colour = "grey80"), 
        panel.grid.minor = element_line(colour = "grey80"),
        axis.ticks = element_line(color = "grey80"))
inverts_without_corr_nmds

ggsave("../figures/inverts_without_corr_nmds.png", inverts_without_corr_nmds, device = "png",
       height = 10, width = 20, dpi = 300)

invertebrates_without_corr_meta_dist_wide <- full_join(invertebrates_without_corr_wide, ppcp_meta_dist) %>%
  mutate(POPULATION_INTENSITY = (SOUTH_SHORE_DIST*(POP_SOUTH_DEV/SOUTH_SHORE_AREA))/SOUTH_DEV_DIST, 
     POP_GROUP = ifelse(Site %in% c(low, mod), "LOW/MOD", NA),
     POP_GROUP = ifelse(Site %in% mod, "MOD", POP_GROUP),
     POP_GROUP = ifelse(Site %in% high, "HIGH", POP_GROUP)) %>%
  drop_na() %>%
  data.frame()

adonis(invertebrates_without_corr_meta_dist_wide[, 3:11] ~ 
       invertebrates_without_corr_meta_dist_wide[,34],
       data = invertebrates_without_corr_meta_dist_wide[invertebrates_without_corr_meta_dist_wide$Site != "LI-1",],
       method = "bray")

