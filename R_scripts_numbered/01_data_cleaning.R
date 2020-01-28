library(dplyr)
library(reshape2)
library(lubridate)
library(stringr)

# Clean PPCP data --------

ppcp_orig <- read.csv("../original_data/PPCP_Baikal_orig_20180524.csv", header = TRUE)

ppcp <- ppcp_orig %>%
  select(Sample_ID, Caffeine, Acetaminophen, X1.7.Dimethylxanthine, Cotinine) %>%
  filter(Sample_ID == 'LI-1'|
           Sample_ID == 'LI-2'|
           Sample_ID == 'LI-3'|
           Sample_ID == 'BK-1' |
           Sample_ID == 'BK-2'| 
           Sample_ID == 'BK-3'|
           Sample_ID == 'BGO-1'| 
           Sample_ID == 'BGO-2'|
           Sample_ID == 'BGO-3'| 
           Sample_ID == 'KD-1'|
           Sample_ID == 'KD-2'|
           Sample_ID == 'EM-1'|
           Sample_ID == 'MS-1'|
           Sample_ID == 'SM-1'|
           Sample_ID == 'OS-1'|
           Sample_ID == 'OS-2'|
           Sample_ID == 'OS-3') %>%
  group_by(Sample_ID) %>%
  mutate(PPCP.SUM = Caffeine+Acetaminophen+X1.7.Dimethylxanthine+Cotinine) %>%
  rename(Paraxanthine = X1.7.Dimethylxanthine, 
         Site = Sample_ID)
head(ppcp)
write.csv(ppcp, "../cleaned_data/ppcp_20190320.csv", row.names = FALSE)

# Nutrient cleaning -----------------------------------------

nutrients_orig <- read.csv("../original_data/baikal_nearshore_nutrient_data_20151009.csv", header = TRUE)

nutrients <- nutrients_orig %>%
  group_by(sample) %>%
  summarize(mean_NH4_mg_dm3 = mean(NH4_mg_dm3), 
            mean_NO3_mg_dm3 = mean(NO3_mg_dm3),
            mean_TP_mg_dm3 = mean(TP_mg_dm3),
            mean_TPO43_mg_dm3 = mean(TPO43_mg_dm3)) %>%
  rename(Site = sample)
head(nutrients)
write.csv(nutrients, "../cleaned_data/nutrients_20190320.csv", row.names = FALSE)

# Clean Chlorophyll a data ------------------------------------------------

chla_orig <- read.csv("../original_data/chlorophyll_20170117.csv", header = TRUE)

chlorophylla <- chla_orig %>%
  select(Station, chl_conc) %>%
  filter(Station == 'LI-1'|
           Station == 'LI-2'|
           Station == 'LI-3'|
           Station == 'BK-1' |
           Station == 'BK-2'| 
           Station == 'BK-3'|
           Station == 'BGO-1'| 
           Station == 'BGO-2'|
           Station == 'BGO-3'| 
           Station == 'KD-1'|
           Station == 'KD-2'|
           Station == 'EM-1'|
           Station == 'MS-1'|
           Station == 'SM-1'|
           Station == 'OS-1'|
           Station == 'OS-2'|
           Station == 'OS-3') %>%
  group_by(Station) %>%
  summarize(mean_chlorophylla = mean(chl_conc)) %>%
  rename(Site = Station)

head(chlorophylla)
write.csv(chlorophylla, "../cleaned_data/chlorophylla_20190320.csv", row.names = FALSE)

# Metadata cleaning for Lat & Long ----------------------------------------
metadata_orig <- read.csv("../original_data/baikal_nearshore_metadata_201508.csv", header = TRUE)

metadata <- metadata_orig %>%
  select(loc_site, lat, long, depth, dist_to_shore, air_temp, surface_temp, mid_temp, bottom_temp) %>%
  rename(Site = loc_site)

head(metadata)
write.csv(metadata, "../cleaned_data/metadata_20190320.csv", row.names = FALSE)

# Macroinbertebrate cleaning ---------------------------------------------------------

inverts_orig <- read.csv("../original_data/macroinvert_community_QAQC_mfm_20171108.csv", header = TRUE)

inverts <- inverts_orig %>%
  select(-X) %>%
  gather(Site, Count, MS1.3:BK1.3) %>%
  rename(Taxon = Invertebrate) %>%
  mutate(Site = gsub(".", "_", Site, fixed = TRUE),
         Taxon = gsub(" ", "_", Taxon),
         Count = ifelse(is.na(Count), 0, Count)) %>%
  separate(Site, c("Location", "Replicate", "Duplicate"), remove = FALSE) %>%
  filter(Taxon != "Propapaidae") %>%
  group_by(Location, Replicate, Duplicate, Taxon) %>%
  summarize(sum_Count = sum(Count)) %>%
  ungroup() %>%
  group_by(Location, Taxon) %>%
  summarize(mean_Count = mean(sum_Count)) %>%
  separate(Taxon, c("Genus", "Species", "Subspecies")) %>%
  mutate(Genus = ifelse(Genus == "E", "Eulimnogammarus", Genus),
         Genus = ifelse(Genus == "Eulimno", "Eulimnogammarus", Genus),
         Genus = ifelse(Genus == "flatworms_", "Flatworms", Genus),
         Genus = ifelse(Genus == "caddisflies", "Caddisflies", Genus),
         Genus = ifelse(Genus == "pallasea", "Pallasea", Genus),
         Genus = ifelse(Genus == "hyallela", "Hyallela", Genus),
         Genus = ifelse(Genus == "poekilo", "Poekilogammarus", Genus),
         Genus = ifelse(Genus == "Poekilo", "Poekilogammarus", Genus),
         Genus = ifelse(Genus == "valvatidae", "Valvatidae", Genus), 
         Species = ifelse(Species == "spp", NA, Species)) %>%
  unite(Taxon, Genus, Species, Subspecies) %>%
  mutate(Taxon = gsub("_NA_NA", "", Taxon),
         Taxon = gsub("_NA", "", Taxon),
         Taxon = ifelse(Taxon == "flatworms_", "Flatworms", Taxon), 
         Taxon = ifelse(Taxon == "Hyallela_cziarnianski_", "Hyallela_cziarnianski", Taxon),
         Taxon = ifelse(Taxon == "Pallasea_cancellus_", "Pallasea_cancellus", Taxon),
         Taxon = ifelse(Taxon == "choronomids_", "Choronomids", Taxon)) %>%
  spread(Taxon, mean_Count) %>%
  select(-Total) %>%
  separate(Location, c("Location", "Number"), sep = -1) %>%
  unite("Site", c("Location", "Number"), sep = "-")
head(inverts)
write.csv(inverts, "../cleaned_data/invertebrates_20190320.csv", row.names = FALSE)


# Periphyton cleaning -----------------------------------------------------

periphyton_orig <- read.csv("../original_data/periphyton_20180917.csv", header = TRUE)

periphyton <- periphyton_orig %>%
  select(-date, -rep, -contains("filament"), -Lyngbya) %>%
  filter(!is.na(diatom)) %>%
  gather(TAXON, COUNT, diatom:desmidales) %>%
  mutate(TAXON = ifelse(TAXON == "tetraporales", "tetrasporales", TAXON)) %>%
  group_by(site, TAXON) %>%
  mutate(MEAN = mean(COUNT)) %>%
  ungroup() %>%
  group_by(site) %>%
  select(-counts, -comments, -COUNT) %>%
  unique() %>%
  spread(TAXON, MEAN) %>%
  rename(Site = site) %>%
  separate(Site, c("Location", "Number"), sep = -1) %>%
  unite("Site", c("Location", "Number"), sep = "-") %>%
  gather(Taxon, Count, desmidales:ulothrix) %>% 
  group_by(Site) %>%
  mutate(Total_count = sum(Count),
         Prop = Count/Total_count) %>%
  select(-Total_count, -Count) %>%
  spread(Taxon, Prop) %>%
  filter(!(Site %in% c("OS-1", "OS-2", "OS-3"))) %>% 
  mutate(PI_group = ifelse(Site %in% high, "High", "NULL"),
         PI_group = ifelse((Site %in% mod) | (Site %in% low), "Mod/Low", PI_group)) %>%
  as.data.frame()
head(periphyton)

write.csv(periphyton, "../cleaned_data/periphyton_20190320.csv", row.names = FALSE)


# Stable Isotope Cleaning -------------------------------------------------

stable_isotopes_orig <- read.csv("../original_data/sia_results_mfm_20170509.csv", header = TRUE)

stable_isotopes <- stable_isotopes_orig %>%
  separate(Identifier, c("Site", "Genus", "Species"), sep = " ") %>%
  mutate(Genus = ifelse(Genus == "E.", "Eulimnogammarus", Genus), 
         Genus = ifelse(Genus == "P.", "Pallasea", Genus), 
         Species = ifelse(Species == "can", "cancellus", Species), 
         Species = ifelse(Species == "ver", "verrucosus", Species),
         Species = ifelse(Species == "vitetus", "vitatus", Species),
         Species = ifelse(Species == "veruossus", "verrucosus", Species),
         Species = ifelse(Species == "cyan", "cyaneus", Species),
         Species = ifelse(Species == "Sp.", "Splash", Species))
head(stable_isotopes)

write.csv(stable_isotopes, "../cleaned_data/stable_isotopes_20190320.csv", row.names = FALSE)


# Fatty Acid cleaning -----------------------------------------------------------

fatty_acid_orig <- read.csv("../original_data/BaikalFAs_wt_20180322.csv", header = TRUE)

fatty_acid <- fatty_acid_orig %>%
  select(-GC_ID, -sample.) %>%
  separate(spp, c("Genus", "Species")) %>%
  mutate(Genus = ifelse(Genus == "E" & Species == "ver", "Eulimnogammarus", Genus), 
         Genus = ifelse(Genus == "E" & Species == "vitatus", "Eulimnogammarus", Genus),
         Genus = ifelse(Genus == "E" & Species == "cyan", "Eulimnogammarus", Genus),
         Genus = ifelse(Genus == "P" & Species == "can", "Pallasea", Genus),
         Genus = ifelse(Genus == "Spl" & Species == "zone", "Splash", Genus),
         Species = ifelse(Genus == "Eulimnogammarus" & Species == "ver", "verrucosus", Species),
         Species = ifelse(Genus == "Eulimnogammarus" & Species == "cyan", "cyaneus", Species),
         Species = ifelse(Genus == "Pallasea" & Species == "can", "cancellus", Species),
         Species = ifelse(Genus == "Spl" & Species == "zone", "Zone", Species)) %>%
  rename(Site = location)
head(fatty_acid)

write.csv(fatty_acid, "../cleaned_data/fatty_acid_20190320.csv", row.names = FALSE)


# Total Lipid Cleaning -------------------------------------------------------------

total_lipid_orig <- read.csv("../original_data/Baikal.total.lipid.mfm.20180322.csv", header = TRUE)

total_lipid <- total_lipid_orig %>%
  select(-sample.num) %>%
  separate(sample.id, c("Site", "SPP"), sep = "\\,") %>%
  separate("SPP", c("Genus", "Species"), sep = "[.]") %>%
  mutate(Genus = ifelse(Genus == " E" , "Eulimnogammarus", Genus), 
         Genus = ifelse(Genus == "E" , "Eulimnogammarus", Genus), 
         Genus = ifelse(Genus == "E" & Species == "ver", "Eulimnogammarus", Genus),
         Genus = ifelse(Genus == "E" & Species == "vitatus", "Eulimnogammarus", Genus),
         Genus = ifelse(Genus == "E" & Species == "cyan", "Eulimnogammarus", Genus),
         Genus = ifelse(Genus == " P" & Species == "can", "Pallasea", Genus),
         Genus = ifelse(Genus == " Spl" & Species == "zone", "Splash", Genus),
         Genus = ifelse(grepl("Drapa", Genus), "Drapa", Genus),
         Genus = ifelse(grepl("Hyalella", Genus), "Hyalella", Genus),
         Genus = ifelse(grepl("Snails", Genus), "Snails", Genus),
         Species = ifelse(Genus == "Eulimnogammarus" & Species == "ver", "verrucosus", Species),
         Species = ifelse(Genus == "Eulimnogammarus" & Species == "ever", "verrucosus", Species),
         Species = ifelse(Genus == "Eulimnogammarus" & grepl("verucossus", Species), "verrucosus", Species),
         Species = ifelse(Genus == "Eulimnogammarus" & grepl("vitatus", Species), "vitatus", Species),
         Species = ifelse(Genus == "Eulimnogammarus" & Species == "cyan", "cyaneus", Species),
         Species = ifelse(Genus == "Eulimnogammarus" & grepl("cyan", Species), "cyaneus", Species),
         Species = ifelse(Genus == "Pallasea" & Species == "can", "cancellus", Species),
         Species = ifelse(Genus == "Spl" & Species == "zone", "Zone", Species)) %>%
  rename(total_lipid_mg_per_g = total.lipid.mg.g)
head(total_lipid)

write.csv(total_lipid, "../cleaned_data/total_lipid_20190320.csv", row.names = FALSE)


# Microplastics Cleaning -----------------------------------------------------------

microplastics_orig <- read.csv("../original_data/microplastics_mfm_20171010.csv", header = TRUE)

microplastics <- microplastics_orig %>%
  select(-comments) %>%
  mutate(VOLUME_FILTERED = (volume * volume_rep)/1000, 
         TOTAL_MP = fragments + fibers + beads,
         DENSITY = TOTAL_MP/VOLUME_FILTERED,
         FRAG_DENSITY = fragments/VOLUME_FILTERED,
         FIBER_DENSITY = fibers/VOLUME_FILTERED,
         BEAD_DENSITY = beads/VOLUME_FILTERED) %>%
  unite(Site, location, site, sep = "-") %>%
  select(Site, rep, TOTAL_MP, VOLUME_FILTERED, DENSITY, FRAG_DENSITY, 
         FIBER_DENSITY, BEAD_DENSITY) %>%
  filter(rep != "C") %>%
  group_by(Site) %>%
  summarize(mean_volume_filtered = mean(VOLUME_FILTERED), 
            mean_total_microplastics = mean(TOTAL_MP),
            mean_microplastic_density = mean(DENSITY),
            mean_fragment_density = mean(FRAG_DENSITY),
            mean_fiber_density = mean(FIBER_DENSITY),
            mean_bead_density = mean(BEAD_DENSITY))
head(microplastics)

write.csv(microplastics, "../cleaned_data/microplastics_20190320.csv", row.names = FALSE)

# Calculated Distance Metrics Cleaning ------------------------------------

distance_orig <- read.csv("../original_data/baikal_site_distances_mfm_20180517.csv", header = TRUE)

distance <- distance_orig %>%
  select(-SOUTH_UP, -NORTH_UP, -X, -X.1, -X.2) %>%
  rename(Site = Sample_ID)

head(distance)

write.csv(distance, "../cleaned_data/distance_20190320.csv", row.names = FALSE)