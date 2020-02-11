# 1. Load packages --------------------------------------------------------

library(tidyverse)
library(reshape2)
library(lubridate)
library(stringr)


# 2. Load and clean PPCP data ---------------------------------------------

ppcp_orig <- read.csv(file = "../original_data/PPCP_Baikal_orig_20180524.csv",
                      header = TRUE)

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
  mutate(PPCP.SUM = Caffeine + Acetaminophen + X1.7.Dimethylxanthine + Cotinine) %>%
  rename(Paraxanthine = X1.7.Dimethylxanthine, Site = Sample_ID)

head(ppcp)

# Export new version of the data
write.csv(x = ppcp, file = "../cleaned_data/ppcp_20190320.csv",
          row.names = FALSE)


# 2. Load and clean nutrient data -----------------------------------------

nutrients_orig <- read.csv(file = "../original_data/baikal_nearshore_nutrient_data_20151009.csv",
                           header = TRUE)

nutrients <- nutrients_orig %>%
  group_by(sample) %>%
  summarize(mean_NH4_mg_dm3 = mean(NH4_mg_dm3),
            mean_NO3_mg_dm3 = mean(NO3_mg_dm3),
            mean_TP_mg_dm3 = mean(TP_mg_dm3),
            mean_TPO43_mg_dm3 = mean(TPO43_mg_dm3)) %>%
  rename(Site = sample)

head(nutrients)

write.csv(x = nutrients, file = "../cleaned_data/nutrients_20190320.csv",
          row.names = FALSE)


# 3. Load and clean chlorophyll a data ------------------------------------

chla_orig <- read.csv(file = "../original_data/chlorophyll_20170117.csv",
                      header = TRUE)

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

write.csv(x = chlorophylla, file = "../cleaned_data/chlorophylla_20190320.csv",
          row.names = FALSE)


# 4. Load and clean lat/long metadata -------------------------------------

metadata_orig <- read.csv(file = "../original_data/baikal_nearshore_metadata_201508.csv",
                          header = TRUE)

metadata <- metadata_orig %>%
  select(loc_site, lat, long, depth, dist_to_shore, air_temp, surface_temp,
         mid_temp, bottom_temp) %>%
  rename(Site = loc_site)

head(metadata)

write.csv(x = metadata, file = "../cleaned_data/metadata_20190320.csv",
          row.names = FALSE)


# 5. Load and clean macroinbertebrate data --------------------------------

inverts_orig <- read.csv(file = "../original_data/macroinvert_community_QAQC_mfm_20171108.csv",
                         header = TRUE)

inverts <- inverts_orig %>%
  select(-X) %>%
  gather(key = Site, value = Count, MS1.3:BK1.3) %>%
  rename(Taxon = Invertebrate) %>%
  mutate(Site = gsub(pattern = ".", replacement = "_", x = Site, fixed = TRUE),
         Taxon = gsub(pattern = " ", replacement = "_", x = Taxon),
         Count = ifelse(test = is.na(Count), yes = 0, no = Count)) %>%
  separate(col = Site, into = c("Location", "Replicate", "Duplicate"),
           remove = FALSE) %>%
  filter(Taxon != "Propapaidae") %>%
  group_by(Location, Replicate, Duplicate, Taxon) %>%
  summarize(sum_Count = sum(Count)) %>%
  ungroup() %>%
  group_by(Location, Taxon) %>%
  summarize(mean_Count = mean(sum_Count)) %>%
  separate(col = Taxon, into = c("Genus", "Species", "Subspecies")) %>%
  mutate(Genus = ifelse(test = Genus == "E",
                        yes = "Eulimnogammarus", no = Genus),
         Genus = ifelse(test = Genus == "Eulimno",
                        yes = "Eulimnogammarus", no = Genus),
         Genus = ifelse(test = Genus == "flatworms_",
                        yes = "Flatworms", no = Genus),
         Genus = ifelse(test = Genus == "caddisflies",
                        yes = "Caddisflies", no = Genus),
         Genus = ifelse(test = Genus == "pallasea",
                        yes = "Pallasea", no = Genus),
         Genus = ifelse(test = Genus == "hyallela",
                        yes = "Hyallela", no = Genus),
         Genus = ifelse(test = Genus == "poekilo",
                        yes = "Poekilogammarus", no = Genus),
         Genus = ifelse(test = Genus == "Poekilo",
                        yes = "Poekilogammarus", no = Genus),
         Genus = ifelse(test = Genus == "valvatidae",
                        yes = "Valvatidae", no = Genus), 
         Species = ifelse(test = Species == "spp",
                          yes = NA, no = Species)) %>%
  unite(col = "Taxon", Genus, Species, Subspecies) %>%
  mutate(Taxon = gsub(pattern = "_NA_NA", replacement = "", x = Taxon),
         Taxon = gsub(pattern = "_NA", replacement = "", x = Taxon),
         Taxon = ifelse(test = Taxon == "flatworms_",
                        yes = "Flatworms", no = Taxon), 
         Taxon = ifelse(test = Taxon == "Hyallela_cziarnianski_",
                        yes = "Hyallela_cziarnianski", no = Taxon),
         Taxon = ifelse(test = Taxon == "Pallasea_cancellus_",
                        yes = "Pallasea_cancellus", no = Taxon),
         Taxon = ifelse(test = Taxon == "choronomids_",
                        yes = "Choronomids", no = Taxon)) %>%
  spread(key = Taxon, value = mean_Count) %>%
  select(-Total) %>%
  separate(col = Location, into = c("Location", "Number"), sep = -1) %>%
  unite(col = "Site", Location, Number, sep = "-")

head(inverts)

write.csv(x = inverts, file = "../cleaned_data/invertebrates_20190320.csv",
          row.names = FALSE)


# 6. Load and clean periphyton data ---------------------------------------

periphyton_orig <- read.csv(file = "../original_data/periphyton_20180917.csv",
                            header = TRUE)

periphyton <- periphyton_orig %>%
  select(-date, -rep, -contains("filament"), -Lyngbya) %>%
  filter(!is.na(diatom)) %>%
  gather(key = TAXON, value = COUNT, diatom:desmidales) %>%
  mutate(TAXON = ifelse(test = TAXON == "tetraporales",
                        yes = "tetrasporales", no = TAXON)) %>%
  group_by(site, TAXON) %>%
  mutate(MEAN = mean(COUNT)) %>%
  ungroup() %>%
  group_by(site) %>%
  select(-counts, -comments, -COUNT) %>%
  unique() %>%
  spread(key = TAXON, value = MEAN) %>%
  rename(Site = site) %>%
  separate(col = Site, into = c("Location", "Number"), sep = -1) %>%
  unite(col = "Site", c("Location", "Number"), sep = "-") %>%
  gather(key = Taxon, value = Count, desmidales:ulothrix) %>% 
  group_by(Site) %>%
  mutate(Total_count = sum(Count),
         Prop = Count/Total_count) %>%
  select(-Total_count, -Count) %>%
  spread(key = Taxon, value = Prop) %>%
  filter(!(Site %in% c("OS-1", "OS-2", "OS-3"))) %>% 
  mutate(PI_group = ifelse(test = Site %in% high,
                           yes = "High", no = "NULL"),
         PI_group = ifelse(test = (Site %in% mod) | (Site %in% low),
                           yes = "Mod/Low", no = PI_group)) %>%
  as.data.frame()

head(periphyton)

write.csv(x = periphyton, file = "../cleaned_data/periphyton_20190320.csv",
          row.names = FALSE)


# 7. Load and clean stable isotope data -----------------------------------

stable_isotopes_orig <- read.csv(file = "../original_data/sia_results_mfm_20170509.csv",
                                 header = TRUE)

stable_isotopes <- stable_isotopes_orig %>%
  separate(col = Identifier, into = c("Site", "Genus", "Species"), sep = " ") %>%
  mutate(Genus = ifelse(test = Genus == "E.",
                        yes = "Eulimnogammarus", Genus), 
         Genus = ifelse(test = Genus == "P.",
                        yes = "Pallasea", no = Genus), 
         Species = ifelse(test = Species == "can",
                          yes = "cancellus", no = Species), 
         Species = ifelse(test = Species == "ver",
                          yes = "verrucosus", no = Species),
         Species = ifelse(test = Species == "vitetus",
                          yes = "vitatus", no = Species),
         Species = ifelse(test = Species == "veruossus",
                          yes = "verrucosus", no = Species),
         Species = ifelse(test = Species == "cyan",
                          yes = "cyaneus", no = Species),
         Species = ifelse(test = Species == "Sp.",
                          yes = "Splash", no = Species))

head(stable_isotopes)

write.csv(x = stable_isotopes, file = "../cleaned_data/stable_isotopes_20190320.csv",
          row.names = FALSE)


# 8. Load and clean fatty acid data ---------------------------------------

fatty_acid_orig <- read.csv(file = "../original_data/BaikalFAs_wt_20180322.csv",
                            header = TRUE)

fatty_acid <- fatty_acid_orig %>%
  select(-GC_ID, -sample.) %>%
  separate(col = spp, into = c("Genus", "Species")) %>%
  mutate(Genus = ifelse(test = Genus == "E" & Species == "ver",
                        yes = "Eulimnogammarus", no = Genus), 
         Genus = ifelse(test = Genus == "E" & Species == "vitatus",
                        yes = "Eulimnogammarus", no = Genus),
         Genus = ifelse(test = Genus == "E" & Species == "cyan",
                        yes = "Eulimnogammarus", no = Genus),
         Genus = ifelse(test = Genus == "P" & Species == "can",
                        yes = "Pallasea", no = Genus),
         Genus = ifelse(test = Genus == "Spl" & Species == "zone",
                        yes = "Splash", no = Genus),
         Species = ifelse(test = Genus == "Eulimnogammarus" & Species == "ver",
                          yes = "verrucosus", no = Species),
         Species = ifelse(test = Genus == "Eulimnogammarus" & Species == "cyan",
                          yes = "cyaneus", no = Species),
         Species = ifelse(test = Genus == "Pallasea" & Species == "can",
                          yes = "cancellus", no = Species),
         Species = ifelse(test = Genus == "Spl" & Species == "zone",
                          yes = "Zone", no = Species)) %>%
  rename(Site = location)

head(fatty_acid)

write.csv(x = fatty_acid, file = "../cleaned_data/fatty_acid_20190320.csv",
          row.names = FALSE)


# 9. Load and clean total lipid data --------------------------------------

total_lipid_orig <- read.csv(file = "../original_data/Baikal.total.lipid.mfm.20180322.csv",
                             header = TRUE)

total_lipid <- total_lipid_orig %>%
  select(-sample.num) %>%
  separate(col = sample.id, into = c("Site", "SPP"), sep = "\\,") %>%
  separate(col = SPP, into = c("Genus", "Species"), sep = "[.]") %>%
  mutate(Genus = ifelse(test = Genus == " E" ,
                        yes = "Eulimnogammarus", no = Genus), 
         Genus = ifelse(test = Genus == "E" ,
                        yes = "Eulimnogammarus", no = Genus), 
         Genus = ifelse(test = Genus == "E" & Species == "ver",
                        yes = "Eulimnogammarus", no = Genus),
         Genus = ifelse(test = Genus == "E" & Species == "vitatus",
                        yes = "Eulimnogammarus", no = Genus),
         Genus = ifelse(test = Genus == "E" & Species == "cyan",
                        yes = "Eulimnogammarus", no = Genus),
         Genus = ifelse(test = Genus == " P" & Species == "can",
                        yes = "Pallasea", no = Genus),
         Genus = ifelse(test = Genus == " Spl" & Species == "zone",
                        yes = "Splash", no = Genus),
         Genus = ifelse(test = grepl("Drapa", Genus),
                        yes = "Drapa", no = Genus),
         Genus = ifelse(test = grepl("Hyalella", Genus),
                        yes = "Hyalella", no = Genus),
         Genus = ifelse(test = grepl("Snails", Genus),
                        yes = "Snails", no = Genus),
         Species = ifelse(test = Genus == "Eulimnogammarus" & Species == "ver",
                          yes = "verrucosus", no = Species),
         Species = ifelse(test = Genus == "Eulimnogammarus" & Species == "ever",
                          yes = "verrucosus", no = Species),
         Species = ifelse(test = Genus == "Eulimnogammarus" & grepl(pattern = "verucossus", x = Species),
                          yes = "verrucosus", no = Species),
         Species = ifelse(test = Genus == "Eulimnogammarus" & grepl(pattern = "vitatus", x = Species),
                          yes = "vitatus", no = Species),
         Species = ifelse(test = Genus == "Eulimnogammarus" & Species == "cyan",
                          yes = "cyaneus", no = Species),
         Species = ifelse(test = Genus == "Eulimnogammarus" & grepl(pattern = "cyan", x = Species),
                          yes = "cyaneus", no = Species),
         Species = ifelse(test = Genus == "Pallasea" & Species == "can",
                          yes = "cancellus", no = Species),
         Species = ifelse(test = Genus == "Spl" & Species == "zone",
                          yes = "Zone", no = Species)) %>%
  rename(total_lipid_mg_per_g = total.lipid.mg.g)
head(total_lipid)

write.csv(x = total_lipid, file = "../cleaned_data/total_lipid_20190320.csv",
          row.names = FALSE)


# 10. Load and clean microplastics data -----------------------------------

microplastics_orig <- read.csv(file = "../original_data/microplastics_mfm_20171010.csv",
                               header = TRUE)

microplastics <- microplastics_orig %>%
  select(-comments) %>%
  mutate(VOLUME_FILTERED = (volume * volume_rep)/1000, 
         TOTAL_MP = fragments + fibers + beads,
         DENSITY = TOTAL_MP/VOLUME_FILTERED,
         FRAG_DENSITY = fragments/VOLUME_FILTERED,
         FIBER_DENSITY = fibers/VOLUME_FILTERED,
         BEAD_DENSITY = beads/VOLUME_FILTERED) %>%
  unite(col = Site, location, site, sep = "-") %>%
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

write.csv(x = microplastics, file = "../cleaned_data/microplastics_20190320.csv",
          row.names = FALSE)


# 11. Load and clean calculated distance metrics --------------------------

distance_orig <- read.csv(file = "../original_data/baikal_site_distances_mfm_20180517.csv",
                          header = TRUE)

distance <- distance_orig %>%
  select(-SOUTH_UP, -NORTH_UP, -X, -X.1, -X.2) %>%
  rename(Site = Sample_ID)

head(distance)

write.csv(x = distance, file = "../cleaned_data/distance_20190320.csv",
          row.names = FALSE)