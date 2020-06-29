# 1. Load packages --------------------------------------------------------

library(tidyverse)
library(reshape2)
library(lubridate)
library(stringr)


# 2. Load and clean PPCP data ---------------------------------------------

# Raw data
ppcp_orig <- read.csv(file = "../original_data/PPCP_Baikal_orig_20180524.csv",
                      header = TRUE)

# Select variables and sites of interest, sum PPCPs
ppcp <- ppcp_orig %>%
  select(Sample_ID, Caffeine, Acetaminophen, X1.7.Dimethylxanthine, Cotinine) %>%
  filter(Sample_ID == "LI-1" |
           Sample_ID == "LI-2" |
           Sample_ID == "LI-3" |
           Sample_ID == "BK-1" |
           Sample_ID == "BK-2" |
           Sample_ID == "BK-3" |
           Sample_ID == "BGO-1" |
           Sample_ID == "BGO-2" |
           Sample_ID == "BGO-3" |
           Sample_ID == "KD-1" |
           Sample_ID == "KD-2" |
           Sample_ID == "EM-1" |
           Sample_ID == "MS-1" |
           Sample_ID == "SM-1" |
           Sample_ID == "OS-1" |
           Sample_ID == "OS-2" |
           Sample_ID == "OS-3") %>%
  group_by(Sample_ID) %>%
  mutate(ppcp_sum = Caffeine + Acetaminophen + X1.7.Dimethylxanthine + Cotinine) %>%
  rename(Paraxanthine = X1.7.Dimethylxanthine, Site = Sample_ID)

# Take a look
head(ppcp)

# Export new version of the data
write.csv(x = ppcp, file = "../cleaned_data/ppcp.csv",
          row.names = FALSE)


# 2. Load and clean nutrient data -----------------------------------------

nutrients_orig <- read.csv(file = "../original_data/baikal_nearshore_nutrient_data_20151009.csv",
                           header = TRUE)

# Take nutrients averages by site
nutrients <- nutrients_orig %>%
  group_by(sample) %>%
  summarize(mean_NH4_mg_dm3 = mean(NH4_mg_dm3),
            mean_NO3_mg_dm3 = mean(NO3_mg_dm3),
            mean_TP_mg_dm3 = mean(TP_mg_dm3),
            mean_TPO43_mg_dm3 = mean(TPO43_mg_dm3)) %>%
  rename(Site = sample)

head(nutrients)

write.csv(x = nutrients, file = "../cleaned_data/nutrients.csv",
          row.names = FALSE)


# 3. Load and clean chlorophyll a data ------------------------------------

chla_orig <- read.csv(file = "../original_data/chlorophyll_20170117.csv",
                      header = TRUE)

# Select sites of interest and take average chl a by site
chlorophylla <- chla_orig %>%
  select(Station, chl_conc) %>%
  filter(Station == "LI-1" |
           Station == "LI-2" |
           Station == "LI-3" |
           Station == "BK-1" |
           Station == "BK-2" |
           Station == "BK-3" |
           Station == "BGO-1" |
           Station == "BGO-2" |
           Station == "BGO-3" |
           Station == "KD-1" |
           Station == "KD-2" |
           Station == "EM-1" |
           Station == "MS-1" |
           Station == "SM-1" |
           Station == "OS-1" |
           Station == "OS-2" |
           Station == "OS-3") %>%
  group_by(Station) %>%
  summarize(mean_chlorophylla = mean(chl_conc)) %>%
  rename(Site = Station)

head(chlorophylla)

write.csv(x = chlorophylla, file = "../cleaned_data/chlorophylla.csv",
          row.names = FALSE)


# 4. Load and clean lat/long metadata -------------------------------------

metadata_orig <- read.csv(file = "../original_data/baikal_nearshore_metadata_201508.csv",
                          header = TRUE)

# Select columns of interest
metadata <- metadata_orig %>%
  select(loc_site, lat, long, depth, dist_to_shore, air_temp, surface_temp,
         mid_temp, bottom_temp) %>%
  rename(Site = loc_site)

head(metadata)

write.csv(x = metadata, file = "../cleaned_data/metadata.csv",
          row.names = FALSE)


# 5. Load and clean macroinbertebrate data --------------------------------

inverts_orig <- read.csv(file = "../original_data/macroinvert_community_QAQC_mfm_20171108.csv",
                         header = TRUE)

# Clean up naming in data and add up counts by taxon
inverts_summarized <- inverts_orig %>%
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
  ungroup()

# Take mean counts by taxon, flesh out taxonomic info, spread to wide format
inverts_wide <- inverts_summarized %>%
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

head(inverts_wide)

write.csv(x = inverts_wide, file = "../cleaned_data/invertebrates.csv",
          row.names = FALSE)


# 6. Load and clean periphyton data ---------------------------------------

periphyton_orig <- read.csv(file = "../original_data/periphyton_20180917.csv",
                            header = TRUE)

# Make long format, take mean counts by taxon
periphyton_summarized <- periphyton_orig %>%
  select(-date, -rep, -contains("filament"), -Lyngbya) %>%
  filter(!is.na(diatom)) %>%
  gather(key = TAXON, value = COUNT, diatom:desmidales) %>%
  mutate(TAXON = ifelse(test = TAXON == "tetraporales",
                        yes = "tetrasporales", no = TAXON)) %>%
  group_by(Site = site, TAXON) %>%
  summarize(MEAN = mean(COUNT)) %>%
  ungroup()

# Make wide version of periphyton using proportion data
periphyton_wide <- periphyton_summarized %>%
  spread(key = TAXON, value = MEAN) %>%
  separate(col = Site, into = c("Location", "Number"), sep = -1) %>%
  unite(col = "Site", Location, Number, sep = "-") %>%
  # gather(key = Taxon, value = Count, desmidales:ulothrix) %>%
  # group_by(Site) %>%
  # mutate(Total_count = sum(Count),
  #        Prop = Count / Total_count) %>%
  # select(-Total_count, -Count) %>%
  # spread(key = Taxon, value = Prop) %>%
  filter(!(Site %in% c("OS-1", "OS-2", "OS-3"))) %>%
  as.data.frame()

head(periphyton_wide)

write.csv(x = periphyton_wide, file = "../cleaned_data/periphyton.csv",
          row.names = FALSE)


# 7. Load and clean stable isotope data -----------------------------------

stable_isotopes_orig <- read.csv(file = "../original_data/sia_results_mfm_20170509.csv",
                                 header = TRUE)

# Parse identifier data into site and taxonomic data
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

write.csv(x = stable_isotopes, file = "../cleaned_data/stable_isotopes.csv",
          row.names = FALSE)


# 8. Load and clean fatty acid data ---------------------------------------

fatty_acid_orig <- read.csv(file = "../original_data/BaikalFAs_wt_20180322.csv",
                            header = TRUE)

# Parse spp column into taxonomic data
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

write.csv(x = fatty_acid, file = "../cleaned_data/fatty_acid.csv",
          row.names = FALSE)


# 9. Load and clean total lipid data --------------------------------------

total_lipid_orig <- read.csv(file = "../original_data/Baikal.total.lipid.mfm.20180322.csv",
                             header = TRUE)

# Parse sample.id column into site and taxonomic data
total_lipid <- total_lipid_orig %>%
  select(-sample.num) %>%
  separate(col = sample.id, into = c("Site", "SPP"), sep = "\\,") %>%
  separate(col = SPP, into = c("Genus", "Species"), sep = "[.]") %>%
  mutate(Genus = ifelse(test = Genus == " E",
                        yes = "Eulimnogammarus", no = Genus),
         Genus = ifelse(test = Genus == "E",
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
         Genus = ifelse(test = grepl(pattern = "Drapa", x = Genus),
                        yes = "Drapa", no = Genus),
         Genus = ifelse(test = grepl(pattern = "Hyalella", x = Genus),
                        yes = "Hyalella", no = Genus),
         Genus = ifelse(test = grepl(pattern = "Snails", x = Genus),
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

write.csv(x = total_lipid, file = "../cleaned_data/total_lipid.csv",
          row.names = FALSE)


# 10. Load and clean microplastics data -----------------------------------

microplastics_orig <- read.csv(file = "../original_data/microplastics_mfm_20171010.csv",
                               header = TRUE)

# Run microplastics post-processing calcs, then average by site
microplastics_uncorrected <- microplastics_orig %>%
  select(-comments) %>%
  unite(col = "Site", location, site, sep = "-") %>%
  filter(rep != "C") 

# Separate out the controls so that they can be removed from
# the experimental counts
microplastics_controls <- microplastics_orig %>%
  select(-comments) %>%
  filter(rep == "C") %>%
  unite(col = "Site", location, site, sep = "-") %>%
  group_by(Site) %>%
  summarize(fiber_controls = mean(fibers),
            fragment_controls = mean(fragments),
            beads_controls = mean(beads))

microplastics_corrected <- left_join(x = microplastics_uncorrected, microplastics_controls,
                                     by = "Site") %>%
  mutate(VOLUME_FILTERED = (volume_rep*volume) / 1000,
         fragments_corrected = fragments - fragment_controls,
         fragments_corrected = ifelse(test = fragments_corrected < 0, 
                                      yes = 0, no = fragments_corrected),
         fibers_corrected = fibers - fiber_controls,
         fibers_corrected = ifelse(test = fibers_corrected < 0, 
                                   yes = 0, no = fibers_corrected),
         beads_corrected = beads - beads_controls,
         beads_corrected = ifelse(test = beads_corrected < 0, 
                                  yes = 0, no = beads_corrected),
         total_microplastics = fragments_corrected + fibers_corrected + beads_corrected, 
         density = total_microplastics / VOLUME_FILTERED,
         fragment_density = fragments_corrected / VOLUME_FILTERED,
         fiber_density = fibers_corrected / VOLUME_FILTERED,
         bead_density = beads_corrected / VOLUME_FILTERED) %>%
  select(Site, rep, total_microplastics, density, 
         fragment_density, fiber_density, bead_density)

head(microplastics_corrected)

write.csv(x = microplastics_corrected, file = "../cleaned_data/microplastics.csv",
          row.names = FALSE)


# 11. Load and clean calculated distance metrics --------------------------

distance_orig <- read.csv(file = "../original_data/baikal_site_distances_mfm_20180517.csv",
                          header = TRUE)

# Select columns of interest
distance <- distance_orig %>%
  select(-SOUTH_UP, -NORTH_UP, -X, -X.1, -X.2) %>%
  rename(Site = Sample_ID)

head(distance)

write.csv(x = distance, file = "../cleaned_data/distance.csv",
          row.names = FALSE)