## This script reads in the raw data and formats the data such that misspelled 
## taxa are correctly spelled, column headers are uniform, and data can be 
## shared easily with collaborators. This script is meant to prime the data 
## for script 01_data_cleaning.R, where data are averaged and prepped for 
## the appropriate analysis in scripts 02-07. 


# 1. Load packages --------------------------------------------------------

library(tidyverse)
library(lubridate)
library(stringr)
library(janitor)
library(sf)
library(spdplyr)


# 2. Load and clean PPCP data ---------------------------------------------

# Raw data
ppcp_orig <- read.csv(file = "../original_data/PPCP_Baikal_orig_20180524.csv",
                      header = TRUE)

# Select variables and sites of interest, sum PPCPs
ppcp <- ppcp_orig %>%
  clean_names() %>%
  filter(sample_id == "LI-1" |
           sample_id == "LI-2" |
           sample_id == "LI-3" |
           sample_id == "BK-1" |
           sample_id == "BK-2" |
           sample_id == "BK-3" |
           sample_id == "BGO-1" |
           sample_id == "BGO-2" |
           sample_id == "BGO-3" |
           sample_id == "KD-1" |
           sample_id == "KD-2" |
           sample_id == "EM-1" |
           sample_id == "MS-1" |
           sample_id == "SM-1" |
           sample_id == "OS-1" |
           sample_id == "OS-2" |
           sample_id == "OS-3") %>%
  group_by(sample_id) %>%
  rename(paraxanthine = x1_7_dimethylxanthine, 
         site = sample_id) %>%
  mutate(collection_year = year(mdy(collection_date)),
         collection_month = month(mdy(collection_date)),
         collection_day = day(mdy(collection_date)),
         analysis_year = year(mdy(analysis_date)),
         analysis_month = month(mdy(analysis_date)),
         analysis_day = day(mdy(analysis_date)))

# Take a look
head(ppcp)

# Export new version of the data
write.csv(x = ppcp, file = "../clean_disaggregated_data/ppcp.csv",
          row.names = FALSE)


# 2. Load and clean nutrient data -----------------------------------------

nutrients_orig <- read.csv(file = "../original_data/baikal_nearshore_nutrient_data_20151009.csv",
                           header = TRUE)

# Take nutrients averages by site
nutrients <- nutrients_orig %>%
  clean_names(case = "snake") %>%
  rename(site = sample)

head(nutrients)

write.csv(x = nutrients, file = "../clean_disaggregated_data/nutrients.csv",
          row.names = FALSE)


# 3. Load and clean chlorophyll a data ------------------------------------

chla_orig <- read.csv(file = "../original_data/chlorophyll_20170117.csv",
                      header = TRUE)

# Select sites of interest and take average chl a by site
chlorophylla <- chla_orig %>%
  select(-Chlorophyll) %>%
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
  clean_names(case = "snake") %>%
  # Correct for incorrect date formatting
  mutate(replicate = as.character(replicate),
         replicate = ifelse(test = grepl(pattern = "Jan", x = replicate), 
                            yes = 1, no = replicate),
         replicate = ifelse(test = grepl(pattern = "Feb", x = replicate), 
                            yes = 2, no = replicate),
         replicate = ifelse(test = grepl(pattern = "Mar", x = replicate), 
                            yes = 3, no = replicate)) %>%
  rename(site = station)

head(chlorophylla)

write.csv(x = chlorophylla, file = "../clean_disaggregated_data/chlorophylla.csv",
          row.names = FALSE)


# 4. Load and clean lat/long metadata -------------------------------------

metadata_orig <- read.csv(file = "../original_data/baikal_nearshore_metadata_201508.csv",
                          header = TRUE)

# Select columns of interest
metadata <- metadata_orig %>%
  rename(common_locaiton = location,
         local_site = site,
         site = loc_site,
         site_description = site_desc,
         depth_m = depth,
         distance_to_shore = dist_to_shore,
         air_temp_celsius = air_temp,
         surface_temp_celsius = surface_temp,
         mid_temp_celsius = mid_temp,
         bottom_temp_celsius = bottom_temp) %>%
  select(year, month, day, time, site, lat, long, site_description, distance_to_shore, 
         air_temp_celsius, surface_temp_celsius, mid_temp_celsius, bottom_temp_celsius, 
         comments, shore_photo, substrate_photo, sponges, brandtia)

head(metadata)

write.csv(x = metadata, file = "../clean_disaggregated_data/metadata.csv",
          row.names = FALSE)


# 5. Load and clean macroinbertebrate data --------------------------------

inverts_orig <- read.csv(file = "../original_data/macroinvert_community_QAQC_mfm_20171108.csv",
                         header = TRUE)

# Clean up naming in data and add up counts by taxon
# During enumeration, choronomids and and oligochaetes were poorly preserved, 
# so we resolved after a few samples were counted to remove them from the dataset. 
# Additionally, hyallela_cziarnianski was misidentified in some samples. Rather than
# potentially erroneously include those counts, we remove this taxon. 
# LI-3 was partitioned into two samples due to the high number of individuals. 
# This step combines that replicate into one row. 
inverts_summarized <- inverts_orig %>%
  select(-X) %>%
  gather(key = Site, value = Count, MS1.3:BK1.3) %>%
  rename(Taxon = Invertebrate) %>%
  mutate(Site = gsub(pattern = ".", replacement = "_", x = Site, fixed = TRUE),
         Taxon = gsub(pattern = " ", replacement = "_", x = Taxon),
         Count = ifelse(test = is.na(Count), yes = 0, no = Count)) %>%
  separate(col = Site, into = c("Location", "Replicate", "Duplicate"),
           remove = FALSE) %>%
  filter(!(Taxon %in% c("Propapaidae", "choronomids_", "hyallela_cziarnianski_"))) %>%
  mutate(Taxon = ifelse(test = grepl(pattern = "Brandtia_latissima", x = Taxon), 
                                     yes = "Brandtia_latissima", no = Taxon),
         Replicate = ifelse(test = Replicate %in% c("1B1", "1B2"),
                            yes = "1", no = Replicate)) %>%
  group_by(Location, Replicate, Duplicate, Taxon) %>%
  summarize(sum_Count = sum(Count)) %>%
  ungroup()

# Take mean counts by taxon, flesh out taxonomic info, spread to wide format
inverts_wide <- inverts_summarized %>%
  group_by(Location, Taxon, Replicate) %>%
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
                        yes = "Flatworms", no = Taxon)) %>%
  spread(key = Taxon, value = mean_Count) %>%
  select(-Total) %>%
  separate(col = Location, into = c("Location", "Number"), sep = -1) %>%
  unite(col = "Site", Location, Number, sep = "-")

head(inverts_wide)

write.csv(x = inverts_wide, file = "../clean_disaggregated_data/invertebrates.csv",
          row.names = FALSE)


# 6. Load and clean periphyton data ---------------------------------------

periphyton_orig <- read.csv(file = "../original_data/periphyton_20180917.csv",
                            header = TRUE)

# Make long format, take mean counts by taxon
periphyton_wide <- periphyton_orig %>%
  rename("Site" = "site",
         "Replicate" = "rep",
         "subsamples_counted" = "counts") %>%
  separate(col = Site, into = c("Location", "Number"), sep = -1) %>%
  unite(col = "Site", Location, Number, sep = "-") %>%
  select(-date) %>%
  clean_names()

head(periphyton_wide)

write.csv(x = periphyton_wide, file = "../clean_disaggregated_data/periphyton.csv",
          row.names = FALSE)


# 7. Load and clean stable isotope data -----------------------------------

stable_isotopes_orig <- read.csv(file = "../original_data/sia_results_mfm_20170509.csv",
                                 header = TRUE)

# Parse identifier data into site and taxonomic data
# In this analysis, "Splash zone" was used as the shorthand 
# for periphyton. We relabel "Splash zone" as periphyton as such. 
stable_isotopes <- stable_isotopes_orig %>%
  separate(col = Identifier, into = c("Site", "Genus", "Species"), sep = " ") %>%
  mutate(Genus = ifelse(test = Genus == "E.",
                        yes = "Eulimnogammarus", Genus),
         Genus = ifelse(test = Genus == "P.",
                        yes = "Pallasea", no = Genus),
         Genus = ifelse(test = Genus == "Sp.",
                        yes = "Periphyton", no = Genus),
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
         Species = ifelse(test = Species == "zone",
                          yes = NA, no = Species))

head(stable_isotopes)

write.csv(x = stable_isotopes, file = "../clean_disaggregated_data/stable_isotopes.csv",
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
         Genus = ifelse(test = Genus == "Spl",
                        yes = "Periphyton", no = Genus),
         Species = ifelse(test = Genus == "Eulimnogammarus" & Species == "ver",
                          yes = "verrucosus", no = Species),
         Species = ifelse(test = Genus == "Eulimnogammarus" & Species == "cyan",
                          yes = "cyaneus", no = Species),
         Species = ifelse(test = Genus == "Pallasea" & Species == "can",
                          yes = "cancellus", no = Species),
         Species = ifelse(test = Species == "zone",
                          yes = NA, no = Species)) %>%
  rename(Site = location)

head(fatty_acid)

write.csv(x = fatty_acid, file = "../clean_disaggregated_data/fatty_acid.csv",
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
         Genus = ifelse(test = Genus == " Spl",
                        yes = "Periphyton", no = Genus),
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
         Species = ifelse(test = Species == "zone",
                          yes = NA, no = Species)) %>%
  rename(total_lipid_mg_per_g = total.lipid.mg.g)

head(total_lipid)

write.csv(x = total_lipid, file = "../clean_disaggregated_data/total_lipid.csv",
          row.names = FALSE)


# 10. Load and clean microplastics data -----------------------------------

microplastics_orig <- read.csv(file = "../original_data/microplastics_mfm_20171010.csv",
                               header = TRUE)

# Here we convert the volume per filtration and the number of times
# that volume was filtered for a total volume filtered. 

microplastics <- microplastics_orig %>%
  select(-date) %>%
  unite(col = "Site", location, site, sep = "-") %>%
  mutate(volume_filtered_mL = volume * volume_rep) %>%
  select(-volume, -volume_rep)

head(microplastics)

write.csv(x = microplastics, file = "../clean_disaggregated_data/microplastics.csv",
          row.names = FALSE)


# 11. Load and clean calculated distance metrics --------------------------

# Step 1. Load data -------------------------------------------------------

# Load the shapefiles and metadata containing information of the development
# polygons and the lake shorelines.

baikal_shapefile <- sf::st_read(dsn = "../clean_disaggregated_data/Baikal_shapefile.kml")


metadata <- read.csv(file = "../clean_disaggregated_data/metadata.csv",
                     header = TRUE, stringsAsFactors = FALSE)


# Step 2: Extract and format development area data ------------------------

loc_areas <- baikal_shapefile %>%
  filter(!grepl("shoreline", Name)) %>%
  st_area() %>%
  enframe(name = NULL) %>%
  rename(development_area_m2 = value) %>%
  mutate(development_area_m2 = as.numeric(development_area_m2),
         development_area_km2 = development_area_m2 / 1000000) %>%
  cbind(baikal_shapefile$Name) %>%
  rename(Site = `baikal_shapefile$Name`) %>%
  mutate(Site = as.character(Site),
         Site = ifelse(test = Site == "Bolshoe Goloustnoe",
                       yes = "BGO", no = Site),
         Site = ifelse(test = Site == "Bolshie Koty",
                       yes = "BK", no = Site),
         Site = ifelse(test = Site == "Listvyanka",
                       yes = "LI", no = Site)) %>%
  filter(Site %in% c("BK", "LI", "BGO"))


# Step 3: Extract Development names ---------------------------------------

shoreline_names <- baikal_shapefile %>%
  filter(grepl(pattern = "shoreline", x = Name)) %>%
  as_tibble()


# Step 4: Extract shoreline length info and join with area ----------------

loc_shoreline_area_length <- baikal_shapefile %>%
  filter(grepl("shoreline", Name)) %>%
  st_length() %>%
  as_tibble() %>%
  rename(development_shoreline_length_m = value) %>%
  mutate(development_shoreline_length_m = as.numeric(development_shoreline_length_m),
         development_shoreline_length_km = development_shoreline_length_m / 1000) %>%
  cbind(shoreline_names$Name) %>%
  rename(Site = `shoreline_names$Name`) %>%
  mutate(Site = as.character(Site),
         Site = ifelse(grepl("Bolshoe Goloustnoe", Site), "BGO", Site),
         Site = ifelse(grepl("Bolshie Koty", Site), "BK", Site),
         Site = ifelse(grepl("Listvyanka", Site), "LI", Site)) %>%
  filter(Site %in% c("BK", "LI", "BGO")) %>%
  full_join(x = ., y = loc_areas, by = c("Site"))


# Step 5: Calculate centroids for each developed site ---------------------

baikal_shapefile$centroids <- st_centroid(x = baikal_shapefile)


# Step 6: Convert metadata into a spatial object --------------------------

site_loc <- metadata %>%
  dplyr::select("Site" = site, lat, long)

site_loc_pts <- st_as_sf(site_loc[, 2:3], coords = c("long", "lat"),
                         crs = 4326)


# Step 7: Join all data together and format into a dataframe --------------

locs_centroids <- st_distance(x = site_loc_pts,
                              y = baikal_shapefile$centroids$geometry[1:3]) %>%
  as_tibble() %>%
  rename(BGO = V1,
         BK = V2,
         LI = V3) %>%
  mutate(BGO = as.numeric(BGO),
         BK = as.numeric(BK),
         LI = as.numeric(LI)) %>%
  cbind(., site_loc) %>%
  gather(key = nearest_neighbor, value = distance, BGO:LI) %>%
  group_by(Site, lat, long) %>%
  mutate(distance_km = distance / 1000,
         population = ifelse(test = nearest_neighbor == "BGO",
                             yes = 600, no = NA),
         population = ifelse(test = nearest_neighbor == "BK",
                             yes = 80, no = population),
         population = ifelse(test = nearest_neighbor == "LI",
                             yes = 5000, no = population)) %>%
  left_join(x = ., y = loc_shoreline_area_length,
            by = c("nearest_neighbor" = "Site")) %>%
  mutate(distance_weighted_population = ((population * development_shoreline_length_km) /
                                           development_area_km2) / distance_km) %>%
  group_by(Site) %>%
  summarize(distance_weighted_population = sum(distance_weighted_population)) %>%
  arrange(distance_weighted_population)

write.csv(x = locs_centroids,
          file = "../clean_disaggregated_data/distance_weighted_population_metrics.csv",
          row.names = FALSE)
