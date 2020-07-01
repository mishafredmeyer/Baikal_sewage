# This script is intending to take raw Google Earth outputs and
# calculate inverse distance weighted population for the entirety of
# our study area.


# Load packages -----------------------------------------------------------

library(sf)
library(tidyverse)
library(spdplyr)


# Step 1. Load data -------------------------------------------------------

# Load the shapefiles and metadata containing information of the development
# polygons and the lake shorelines.

baikal_shapefile <- sf::st_read(dsn = "../original_data/Baikal_shapefile.kml")


metadata <- read.csv(file = "../cleaned_data/metadata.csv",
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
  dplyr::select(Site, lat, long)

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
          file = "../cleaned_data/distance_weighted_population_metrics.csv",
          row.names = FALSE)
