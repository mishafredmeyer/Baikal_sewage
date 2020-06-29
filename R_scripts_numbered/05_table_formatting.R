## This script aggregates data to produce tables used within the 
## associated manuscript. 

library(dplyr)
library(tidyr)
library(stringr)


# 1. Load the data --------------------------------------------------------

metadata <- read.csv(file = "../cleaned_data/metadata.csv", 
                     header = TRUE)

ppcp <- read.csv(file = "../cleaned_data/ppcp.csv", 
                 header = TRUE)

nutrients <- read.csv(file = "../cleaned_data/nutrients.csv", 
                      header = TRUE)

microplastics <- read.csv(file = "../cleaned_data/microplastics.csv", 
                          header = TRUE)

distance <- read.csv(file = "../cleaned_data/distance_weighted_population_metrics.csv", 
                     header = TRUE)


# 2. Site metadata table --------------------------------------------------

metadata_formatted <- metadata %>%
  select(-mid_temp, -bottom_temp) %>%
  rename(Latitude = lat,
         Longitude = long,
         Depth_m = depth,
         Distance_to_shore_m = dist_to_shore,
         Air_Temperature_C = air_temp,
         Surface_Temperature_C = surface_temp) %>%
  mutate(Adjacent_Population = ifelse(test = Site %in% c("BGO-2", "BGO-3"), 
                                      yes = 600, no = 0),
         Adjacent_Population = ifelse(test = Site %in% c("BK-1", "BK-2", "BK-3"), 
                                      yes = 200, no = Adjacent_Population),
         Adjacent_Population = ifelse(test = Site %in% c("LI-1", "LI-2", "LI-3"), 
                                      yes = 2000, no = Adjacent_Population),
         Adjacent_Population = ifelse(test = Site %in% c("OS-1", "OS-2", "OS-3"), 
                                      yes = NA, no = Adjacent_Population))

## This table is meant to be associated with Table 1 in the main 
## body of the manuscript. 

write.csv(x = metadata_formatted, file = "../tables/metadata_table1.csv", 
          row.names = FALSE)

# 3. PPCP table -----------------------------------------------------------

ppcp_formatted <- ppcp %>%
  rename(Total_PPCP_Concentration = ppcp_sum)

write.csv(x = ppcp_formatted, file = "../tables/ppcp.csv", 
          row.names = FALSE)


# 4. Nutrients table ------------------------------------------------------

nutrients_formatted <- nutrients %>%
  select(-mean_TP_mg_dm3) %>%
  rename(NH4_mg_dm3 = mean_NH4_mg_dm3,
         NO3_mg_dm3 = mean_NO3_mg_dm3, 
         PO4_mg_dm3 = mean_TPO43_mg_dm3)

write.csv(x = nutrients_formatted, file = "../tables/nutrients.csv", 
          row.names = FALSE)


# 5. Microplastics table --------------------------------------------------

microplastics_formatted <- microplastics %>%
  group_by(Site) %>%
  summarize(mean_microplastic_density = mean(density),
            mean_fragment_density = mean(fragment_density),
            mean_fiber_density = mean(fiber_density),
            mean_bead_density = mean(bead_density)) %>%
  rename(Microplastic_density_microplastics_per_L = mean_microplastic_density, 
         Fragment_density_microplastics_per_L = mean_fragment_density,
         Fiber_density_microplastics_per_L = mean_fiber_density, 
         Bead_density_microplastics_per_L = mean_bead_density)

write.csv(x = microplastics_formatted, file = "../tables/microplastics.csv", 
          row.names = FALSE)


# 6. Combined table -------------------------------------------------------

low <- c("BGO-1", "BGO-2", "BGO-3", "KD-1", "KD-2", "MS-1", "OS-1")
mod <- c( "BK-2", "BK-3", "SM-1", "OS-3")
high <- c("BK-1", "EM-1", "LI-3", "LI-2", "LI-1", "OS-2")

## This table is associated with Table 3 in the manuscript

meta_nut_ppcp_mp <- full_join(metadata_formatted, nutrients_formatted) %>%
  full_join(., ppcp_formatted) %>%
  full_join(., distance) %>%
  mutate(Categorical_distance_weighted_population = ifelse(test = Site %in% low, 
                                                           yes = "Low", no = NA),
         Categorical_distance_weighted_population = ifelse(test = Site %in% mod, 
                                                           yes = "Mod", 
                                                           no = Categorical_distance_weighted_population),
         Categorical_distance_weighted_population = ifelse(test = Site %in% high, 
                                                           yes = "High", 
                                                           no = Categorical_distance_weighted_population)) %>%
  select(Site, NH4_mg_dm3:Cotinine, distance_weighted_population, Categorical_distance_weighted_population)

write.csv(x = meta_nut_ppcp_mp, file = "../tables/combined_table3.csv", 
          row.names = FALSE)

