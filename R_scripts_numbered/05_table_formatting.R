library(dplyr)
library(reshape2)
library(lubridate)
library(stringr)


# 1. Load the data --------------------------------------------------------

metadata <- read.csv("../cleaned_data/metadata_20190320.csv", header = TRUE)
ppcp <- read.csv("../cleaned_data/ppcp_20190320.csv", header = TRUE)
nutrients <- read.csv("../cleaned_data/nutrients_20190320.csv", header = TRUE)
microplastics <- read.csv("../cleaned_data/microplastics_20190320.csv", header = TRUE)
stable_isotopes <- read.csv("../cleaned_data/stable_isotopes_20190320.csv", header = TRUE)
distance <- read.csv("../cleaned_data/distance_20190320.csv", header = TRUE)


# 2. Site metadata table --------------------------------------------------

metadata_formatted <- metadata %>%
  rename(Latitude = lat,
         Longitude = long,
         Depth_m = depth,
         Distance_to_shore_m = dist_to_shore,
         Air_Temperature_C = air_temp,
         Surface_Temperature_C = surface_temp,
         Midpoint_Temperature_C = mid_temp,
         Bottom_Temperature_C = bottom_temp) %>%
  mutate(Adjacent_Population = ifelse(Site %in% c("BGO-2", "BGO-3"), 200, 0),
         Adjacent_Population = ifelse(Site %in% c("BK-1", "BK-2", "BK-3"), 80, Adjacent_Population),
         Adjacent_Population = ifelse(Site %in% c("LI-1", "LI-2", "LI-3"), 3000, Adjacent_Population),
         Adjacent_Population = ifelse(Site %in% c("OS-1", "OS-2", "OS-3"), NA, Adjacent_Population))
  
write.csv(metadata_formatted, "../tables/metadata_table1.csv", row.names = FALSE)


# 3. PPCP table -----------------------------------------------------------

ppcp_formatted <- ppcp %>%
  rename(Total_PPCP_Concentration = PPCP.SUM)

write.csv(ppcp_formatted, "../tables/ppcp_20190405.csv", row.names = FALSE)


# 4. Nutrients table ------------------------------------------------------

nutrients_formatted <- nutrients %>%
  select(-mean_TP_mg_dm3) %>%
  rename(NH4_mg_dm3 = mean_NH4_mg_dm3,
         NO3_mg_dm3 = mean_NO3_mg_dm3, 
         PO4_mg_dm3 = mean_TPO43_mg_dm3)

write.csv(nutrients_formatted, "../tables/nutrients_20190405.csv", row.names = FALSE)


# 5. Microplastics table --------------------------------------------------

microplastics_formatted <- microplastics %>%
  select(Site, mean_microplastic_density, mean_fragment_density, mean_fiber_density, 
         mean_bead_density) %>%
  rename(Microplastic_density_microplastics_per_L = mean_microplastic_density, 
         Fragment_density_microplastics_per_L = mean_fragment_density,
         Fiber_density_microplastics_per_L = mean_fiber_density, 
         Bead_density_microplastics_per_L = mean_bead_density)

write.csv(microplastics_formatted, "../tables/microplastics_20190405.csv", row.names = FALSE)


# 6. Combined table -------------------------------------------------------

low <- c("BGO-1", "BGO-2", "KD-1", "KD-2")
mod <- c("BGO-3", "BK-2", "BK-3", "MS-1")
high <- c("BK-1", "SM-1", "EM-1", "LI-3", "LI-2")

meta_nut_ppcp_mp <- full_join(metadata_formatted, nutrients_formatted) %>%
  full_join(., ppcp_formatted) %>%
  full_join(., microplastics_formatted) %>%
  full_join(., distance) %>%
  mutate(Distance_weighted_population = (SOUTH_SHORE_DIST*(POP_SOUTH_DEV/SOUTH_SHORE_AREA))/SOUTH_DEV_DIST,
         Categorical_distance_weighted_population = ifelse(Site %in% low, "Low", NA),
         Categorical_distance_weighted_population = ifelse(Site %in% mod, "Mod", Categorical_distance_weighted_population),
         Categorical_distance_weighted_population = ifelse(Site %in% high, "High", Categorical_distance_weighted_population)) %>%
  select(Site, NH4_mg_dm3:Cotinine, Distance_weighted_population, Categorical_distance_weighted_population)

write.csv(meta_nut_ppcp_mp, "../tables/combined_table2.csv", row.names = FALSE)


# 7. Stable Isotopes ------------------------------------------------------

stable_isotopes_formatted <- stable_isotopes %>%
  unite("Taxon", Genus, Species) %>%
  gather(Isotope, Value, C13:N15) %>%
  group_by(Taxon, Isotope) %>%
  summarize(min = min(Value),
            Q1 = quantile(Value, probs = 0.25),
            Mean = mean(Value),
            Median = median(Value),
            Q3 = quantile(Value, probs = 0.75),
            max = max(Value)) 

write.csv(stable_isotopes_formatted, "../tables/stable_isotopes_individual_table.csv", row.names = FALSE)
