## Checking name consistency in the data


# Load necessary libraries
library(tidyverse) 

# Load data
fish_data_2009 <- readr::read_rds(here::here("processed_data", "clean_fish_2009.rds"))
locations_2009 <- readr::read_rds(here::here("processed_data", "clean_location_2009.rds"))
dive_deets_2009 <- readr::read_rds(here::here("processed_data", "clean_dive_details_2009.rds"))
sst_2009 <- read.csv(here::here("processed_data", "locations_with_sst_2009.csv"))
chla_2009 <- read.csv(here::here("processed_data", "locations_with_chla_2009.csv"))
h_gravity_2009 <- read.csv(here::here("processed_data", "locations_with_gravity.csv"))

# --- 1. Extract Unique Station Names ---
# Use the pull() and unique() functions to extract the unique names 
# as character vectors from the specified column.

# Extract from clean_fish_2018
fish_stations <- fish_data_2009 %>%
  pull(station) %>%
  unique()
fish_stations

# Extract from clean_location_2018
location_stations <- locations_2009 %>%
  pull(station) %>%
  unique()
location_stations

# Extract from clean_dive_details_2018
dive_stations <- dive_deets_2009%>%
  filter(country == 'comoros') %>% 
  pull(station) %>%
  unique()
dive_stations


## Find stations only in one of the files
only_in_fish <- setdiff(fish_stations, union(location_stations, dive_stations))
only_in_location <- setdiff(location_stations, union(fish_stations, dive_stations))
only_in_dive <- setdiff(dive_stations, union(location_stations, fish_stations))

length(only_in_fish)
only_in_fish
length(only_in_location)
only_in_location
length(only_in_dive)



# Extract Unique Site Names ---
# Use the pull() and unique() functions to extract the unique names 
# as character vectors from the specified column.

# Extract from clean_fish_2018
fish_sites <- fish_data_2009 %>%
  pull(site) %>%
  unique()
fish_sites

# Extract from clean_location_2018
location_sites <- locations_2009 %>%
  pull(site) %>%
  unique()
location_sites

# Extract from clean_dive_details_2018
dive_sites <- dive_deets_2009%>%
  filter(country == 'comoros') %>% 
  pull(site) %>%
  unique()
dive_sites


## Find sites only in one of the files
only_in_fish <- setdiff(fish_sites, union(location_sites, dive_sites))
only_in_location <- setdiff(location_sites, union(fish_sites, dive_sites))
only_in_dive <- setdiff(dive_sites, union(location_sites, fish_sites))

length(only_in_fish)
only_in_fish
length(only_in_location)
only_in_location
length(only_in_dive)
only_in_dive

# # Extract Unique Location Names ---
# Use the pull() and unique() functions to extract the unique names 
# as character vectors from the specified column.

# Extract from clean_fish_2018
fish_locations <- fish_data_2009 %>%
  pull(location) %>%
  unique()
fish_locations


# Extract from clean_location_2018
location_locations <- locations_2009 %>%
  pull(location) %>%
  unique()
location_locations

# Extract from clean_dive_details_2018
dive_locations <- dive_deets_2009 %>%
  filter(country == 'comoros') %>% 
  pull(location) %>%
  unique()
dive_locations


## Find locations only in one of the files
only_in_fish <- setdiff(fish_locations, union(location_locations, dive_locations))
only_in_location <- setdiff(location_locations, union(fish_locations, dive_locations))
only_in_dive <- setdiff(dive_locations, union(location_locations, fish_locations))

length(only_in_fish)
only_in_fish
length(only_in_location)
only_in_location
length(only_in_dive)
only_in_dive

