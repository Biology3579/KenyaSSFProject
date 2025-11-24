##Exploration of data
##
##

# Load necessary libraries
library(tidyverse) 

# --- 1. Extract Unique Station Names ---
# Use the pull() and unique() functions to extract the unique names 
# as character vectors from the specified column.

# Extract from clean_fish_2018
fish_stations_comoros <- clean_fish_2018 %>%
  filter(country == 'comoros') %>% 
  pull(station) %>%
  unique()
fish_stations_comoros

# Extract from clean_location_2018
location_stations_comoros <- clean_location_2018 %>%
  filter(country == 'comoros') %>% 
  pull(station) %>%
  unique()
location_stations

# Extract from clean_dive_details_2018
dive_stations <- clean_dive_details_2018 %>%
  filter(country == 'comoros') %>% 
  pull(station) %>%
  unique()
dive_stations


# --- 2. Compare the Station Name Sets ---

## Find stations present in ALL three files (The Overlap)
all_overlap <- intersect(fish_stations, intersect(location_stations, dive_stations))

## Find stations UNIQUE to the fish file (Only in fish, not in location or dive)
only_in_fish <- setdiff(fish_stations, union(location_stations, dive_stations))

## Find stations that are MISSING from the fish file (Present in location OR dive, but not fish)
missing_from_fish <- setdiff(union(location_stations, dive_stations), fish_stations)


# --- 3. View Results (Print Statements) ---

# Total unique stations across all files
total_unique_stations <- union(fish_stations, union(location_stations, dive_stations))

cat("--- Summary of Unique Station Names ---\n")
cat("Total Unique Stations (Combined):", length(total_unique_stations), "\n")
cat("Stations in ALL Three Files (Overlap):", length(all_overlap), "\n\n")

cat("Stations ONLY in clean_fish_2018:", only_in_fish, "\n")
cat("Stations MISSING from clean_fish_2018 (Present elsewhere):", missing_from_fish, "\n")

# To see a table of overlaps, you can create a tibble
overlap_table <- tibble(
  `Source File` = c("Fish", "Location", "Dive Details"),
  `Count` = c(length(fish_stations), length(location_stations), length(dive_stations))
)

print(overlap_table)

## Tanzania

# --- 1. Extract Unique Station Names ---
# Use the pull() and unique() functions to extract the unique names 
# as character vectors from the specified column.

# Extract from clean_fish_2018
fish_stations_tanzania <- clean_fish_2018 %>%
  filter(country == 'tanzania') %>% 
  pull(station) %>%
  unique()
fish_stations_tanzania

# Extract from clean_location_2018
location_stations_tanzania <- clean_location_2018 %>%
  filter(country == 'tanzania') %>% 
  pull(station) %>%
  unique()
location_stations_tanzania

# Extract from clean_dive_details_2018
dive_stations_tanzania <- clean_dive_details_2018 %>%
  filter(country == 'tanzania') %>% 
  pull(station) %>%
  unique()
dive_stations_tanzania


# --- 2. Compare the Station Name Sets ---

## Find stations present in ALL three files (The Overlap)
all_overlap <- intersect(fish_stations, intersect(location_stations, dive_stations))
all_overlap

## Find stations UNIQUE to the fish file (Only in fish, not in location or dive)
only_in_fish <- setdiff(fish_stations, union(location_stations, dive_stations))
only_in_fish

## Find stations that are MISSING from the fish file (Present in location OR dive, but not fish)
missing_from_fish <- setdiff(union(location_stations, dive_stations), fish_stations)
missing_from_fish


# --- 3. View Results (Print Statements) ---

# Total unique stations across all files
total_unique_stations <- union(fish_stations, union(location_stations, dive_stations))

cat("--- Summary of Unique Station Names ---\n")
cat("Total Unique Stations (Combined):", length(total_unique_stations), "\n")
cat("Stations in ALL Three Files (Overlap):", length(all_overlap), "\n\n")

cat("Stations ONLY in clean_fish_2018:", only_in_fish, "\n")
cat("Stations MISSING from clean_fish_2018 (Present elsewhere):", missing_from_fish, "\n")

# To see a table of overlaps, you can create a tibble
overlap_table <- tibble(
  `Source File` = c("Fish", "Location", "Dive Details"),
  `Count` = c(length(fish_stations), length(location_stations), length(dive_stations))
)

print(overlap_table)