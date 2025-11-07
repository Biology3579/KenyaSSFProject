## ---------------------------
##
## Script name: cleaning_data.r
##
## Purpose of script: A file to prepare the raw_data for analysis
##
## Author: Candela Ferrer Díez
##
## Date Created: 26-10-2025 
##
##
## ---------------------------
##
## Notes:
##   Use renv::restore() to restore packages
##   The raw data involves data from different sites in Kenya monitored between 2018-2023.  
##   In this script I will clean and process this data ready for analysis
##
## ---------------------------



# ---------------------------
# LIBRARIES AND FUNCTIONS
# ---------------------------
library(tidyverse)  # Data manipulation (dplyr, ggplot2), pipes (%>%), string functions (str_trim, str_to_title)
library(readxl)     # Read Excel files (.xlsx, .xls)
library(here)       # Reproducible file paths that work across different computers
library(janitor)    # Clean column names to lowercase snake_case
library(lubridate)  # Parse and manipulate dates easily
library(dplyr)      # 

source(here("functions", "exploring_cleaning_processing.R")) ## to source own functions from cleaning.R file.

# ---------------------------
# LOAD DATA
# ---------------------------

#Load, filter and save fish census data
(raw_fish_data <- read_excel(
  here("raw_data", "CORDIO Fish_abund_2018_2023 for Laura_19122023.xlsx"),
  sheet = 1,
  col_types = c("text", rep("guess", 14))
) %>%
    mutate(
      Date = suppressWarnings(case_when(
        grepl("^\\d{2}/\\d{2}/\\d{4}$", Date) ~ dmy(Date),
        grepl("^\\d+$", Date) ~ as.Date(as.numeric(Date), origin = "1899-12-30"),
        TRUE ~ as.Date(NA)
      ))
    ) %>%
  # Fix VumaN_01 dates - all should be 15/02/2021
  mutate(Date = if_else(
    Station == "VumaN_01" & Date >= as.Date("2024-02-15") & Date <= as.Date("2026-02-15"),
    as.Date("2021-02-15"),
    Date)) %>% 
    filter(Country == "Kenya") %>%
  write_csv(here("raw_data", "raw_fish_data.csv")))


# Load, filter and save location data
(raw_location <- read_excel(
  here("raw_data", "CORDIO Fish_abund_2018_2023 for Laura_19122023.xlsx"), 
  sheet = 2,
  col_types = c("text", rep("guess", 13))
  ) %>% 
    mutate(
      Date = suppressWarnings(case_when(
        grepl("^\\d{2}/\\d{2}/\\d{4}$", Date) ~ dmy(Date),
        grepl("^\\d+$", Date) ~ as.Date(as.numeric(Date), origin = "1899-12-30"),
        TRUE ~ as.Date(NA)
      ))
  ) %>%
     # Manually fix the IWE_01 date - change 3022 to 2022
  mutate(Date = if_else(Station == "IWE_01" & year(Date) == 3022,
                        as.Date("2022-03-04"),
                        Date)) %>% 
    filter(Country == "Kenya") %>%
  write_csv(here("raw_data", "raw_location.csv")))

# Load, filter and diver details data
(raw_dive_details <- read_excel(
  here("raw_data", "CORDIO Fish_abund_2018_2023 for Laura_19122023.xlsx"), 
  sheet = 3,
  col_types = c("text", rep("guess", 14))) %>% 
    mutate(
      Date = suppressWarnings(case_when(
        grepl("^\\d{2}/\\d{2}/\\d{4}$", Date) ~ dmy(Date),
        grepl("^\\d+$", Date) ~ as.Date(as.numeric(Date), origin = "1899-12-30"),
        TRUE ~ as.Date(NA)))) %>%
    # Manually fix the IWE_01 date - change 3022 to 2022
    mutate(Date = if_else(Station == "IWE_01" & year(Date) == 3022,
                          as.Date("2022-03-04"),
                          Date)) %>% 
    filter(Country == "Kenya") %>%
  write_csv(here("raw_data", "raw_dive_details.csv")))

# ---------------------------
# EXPLORE DATA
# ---------------------------

# Explore all datasets
explore_data(raw_location, "LOCATION")
explore_data(raw_fish_data, "FISH ABUNDANCE")
explore_data(raw_dive_details, "DIVE DETAILS")

# ---------------------------
# CLEAN DATA
# ---------------------------

# Clean location data
clean_location <- raw_location %>% 
  clean_names() %>% 
  standardize_site_format() %>% 
  mutate(site = case_when(
    # Fix abbreviations
    site == "kilifin" ~ "kilifi_north",
    
    # Fix spelling variants
    site == "diani_reef_2" ~ "diani_reef",
    
    # Fix Coral Garden variants
    site == "coral_garden_malindi" ~ "malindi_coral_garden",
    site == "coral_garden" ~ "mombasa_coral_garden",
    site == "coral_garden_watamu" ~ "watamu_coral_garden",
    
    TRUE ~ site
  )) %>% 
  # Remove sites not monitored (using the cleaned name 'richard_bennet')
  filter(site != "richard_bennet") %>% 
  clean_data(
    # Date is a character string in DD/MM/YYYY format, so use 'dmy'
    date_cols = NULL,
    numeric_cols = c(
      "latitude", 
      "longitude", 
      "exposure", 
      "slope" # Treat as numeric, as it represents a degree measurement
    ),
    
    factor_cols = c(
      "country", 
      "site", 
      "station", 
      "location", 
      "orientation", 
      "geology", 
      "reef_geomorphology", # Now cleaned by clean_names()
      "reef_type"           # Now cleaned by clean_names()
    ))

# Clean fish abundance data
clean_fish_data <- raw_fish_data %>% 
  clean_names() %>% 
  standardize_site_format() %>% 
  mutate(site = case_when(
    site == "andromache" ~ "androumache",
    site == "new_coral_garden" ~ "malindi_coral_garden",
    site == "coral_garden" & str_detect(station, "WCG") ~ "watamu_coral_garden",
    site == "coral_garden" & str_detect(station, "CRG") ~ "mombasa_coral_garden",
    site == "kisite_east" ~ "kisite",
    TRUE ~ site
  )) %>% 
  clean_data(
    date_cols = NULL,
    numeric_cols = c(
      "transect_no", # Treating this as a numeric ID for now, though Factor is also common
      "number", 
      "size_class",
      "l", 
      "a", 
      "b", 
      "weight_g", 
      "total_wt_g"),
    factor_cols = c(
      "station", 
      "site", 
      "country", 
      "family", 
      "species", 
      "trophic_groups"))

# Clean dive details
clean_dive_details <- raw_dive_details %>% 
  clean_names() %>% 
  standardize_site_format() %>% 
  mutate(site = case_when(
    # Fix spelling/name writing problems
    site == "vuma" ~ "wembe",
    site == "coral_garden_malindi" ~ "malindi_coral_garden",
    site == "coral_garden_watamu" ~ "watamu_coral_garden",
    site == "coral_garden_mombasa" ~ "mombasa_coral_garden",
    site == "kilifi_plantation" ~ "takaungu",
    site == "jumba" ~ "jumbo_ruins",
    
    TRUE ~ site
  )) %>% 
  # Remove sites not monitored (using the cleaned name 'richard_bennet')
  filter(site != "richard_bennet") %>% 
  clean_data(
    # Date is a character string in DD/MM/YYYY format, so use 'dmy'
    date_cols = NULL,
    
    numeric_cols = c(
      "latitude", 
      "longitude", 
      "rugosity", 
      "viz",             # Should be numeric, even though raw is <chr>
      "temp",            # Should be numeric, even though raw is <chr>
      "max_depth_m",     # Should be numeric, even though raw is <chr>
      "min_depth_m",     # Should be numeric, even though raw is <chr>
      "div_maxd",        # Assuming this is another depth measure
      "total_dive_time_minutes" # Now cleaned by clean_names()
    ),
    
    factor_cols = c(
      "station", 
      "site", 
      "location", 
      "country",
      "habitat"
    ))

# ---------------------------
# VERIFY CLEANING
# ---------------------------

cat("\n========================================\n")
cat("CLEANED DATA SUMMARY\n")
cat("========================================\n\n")

cat("Location data:\n")
print(summary(clean_location))

cat("\nFish data:\n")
print(summary(clean_fish_data))

cat("\nDive details:\n")
print(summary(clean_dive_details))

# Check site consistency after cleaning
cat("\n--- Sites after cleaning ---\n")
cat("Unique sites in location:", length(unique(clean_location$site)), "\n")
cat("Unique sites in fish:", length(unique(clean_fish_data$site)), "\n")
cat("Unique sites in dive:", length(unique(clean_dive_details$site)), "\n")

# ---------------------------
# SAVE CLEANED DATA
# ---------------------------

write_csv(clean_location, here("cleaned_data", "clean_location.csv"))
write_csv(clean_fish_data, here("cleaned_data", "clean_fish_data.csv"))
write_csv(clean_dive_details, here("cleaned_data", "clean_dive_details.csv"))

cat("\n✓ Cleaned data saved to cleaned_data/ folder\n")



#Load Libraries
library(tidyverse) #for data manipulation and visualization
library(readxl) # To read an excel file (within tidyverse)
library(here) # Reproducible working directories

# Load, filter and save fish abundance data


# ---------------------------
# Explore the data
# ---------------------------

# Fish abundance data exploration
cat("\n=== FISH DATA OVERVIEW ===\n\n")
cat("\nUnique station in fish abundance data:\n")
print(unique(raw_fish_data$Station))
cat("Total unique stations:", length(unique(raw_fish_data$Station)), "\n\n")


# Location data exploration
cat("=== LOCATION DATA OVERVIEW ===\n\n")
cat("\nUnique station in location data:\n")
print(unique(raw_location$Station))
cat("Total unique stations:", length(unique(raw_location$Station)), "\n\n")

# Diver details data exploration
cat("=== DIVE DETAILS DATA OVERVIEW ===\n\n")
cat("\nUnique station in location data:\n")
print(unique(raw_dive_details$Station))
cat("Total unique stations:", length(unique(raw_dive_details$Station)), "\n\n")

# ---------------------------
# Check station name consistency between dataframes
# ---------------------------
cat("\n=== STATION NAME CONSISTENCY CHECK ===\n\n")

# Get unique Stations from each dataframe
Stations_location <- unique(raw_location$Station)
Stations_fish <- unique(raw_fish_data$Station)
Stations_dive <- unique(raw_dive_details$Station)

Stations_location
Stations_fish
Stations_dive

# Stations in location but not in fish abundance
Stations_only_location <- setdiff(Stations_location, Stations_fish)
cat("Stations in LOCATION data but NOT in FISH ABUNDANCE data:\n")
if(length(Stations_only_location) > 0) {
  print(Stations_only_location)
} else {
  cat("None\n")
}

# Stations in fish abundance but not in location
Stations_only_fish <- setdiff(Stations_fish, Stations_location)
cat("\nStations in FISH ABUNDANCE data but NOT in LOCATION data:\n")
if(length(Stations_only_fish) > 0) {
  print(Stations_only_fish)
} else {
  cat("None\n")
}

# Stations in dive details but not in location
Stations_only_dive <- setdiff(Stations_dive, Stations_location)
cat("\nStations in DIVE DETAILS data but NOT in LOCATION data:\n")
if(length(Stations_only_dive) > 0) {
  print(Stations_only_dive)
} else {
  cat("None\n")
}

# Stations in all three dataframes
Stations_all_three <- Reduce(intersect, list(Stations_location, Stations_fish, Stations_dive))
cat("\nStations in ALL THREE dataframes:\n")
print(Stations_all_three)
cat("Total Stations in all three:", length(Stations_all_three), "\n")

# Stations in at least two dataframes
Stations_location_fish <- intersect(Stations_location, Stations_fish)
Stations_location_dive <- intersect(Stations_location, Stations_dive)
Stations_fish_dive <- intersect(Stations_fish, Stations_dive)

cat("\nStations in LOCATION and FISH (but may not be in DIVE):", length(Stations_location_fish), "\n")
cat("Stations in LOCATION and DIVE (but may not be in FISH):", length(Stations_location_dive), "\n")
cat("Stations in FISH and DIVE (but may not be in LOCATION):", length(Stations_fish_dive), "\n")

# Summary
cat("\n=== SUMMARY ===\n")
cat("Location data Stations:", length(Stations_location), "\n")
cat("Fish abundance data Stations:", length(Stations_fish), "\n")
cat("Dive details data Stations:", length(Stations_dive), "\n")
cat("Stations in all three:", length(Stations_all_three), "\n")
cat("Stations only in location:", length(Stations_only_location), "\n")
cat("Stations only in fish abundance:", length(Stations_only_fish), "\n")
cat("Stations only in dive details:", length(Stations_only_dive), "\n")

# Check for potential naming issues (e.g., extra spaces, case differences)
if(length(Stations_only_location) > 0 | length(Stations_only_fish) > 0 | length(Stations_only_dive) > 0) {
  cat("\n=== CHECKING FOR NAMING ISSUES ===\n")
  cat("Checking for case sensitivity, whitespace, or spelling differences...\n\n")
  
  # Check for case-insensitive matches
  all_Stations <- unique(c(Stations_location, Stations_fish, Stations_dive))
  all_Stations_lower <- tolower(all_Stations)
  duplicates_when_lowercase <- all_Stations[duplicated(all_Stations_lower) | duplicated(all_Stations_lower, fromLast = TRUE)]
  
  if(length(duplicates_when_lowercase) > 0) {
    cat("Potential case sensitivity issues found:\n")
    print(duplicates_when_lowercase)
  } else {
    cat("No case sensitivity issues detected.\n")
  }
}

unique_locations_fish_data<- raw_fish_data$Site %>%
  unique() %>%
  sort()
unique_locations_fish_data

unique_locations_location<- raw_location$Site %>%
  unique() %>%
  sort()
unique_locations_location

unique_locations_2018 <- fish_abun_18_23_raw_location$Location %>%
  unique() %>%
  sort()

unique_locations_2009
unique_locations_2018

fish_abun_09_16_raw_LS <- read_excel(here("raw_data", "DATA-fish abund-2009-2016_20Dec2022.xlsx"), sheet = 2) #Make sure to specify the specific sheet!!
fish_abun_09_16_raw_TS <- read_excel(here("raw_data", "DATA-fish abund-2009-2016_20Dec2022.xlsx"), sheet = 3)
fish_abun_09_16_raw_location <- read_excel(here("raw_data", "DATA-fish abund-2009-2016_20Dec2022.xlsx"), sheet = 4)
fish_abun_09_16_raw_Dive_details <- read_excel(here("raw_data", "DATA-fish abund-2009-2016_20Dec2022.xlsx"), sheet = 5) %>% 
  select(-c("Dive time", "Other notes...11", "Other notes...12")) # Removing the lsited items - do not need these 
fish_abun_09_16_raw_Location_descr <- read_excel(here("raw_data", "DATA-fish abund-2009-2016_20Dec2022.xlsx"), sheet = 6)
fish_abun_18_23_raw_TS <- read_excel(here("raw_data", "CORDIO Fish_abund_2018_2023 for Laura_19122023.xlsx"), sheet = 1)
fish_abun_18_23_raw_location <- read_excel(here("raw_data", "CORDIO Fish_abund_2018_2023 for Laura_19122023.xlsx"), sheet = 2)
fish_abun_18_23_raw_Dive_details <- read_excel(here("raw_data", "CORDIO Fish_abund_2018_2023 for Laura_19122023.xlsx"), sheet = 3)
chl_09_15_raw <- read_excel(here("raw_data", "MS_Fish_sites_Master_Chl_09_15_200623.xlsx"), sheet = 3)
sst_09_15_raw <- read_excel(here("raw_data", "MS_Fish_sites_Master_sst_09_15_140623.xlsx"))


#Exploring the data pipe
# data %>%
  #view() %>% # Open data in a new tab
  #dim() %>% # Explore dimensions of data
  #glimpse() %>% # Explore structure of data
  #head() %>% # Explore first 6 rows
  #tail() # Explore last 6 rows








# A function to make sure the column names are cleaned up, 
# eg lower case and snake case
clean_column_names <- function(penguins_data) {
  penguins_data %>%
    clean_names()
}

# A function to remove columns based on a vector of column names
remove_columns <- function(penguins_data, column_names) {
  penguins_data %>%
    select(-starts_with(column_names))
}

# A function to make sure the species names are shortened
shorten_species <- function(penguins_data) {
  penguins_data %>%
    mutate(species = case_when(
      species == "Adelie Penguin (Pygoscelis adeliae)" ~ "Adelie",
      species == "Chinstrap penguin (Pygoscelis antarctica)" ~ "Chinstrap",
      species == "Gentoo penguin (Pygoscelis papua)" ~ "Gentoo"
    ))
}

# A function to remove any empty columns or rows
remove_empty_columns_rows <- function(penguins_data) {
  penguins_data %>%
    remove_empty(c("rows", "cols"))
}


# A function to remove rows which contain NA values
remove_NA <- function(penguins_data) {
  penguins_data %>%
    na.omit()
}




