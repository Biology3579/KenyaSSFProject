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
##   The raw data involves sites all over the WIO. 
##   In this analysis I focus on data from Kenya only, so I will clean and process this data.
##
## ---------------------------

#Load Libraries
library(tidyverse) #for ...
library(readxl) # To read an excel file (within tidyverse)
library(here) # Reproducible working directories

#Load separate data files

# Locations

fish_abun_09_16_raw_Master_Location <- read_excel(
  here("raw_data", "DATA-fish abund-2009-2016_20Dec2022.xlsx"), # Original data file
  sheet = 4) %>%  # Excel sheet number (as listed in original file)
  filter(Country == "Kenya") # Filter for Kenya only 
fish_abun_18_23_raw_Master_Location <- read_excel(
  here("raw_data", "CORDIO Fish_abund_2018_2023 for Laura_19122023.xlsx"), 
  sheet = 2) %>% 
filter(Country == "Kenya") # Filter for Kenya only 

unique_sites_2009 <- fish_abun_09_16_raw_Master_Location$Site %>%
  unique() %>%
  sort()

unique_sites_2018 <- fish_abun_18_23_raw_Master_Location$Site %>%
  unique() %>%
  sort()

unique_sites_2009
unique_sites_2018

unique_locations_2009 <- fish_abun_09_16_raw_Master_Location$Location %>%
  unique() %>%
  sort()

unique_locations_2018 <- fish_abun_18_23_raw_Master_Location$Location %>%
  unique() %>%
  sort()

unique_locations_2009
unique_locations_2018

fish_abun_09_16_raw_LS <- read_excel(here("raw_data", "DATA-fish abund-2009-2016_20Dec2022.xlsx"), sheet = 2) #Make sure to specify the specific sheet!!
fish_abun_09_16_raw_TS <- read_excel(here("raw_data", "DATA-fish abund-2009-2016_20Dec2022.xlsx"), sheet = 3)
fish_abun_09_16_raw_Master_Location <- read_excel(here("raw_data", "DATA-fish abund-2009-2016_20Dec2022.xlsx"), sheet = 4)
fish_abun_09_16_raw_Dive_details <- read_excel(here("raw_data", "DATA-fish abund-2009-2016_20Dec2022.xlsx"), sheet = 5) %>% 
  select(-c("Dive time", "Other notes...11", "Other notes...12")) # Removing the lsited items - do not need these 
fish_abun_09_16_raw_Location_descr <- read_excel(here("raw_data", "DATA-fish abund-2009-2016_20Dec2022.xlsx"), sheet = 6)
fish_abun_18_23_raw_TS <- read_excel(here("raw_data", "CORDIO Fish_abund_2018_2023 for Laura_19122023.xlsx"), sheet = 1)
fish_abun_18_23_raw_Master_Location <- read_excel(here("raw_data", "CORDIO Fish_abund_2018_2023 for Laura_19122023.xlsx"), sheet = 2)
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



