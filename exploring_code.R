





library(dplyr)
library(stringr)
library(here)

renv::

#Load data

clean_fish_2009 <- read.csv(here("processed_data","clean_fish_2009.csv"))
clean_fish_2018 <- read.csv(here("processed_data","clean_fish_2018.csv"))
fish_groups <-  read.csv(here("raw_data","fish_groupings_samoilys_2018.csv"))

# Check unique standardized Species names (using Species_new)
unique_species_2009 <- clean_fish_2009 %>%
  filter(country == "kenya") %>% 
  mutate(
    standard_species = str_to_lower(str_trim(species_new))
  ) %>%
  distinct(standard_species) %>%
  arrange(standard_species)

# Check unique standardized Family names
unique_family_2009 <- clean_fish_2009 %>%
  filter(country == "kenya") %>% 
  mutate(
    standard_family = str_to_lower(str_trim(family))
  ) %>%
  distinct(standard_family) %>%
  arrange(standard_family)

# View the counts (you will need to scroll through this output)
print(unique_species_2009)
print(unique_family_2009)




# Check unique standardized Species names
unique_species_2018 <- clean_fish_2018 %>%
  filter(country == "kenya") %>% 
  mutate(
    standard_species = str_to_lower(str_trim(species))
  ) %>%
  distinct(standard_species) %>%
  arrange(standard_species)

# Check unique standardized Family names
unique_family_2018 <- clean_fish_2018 %>%
  filter(country == "kenya") %>% 
  mutate(
    standard_family = str_to_lower(str_trim(family))
  ) %>%
  distinct(standard_family) %>%
  arrange(standard_family)

# View the counts (you will need to scroll through this output)
print(unique_species_2018)
print(unique_family_2018)



# Which families in the 2009 data are MISSING from the groupings table?
fish_abund_2009_raw %>%
  mutate(standard_family = str_to_lower(str_trim(Family))) %>%
  anti_join(
    fish_groupings %>% mutate(standard_family = str_to_lower(str_trim(Family))),
    by = "standard_family"
  ) %>%
  distinct(standard_family)

# Which families in the 2018 data are MISSING from the groupings table?
fish_abund_2018_raw %>%
  mutate(standard_family = str_to_lower(str_trim(Family))) %>%
  anti_join(
    fish_groupings %>% mutate(standard_family = str_to_lower(str_trim(Family))),
    by = "standard_family"
  ) %>%
  distinct(standard_family)

# Species found in BOTH 2009 and 2018
common_species <- intersect(unique_species_2009, unique_species_2018)
common_families <- intersect(unique_family_2009, unique_family_2018)

# Species ONLY in 2009
unique_to_2009 <- setdiff(unique_species_2009, unique_species_2018)

# Species ONLY in 2018
unique_to_2018 <- setdiff(unique_species_2018, unique_species_2009)

# Total unique species across both datasets
total_species <- union(unique_species_2009, unique_species_2018)



length(total_species)
length(common_species)
length(common_families)
length(unique_to_2009)
length(unique_to_2018)
