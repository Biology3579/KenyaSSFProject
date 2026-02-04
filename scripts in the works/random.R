# Testing for eco-region 


#Load require libraries
library(tidyverse)
library(arm)
library(here)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggfortify) # extends ggplot2 for some more ease with stats
library(MASS) # fits regression models etc.






# 2009-2015 data

## Load data for 2009-2015


#| label: load-data

fish_data_2009 <- readr::read_rds(here::here("processed_data", "clean_fish_2009.rds"))
locations_2009 <- readr::read_rds(here::here("processed_data", "clean_location_2009.rds"))
sst_2009 <- read.csv(here::here("processed_data", "locations_with_sst_2009.csv"))
chla_2009 <- read.csv(here::here("processed_data", "locations_with_chla_2009.csv"))
h_gravity_2009 <- read.csv(here::here("processed_data", "locations_with_gravity.csv"))




library(sf)
library(ggplot2)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)


world <- ne_countries(scale = "medium", returnclass = "sf")
focus_region <- world %>%
  filter(admin %in% c("Kenya", "Tanzania", "Comoros"))

# Create the map with a tighter zoom
map_plot_2009 <- ggplot() +
  
  # Plot all countries in Africa as the base layer
  geom_sf(data = world[world$continent == "Africa",], 
          fill = "gray95", color = "gray50", size = 0.3) + 
  
  geom_sf(data = focus_region, fill = "lightblue", color = "gray30", size = 0.3) +
  
  geom_point(data = locations_2009,
             aes(x = longitude, y = latitude, ),
             alpha = 0.6, color = "#440154") +
  
  scale_size_continuous(name = "Number of\nObservations", 
                        range = c(0.5, 5),
                        breaks = c(10, 25, 50, 75, 100, 150, 200)) +
  
  # Specify region of focus - Kenya in frame
  coord_sf(xlim = c(36, 49), ylim = c(1, -15)) + 
  
  labs(title = "Marine Survey Sites in Kenya",
       subtitle = "Point size represents number of fish observations per site",
       x = "Longitude", y = "Latitude") +
  
  theme_minimal() +
  theme(panel.grid.major = element_line(color = "gray90", size = 0.2),
        panel.background = element_rect(fill = "aliceblue"),
        plot.title = element_text(face = "bold", size = 14),
        legend.position = "right")

print(map_plot_2009)


## Total Biomass

This refers to an aggregate of all species

### Calculate total biomass per station, site, and country


# Mean biomass per station
biomass_per_station <- fish_data_2009 %>%
  group_by(station) %>%
  summarise(
    mean_tot_biomass = mean(tot_wt_g, na.rm = TRUE),
    site = first(site),
    country = first(country), 
    .groups = "drop")

# Mean biomass per site
biomass_per_site <- biomass_per_station %>%
  group_by(site) %>%
  summarise(
    mean_biomass = mean(mean_tot_biomass),
    country = first(country))

# Mean biomass per country
biomass_per_country <- biomass_per_site %>%
  group_by(country) %>%
  summarise(
    mean_biomass = mean(mean_biomass, na.rm = TRUE))


### Calculate env. variable means


# Calculate mean chla per site
chla_2009_sites <- chla_2009 %>% 
  group_by(site) %>% 
  summarise(
    mean_annual_chla = mean(chla_annual_mean))

# Calculate mean SST per site
sst_2009_sites <- sst_2009 %>% 
  group_by(site) %>% 
  summarise(
    mean_annual_sst = mean(sst_annual_mean))

# Calculate mean human gravity per site
h_gravity_2009_sites <- h_gravity_2009 %>% 
  group_by(site) %>% 
  summarise(
    mean_gravity = mean(Grav_tot))

### Linear regressions

#### Biomass vs Country (Ecoregion)

# Graph the data - Barchart
biomass_country_fig1 <- biomass_per_country %>% 
  ggplot(aes(x=country, y=mean_biomass))+
  geom_bar(stat = "identity")

biomass_country_fig1

# Build lm 
lm_biomass_country <- lm(mean_biomass~country,data = biomass_per_country)

# Analyse the model 
anova(lm_biomass_country)
summary(lm_biomass_country)