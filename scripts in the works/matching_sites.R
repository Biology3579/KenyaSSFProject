#finding matched sites
#draw a circle encompassing the stations that make up the site 


# libraries
library(tidyverse)
library(dplyr)
library(ggplot2)
library(maps)
library(sf)
library(RColorBrewer) # or use viridis for colorblind-friendly palettes

#load the data 
locations_2009 <- read_rds(here::here("processed_data", "clean_location_2009.rds"))
locations_2018 <- read_rds(here::here("processed_data", "clean_location_2018.rds"))


# Visually explore the location of the stations by country

# Get map data for Kenya and surrounding region

#Map colours 
country_cols <- c(
  "kenya"    = "#1b9e77",  # teal
  "tanzania" = "#d95f02",  # orange
  "comoros"  = "#7570b3"   # purple
)

world <- ne_countries(scale = "medium", returnclass = "sf")
focus_region <- world %>%
  filter(admin %in% c("Kenya", "Tanzania", "Comoros"))

# Create the map with a tighter zoom
map_plot<- ggplot() +
  
  # Plot all countries in Africa as the base layer
  geom_sf(data = world[world$continent == "Africa",], 
          fill = "gray95", color = "gray50", size = 0.3) + 
  
  geom_sf(data = focus_region, fill = "lightblue", color = "gray30", size = 0.3) +
  
  geom_point(data = locations_2009$station, groupby(country),
             aes(x = longitude, y = latitude, ),
             alpha = 0.6, color = "#440154") +
  # each country gets a different 
  geom_point(data = , 
             aes(x = longitude, y = latitude, size = n_observations),
             alpha = 0.6, color = "#440154") +
  
  scale_size_continuous(name = "Number of\nObservations", 
                        range = c(0.5, 5),
                        breaks = c(10, 25, 50, 75, 100, 150, 200)) +
  
  # Specify region of focus - Kenya in frame
  coord_sf(xlim = c(38.9, 41.8), ylim = c(-5, -1.5)) + 
  
  labs(title = "Marine Survey Sites in Kenya",
       subtitle = "Point size represents number of fish observations per site",
       x = "Longitude", y = "Latitude") +
  
  theme_minimal() +
  theme(panel.grid.major = element_line(color = "gray90", size = 0.2),
        panel.background = element_rect(fill = "aliceblue"),
        plot.title = element_text(face = "bold", size = 14),
        legend.position = "right")

print(map_plot_2009)


# find true sites
# this finds stations that are within 4km of each other to define the sites 
radius <- 1000  # metres

#fucntion to get overlaps
do_overlap <- function(lat1, lon1, lat2, lon2, radius){
  dist <- coord_distance(lat1, lat2, lon1, lon2)
  dist < (2 * radius)
}

# initialise station-level clustering table
station_areas_2009 <- stations_2009 %>%
  mutate(
    radius = radius,
    cluster_id = 0
  )

#cluster circles tat overlap
next_id <- 1
n <- nrow(station_areas_2009)

for (i in 1:(n - 1)) {
  for (j in (i + 1):n) {
    
    # assign a cluster to i if it doesn't have one
    if (station_areas_2009$cluster_id[i] == 0) {
      station_areas_2009$cluster_id[i] <- next_id
      current_id <- next_id
      next_id <- next_id + 1
    } else {
      current_id <- station_areas_2009$cluster_id[i]
    }
    
    lat1 <- station_areas_2009$latitude[i]
    lon1 <- station_areas_2009$longitude[i]
    lat2 <- station_areas_2009$latitude[j]
    lon2 <- station_areas_2009$longitude[j]
    
    if (do_overlap(lat1, lon1, lat2, lon2, radius)) {
      station_areas_2009$cluster_id[j] <- current_id
    }
  }
}


# find overlapping sites

#find the centre of the coordinates
centroid_pts <- function(lats, longs){
  sf::sf_use_s2(TRUE)
  
  # Build an sf POINT object from lon/lat vectors to be able to measure distance in metres
  pts <- st_as_sf(
    data.frame(lon = longs, lat = lats),
    coords = c("lon", "lat"),
    crs = 4326
  )
  # Compute centroid of the unioned set of points (geodesic)
  centroid <- st_centroid(st_union(pts))
  
  # Extract lon/lat coordinates
  coords <- st_coordinates(centroid)
  
  # Return named lat/lon vector
  coords
}

# compute distance in metres between two points (lat/lon)
coord_distance <- function(lat1, lat2, lon1, lon2) {
  # ensure s2 is on
  sf::sf_use_s2(TRUE)

  # number of entries
  n <- length(lat2)

  # repeat reference coordinates if needed
  if (length(lat1) == 1) lat1 <- rep(lat1, n)
  if (length(lon1) == 1) lon1 <- rep(lon1, n)

  # result vector
  distances <- numeric(n)

  # loop element-wise
  for (i in seq_len(n)) {

    p1 <- st_sfc(st_point(c(lon1[i], lat1[i])), crs = 4326)
    p2 <- st_sfc(st_point(c(lon2[i], lat2[i])), crs = 4326)

    # st_distance returns a units object; convert to numeric
    distances[i] <- as.numeric(st_distance(p1, p2))
  }

  distances
}

#1km radius ----

# compute site areas 
#  this draws a circle around the mean of multiple stations and finds the furthest station to make the radius for the circle
site_areas_2009 <- stations_2009 %>% 
  group_by(site) %>% 
  summarise(
    centre_latitude = centroid_pts(latitude, longitude)[2],
    centre_longitude = centroid_pts(latitude, longitude)[1],
    radius = max(500, max(coord_distance(centre_latitude, latitude, centre_longitude, longitude ))),
    country = country[1],
    geology = geology[1],
    reef_geomorphology = reef_geomorphology[1],
    reef_type = reef_type[1],
    cluster_id = 0
    ) 

# 2018 ----
# no need to filter (only contains Kenya, Tanzania and Comoros)
site_areas_2018 <- locations_2018 %>% 
  group_by(site) %>% 
  summarise(
    centre_latitude = centroid_pts(latitude, longitude)[2],
    centre_longitude = centroid_pts(latitude, longitude)[1],
    radius = max(500, max(coord_distance(centre_latitude, latitude, centre_longitude, longitude ))),
    country = country[1],
    geology = geology[1],
    reef_geomorphology = reef_geomorphology[1],
    reef_type = reef_type[1],
    cluster_id = 0
  ) 

# find overlap between sites in 2009 and 2018 
# comparing to sites in 2018 (as the reference point)
# 
# function to determine if two circles overlap
do_overlap <- function(lat1, lon1, rad1, lat2, lon2, rad2 ){
  dist = coord_distance(lat1, lat2, lon1, lon2)
  if ((dist < (rad1+rad2)) == TRUE){
    return (TRUE)
  }
  else{
    return (FALSE)
  }
  
}

# clusters stations within the 2009 dataset (to check if any overlap within 1km of each other)
next_id <- 1
for (i in 1:(length(site_areas_2009$cluster_id)-1)){
  for (j in (i+1):length(site_areas_2009$cluster_id)){
    # compare rows i and j (do they overlap)
    current_id = next_id
    if (site_areas_2009$cluster_id[i] == 0){
      site_areas_2009$cluster_id[i] = current_id
      next_id = next_id + 1
    }
    else{
      current_id = site_areas_2009$cluster_id[i]
    }
    
    lat1 = site_areas_2009$centre_latitude[i]
    lon1 = site_areas_2009$centre_longitude[i]
    rad1 = site_areas_2009$radius[i]
    lat2 = site_areas_2009$centre_latitude[j]
    lon2 = site_areas_2009$centre_longitude[j]
    rad2 = site_areas_2009$radius[j]
    
    if (do_overlap(lat1, lon1, rad1, lat2, lon2, rad2 ) == TRUE){
      site_areas_2009$cluster_id[j] = current_id
    }
  }
}



# clusters sites within the 2018 dataset (to check if any overlap within 1km of each other)
next_id <- 1
for (i in 1:(length(site_areas_2018$cluster_id)-1)){
  for (j in (i+1):length(site_areas_2018$cluster_id)){
    # compare rows i and j (do they overlap)
    current_id = next_id
    if (site_areas_2018$cluster_id[i] == 0){
      site_areas_2018$cluster_id[i] = current_id
      next_id = next_id + 1
    }
    else{
      current_id = site_areas_2018$cluster_id[i]
    }
    
    lat1 = site_areas_2018$centre_latitude[i]
    lon1 = site_areas_2018$centre_longitude[i]
    rad1 = site_areas_2018$radius[i]
    lat2 = site_areas_2018$centre_latitude[j]
    lon2 = site_areas_2018$centre_longitude[j]
    rad2 = site_areas_2018$radius[j]
    
    if (do_overlap(lat1, lon1, rad1, lat2, lon2, rad2 ) == TRUE){
      site_areas_2018$cluster_id[j] = current_id
    }
  }
}

# comparing 2009 to 2018 
for (i in 1:length(site_areas_2009$cluster_id)){
  for (j in 1:length(site_areas_2018$cluster_id)){
    # compare rows i and j (do they overlap)
    
    lat1 = site_areas_2009$centre_latitude[i]
    lon1 = site_areas_2009$centre_longitude[i]
    rad1 = site_areas_2009$radius[i]
    lat2 = site_areas_2018$centre_latitude[j]
    lon2 = site_areas_2018$centre_longitude[j]
    rad2 = site_areas_2018$radius[j]
    
    if (do_overlap(lat1, lon1, rad1, lat2, lon2, rad2 ) == TRUE){
      if ((site_areas_2009$cluster_id[i] != 0) & (site_areas_2009$cluster_id[i] != site_areas_2018$cluster_id[j])){
        message("WARNING: close to two clusters:", i)
        message(j)
      }
      site_areas_2009$cluster_id[i] = site_areas_2018$cluster_id[j]
    }
  }
}


