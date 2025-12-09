#finding matched sites
#draw a circle encompassing the stations that make up the site 


# libraries
library(tidyverse)
library(dplyr)

#load the data 
locations_2009 <- read.csv(here::here("processed_data", "clean_location_2009.csv"))
locations_2018 <- read.csv(here::here("processed_data", "clean_location_2018.csv"))

# 2009 ----

## filter countries in 2009 to Kenya, Comoros and Tanzania
filtered_2009 <- locations_2009 %>% 
  filter(country =="kenya"  |country == "comoros" | country =="tanzania" ) 

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

# compute site areas
#  this draws a circle around the mean of multiple stations and finds the furthest station to make the radius for the circle
site_areas_2009 <- filtered_2009 %>% 
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
