# Map Figures
# 
# Packages
library(sf)
library(ggplot2)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)

# Load location data 
location_2009 <- readr::read_rds(here::here("processed_data", "clean_location_connectivity.rds"))
names(location_2009)
unique(location_2009$mpa_status)

# Make sites into an sf object
locations_sf <- st_as_sf(
  location_2009,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)



# Load Africa Boundaries
world <- ne_countries(scale = "medium", returnclass = "sf")
africa <- world %>%
  filter(continent == "Africa")

# Define area of interest (AOI)
wio_box <- st_as_sfc(
  st_bbox(
    c(
      xmin = 37,
      xmax = 52,
      ymin = -16.5,
      ymax = -4
    ),
    crs = st_crs(4326)
  )
)

# Plot the whole of africa
# Plot Africa with the WIO box
ggplot() +
  # Africa layer
  geom_sf(data = africa,
          fill = "grey85",
          colour = "grey60",
          linewidth = 0.01) +
  # WIO box 
  geom_sf(data = wio_box,
          fill = NA,
          colour = "black",
          linewidth = 0.8) +
  # Sites layer
  geom_sf(data = locations_sf,
          aes(shape = mpa_status, fill = mpa_status),
          colour = "black",
          size = 2.5,
          stroke = 0.4) +
  scale_shape_manual(values = c(
    "none" = 21,
    "low" = 22,
    "medium" = 24
  )) +
  scale_fill_manual(values = c(
    "none" = "white",
    "low" = "grey60",
    "medium" = "black"
  )) +
  theme_void()


# Crop Africa map to the region of interest
africa_crop <- st_crop(africa, st_bbox(wio_box))

#Plot the cropped region
ggplot() +
  
  geom_sf(data = africa_crop,
          fill = "grey85",
          colour = "grey60",
          linewidth = 0.2) +
  
  geom_sf(data = locations_sf,
          aes(shape = mpa_status,
              fill = mpa_status),
          colour = "black",
          size = 3,
          stroke = 0.4) +
  
  scale_shape_manual(values = c(
    "none" = 21,
    "low" = 22,
    "medium" = 24
  )) +
  
  scale_fill_manual(values = c(
    "none" = "white",
    "low" = "grey60",
    "medium" = "black"
  )) +
  
  theme_void()

