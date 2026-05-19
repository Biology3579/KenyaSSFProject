library(sf)
library(ggplot2)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)
library(cowplot)
library(grid)

# Load location data
location_2009 <- readr::read_rds(here::here("processed_data", "clean_location_connectivity.rds"))

locations_sf <- st_as_sf(
  location_2009,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)

# Africa boundaries
world <- ne_countries(scale = "medium", returnclass = "sf")
africa <- world %>%
  filter(continent == "Africa")

# Main WIO box
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

# Projection for the WIO panel
wio_crs <- "+proj=aea +lat_1=-6 +lat_2=-20 +lat_0=-12 +lon_0=44.5 +datum=WGS84 +units=m +no_defs"

# Transform layers
africa_wio <- st_transform(africa, wio_crs)
sites_wio   <- st_transform(locations_sf, wio_crs)
box_wio     <- st_transform(wio_box, wio_crs)
box_lim     <- st_bbox(box_wio)

# Zoom box for southern Tanzania / northern Mozambique
tanzania_moz_box_ll <- st_as_sfc(
  st_bbox(
    c(
      xmin = 38.8,
      xmax = 42.2,
      ymin = -12.8,
      ymax = -6.0
    ),
    crs = st_crs(4326)
  )
)

tanzania_moz_box <- st_transform(tanzania_moz_box_ll, wio_crs)
tanmoz_lim <- st_bbox(tanzania_moz_box)

# Main regional map
main_map <- ggplot() +
  geom_sf(data = africa_wio,
          fill = "grey85",
          colour = "grey70",
          linewidth = 0.2) +
  geom_sf(data = sites_wio,
          shape = 21,
          fill = "black",
          colour = "white",
          size = 2.8,
          stroke = 0.4) +
  geom_sf(data = tanzania_moz_box,
          fill = NA,
          colour = "grey40",
          linewidth = 0.5,
          linetype = "dashed") +
  coord_sf(
    xlim = c(box_lim["xmin"], box_lim["xmax"]),
    ylim = c(box_lim["ymin"], box_lim["ymax"]),
    expand = FALSE
  ) +
  annotate("text", x = -410000, y = 220000,  label = "Tanzania",   size = 3, colour = "grey30") +
  annotate("text", x = -480000, y = -180000, label = "Mozambique", size = 3, colour = "grey30") +
  annotate("text", x = 200000,  y = 100000,   label = "Comoros",    size = 3, colour = "grey30") +
  annotate("text", x = 600000,  y = -200000,  label = "Madagascar", size = 3, colour = "grey30") +
  theme_void()

# Africa locator inset
africa_map <- ggplot() +
  geom_sf(data = africa_wio, fill = "grey85", colour = "grey70", linewidth = 0.15) +
  geom_sf(data = box_wio, fill = NA, colour = "black", linewidth = 0.8) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.3))

# Tight inset for southern Tanzania / northern Mozambique
tanzania_moz_inset <- ggplot() +
  geom_sf(data = africa_wio,
          fill = "grey85",
          colour = "grey70",
          linewidth = 0.2) +
  geom_sf(data = sites_wio,
          shape = 21,
          fill = "black",
          colour = "white",
          size = 3,
          stroke = 0.4) +
  coord_sf(
    xlim = c(tanmoz_lim["xmin"], tanmoz_lim["xmax"]),
    ylim = c(tanmoz_lim["ymin"], tanmoz_lim["ymax"]),
    expand = FALSE
  ) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.3),
    plot.margin = margin(2, 2, 2, 2)
  )

# Final composition
final_map <- ggdraw() +
  draw_plot(main_map, x = 0, y = 0, width = 1, height = 1) +
  draw_plot(tanzania_moz_inset, x = 0.40, y = 0.46, width = 0.34, height = 0.36) +
  draw_plot(africa_map, x = 0.79, y = 0.05, width = 0.18, height = 0.22) +
  draw_label("A", x = 0.42, y = 0.96, fontface = "bold", size = 12)

final_map
