# ---------------------------
#
# Script name: extracting.r
#
# Purpose of script:
#   Small, reusable functions to parse satellite bin files (NASA L3b),
#   convert bin numbers to lat/lon, process single files, and combine results.
#
# Author: Candela Ferrer Diez
#
# Date Created: 2025-12-01
#
# Notes:
#   - This file contains only function definitions and library() calls.
#   - No top-level processing or example data to avoid side-effects when sourced.
#
# ---------------------------

# Extracting chl_a and SST data ----

## Libraries ----
library(tidyverse)   # 
library(data.table)  # fast table operations (fread/fwrite, rbindlist)
library(arrow)       # write_parquet
library(sf)          # spatial matching utilities (st_nearest_feature)
library(dplyr)       # light data manipulation helpers (used in a few helpers)
library(ncdf4)

## Processing extarcted data files ----

### A function to parse the date from filename ----
# Extracts first 8-digit sequence (YYYYMMDD) from the NetCDF file and returns it as a date

parse_date <- function(filename) {
  date_str <- stringr::str_extract(filename, "\\d{8}")
  if (is.na(date_str)) return(NA)
  as.Date(date_str, format = "%Y%m%d")
}


## A function to extract metric (e.g chla or SST) and coordinate information from extracted satellite files ----
extract_sat_data <- function(filepath,
                                 lon_name = "lon",
                                 lat_name = "lat",
                                 metric_name = "chlor_a",
                                 filter = FALSE,
                                 keep_idxs = NA) {
  nc <- nc_open(filepath)
  on.exit(nc_close(nc))
  
  lon  <- ncvar_get(nc, lon_name)
  lat  <- ncvar_get(nc, lat_name)
  metric <- ncvar_get(nc, metric_name)

  #Get all indices (we are not removing NAs as we don't want to change the shape)
  idx <- which(is.na(metric)|!is.na(metric), arr.ind = TRUE)
  if (filter){
    if (all(is.na(keep_idxs))){ #If we havent provided filter indices, filter by AOI
      lon_max = 55
      lon_min = 37
      lat_min = -19
      lat_max = 1
      
      # Extract coordinates for the selected indices
      lons <- lon[idx[, 1]]
      lats <- lat[idx[, 2]]
      
      # Logical mask for bounding box
      keep <- lons >= lon_min & lons <= lon_max &
        lats >= lat_min & lats <= lat_max
      

      # Filter indices
      idx <- idx[keep, ,drop = FALSE]
      
    }
    else{
      idx <- keep_idxs
    }
    
  }



  # map indices to coordinates
  dt <- data.table(
    lon = lon[idx[, 1]],  # first column in index → longitude
    lat = lat[idx[, 2]],  # second column in index → latitude
    value = metric[idx]
  )[
    # add metadata
    , `:=`(
      date = parse_date(basename(filepath)),
      source_file = basename(filepath)
    )
  ]
  setnames(dt, "value", metric_name)
  message(filepath)
  dt
}


## A function that assigns to each station, the nearest satellite cell (store lat/lon and distance to cell) ----
assign_nearest_cell <- function(raw_sat_data, clean_locations){
  
  unique_cells <- raw_sat_data %>% distinct(lat, lon)
  
  # convert to sf
  stations_as_sf     <- st_as_sf(clean_locations, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
  unique_cells_as_sf <- st_as_sf(unique_cells, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  
  # find nearest cell for each station
  nearest_indices <- st_nearest_feature(stations_as_sf, unique_cells_as_sf)
  
  clean_locations$closest_cell_lat <- unique_cells$lat[nearest_indices]
  clean_locations$closest_cell_lon <- unique_cells$lon[nearest_indices]
  
  dists <- st_distance(stations_as_sf, unique_cells_as_sf[nearest_indices,], by_element = TRUE)
  clean_locations$cell_distance_km <- as.numeric(dists) / 1000
  
  clean_locations
}


## A function to compute the annual mean ----

compute_annual_mean <- function(clean_locations, raw_sat_data, metric_name, metric_annual_mean){
  for (i in 1:nrow(clean_locations)){
    
    ref_date         <- clean_locations$date[i]
    closest_cell_lat <- clean_locations$closest_cell_lat[i]
    closest_cell_lon <- clean_locations$closest_cell_lon[i]
    
    start_month <- floor_date(ref_date, "month") %m-% months(11)
    end_month   <- ceiling_date(ref_date, "month") - days(1)
    
    relevant_cells <- raw_sat_data %>%
      filter(
        lat == closest_cell_lat,
        lon == closest_cell_lon,
        .data$date >= start_month,
        .data$date <= end_month
      )
    
    # if no data in window, find nearest cell that has data
    if (nrow(relevant_cells) == 0) {
      
      # get all cells that DO have data in the window
      cells_with_data <- raw_sat_data %>%
        filter(.data$date >= start_month, .data$date <= end_month) %>%
        distinct(lat, lon)
      
      if (nrow(cells_with_data) > 0) {
        # find which of those is geographically nearest to this station
        station_sf    <- st_as_sf(clean_locations[i, ], coords = c("longitude", "latitude"), crs = 4326)
        fallback_sf   <- st_as_sf(cells_with_data, coords = c("lon", "lat"), crs = 4326)
        nearest_idx   <- st_nearest_feature(station_sf, fallback_sf)
        
        fallback_lat  <- cells_with_data$lat[nearest_idx]
        fallback_lon  <- cells_with_data$lon[nearest_idx]
        
        relevant_cells <- raw_sat_data %>%
          filter(
            lat == fallback_lat,
            lon == fallback_lon,
            .data$date >= start_month,
            .data$date <= end_month
          )
        
        message(sprintf(
          "Station %s: no data at closest cell, fell back to nearest cell with data (%.4f, %.4f)",
          clean_locations$station[i], fallback_lat, fallback_lon
        ))
      }
    }
    
    clean_locations[[metric_annual_mean]][i] <- mean(relevant_cells[[metric_name]], na.rm = TRUE)
  }
  
  clean_locations
}


# End of extracting.r
