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
                                 metric_name = "chlor_a") {
  nc <- nc_open(filepath)
  on.exit(nc_close(nc))
  
  lon  <- ncvar_get(nc, lon_name)
  lat  <- ncvar_get(nc, lat_name)
  metric <- ncvar_get(nc, metric_name)
  
  # get indices of non-zero metric values
  idx <- which(metric != 0, arr.ind = TRUE)
  
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
  dt
}


## A function that assigns to each station, the nearest satellite cell (store lat/lon and distance to cell) ----
assign_nearest_cell <- function(raw_sat_data, clean_locations){
  
  unique_cells <- raw_sat_data %>% distinct(lat, lon)
  
  # convert to sf (assumes columns as specified)
  stations_as_sf  <- st_as_sf(clean_locations, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
  unique_cells_as_sf  <- st_as_sf(unique_cells, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  
  # find index of nearest feature in first_sf for each point in clean_sf
  nearest_indices <- st_nearest_feature(stations_as_sf, unique_cells_as_sf)
  # nearest_idx[i] is the row index in first_sf closest to clean_sf[i,]
  
  # extract corresponding lat/lon
  clean_locations$closest_cell_lat <- unique_cells$lat[nearest_indices]
  clean_locations$closest_cell_lon <- unique_cells$lon[nearest_indices]
  
  # compute accurate geodetic distances (in meters) for those pairs
  # st_distance on geographic coordinates uses s2/geodetic when available; specify by_element=TRUE
  dists <- st_distance(stations_as_sf, unique_cells_as_sf[nearest_indices,], by_element = TRUE)
  
  # convert to numeric (in meters)
  stations_as_sf$cell_distance_km <- as.numeric(dists) / 1000
  
  # put back into the original dataframe if needed
  clean_locations$closest_cell   <- stations_as_sf$closest_cell
  clean_locations$cell_distance_km  <- stations_as_sf$cell_distance_km
  
  clean_locations
  
}


## A function to compute the annual mean ----

compute_annual_mean <- function(clean_locations, raw_sat_data, metric_name, metric_annual_mean){
  for (i in 1:length(clean_locations$station)){
    #get month and year and id
    date <- clean_locations$date[i]
    closest_cell_lat = clean_locations$closest_cell_lat[i]
    closest_cell_lon = clean_locations$closest_cell_lon[i]
    
    # start = first day of the month 11 months before the reference month
    start_month <- floor_date(date, "month") %m-% months(11)
    
    # end = last day of the reference month
    end_month <- ceiling_date(date, "month") - days(1)

    relevant_cells <- raw_sat_data %>%
      filter(
        lat == closest_cell_lat,
        lon == closest_cell_lon,
        date >= start_month,
        date <= end_month
      ) 
    
    #find their chlorophyll average
    mean <- mean(relevant_cells[[metric_name]])

    #add it to the table
    clean_locations[[metric_annual_mean]][i] = mean
  }
  
  clean_locations
}


# End of extracting.r
