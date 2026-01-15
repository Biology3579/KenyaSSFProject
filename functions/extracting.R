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
library(hdf5r)       # read HDF5-style files (NASA L3b)
library(arrow)       # write_parquet
library(sf)          # spatial matching utilities (st_nearest_feature)
library(dplyr)       # light data manipulation helpers (used in a few helpers)

## Processing extarcted data files ----

### A function to parse the date from filename ----
# Extracts first 8-digit sequence (YYYYMMDD) from the HDF5-NetCDF file and returns it as a date

parse_date <- function(filename) {
  date_str <- stringr::str_extract(filename, "\\d{8}")
  if (is.na(date_str)) return(NA)
  as.Date(date_str, format = "%Y%m%d")
}


### A function to NASA bin numbers to lat/lon (vectorized) ----
# Every MODIS L3b bin has a fixed location on the globe, 
# Each row in BinList corresponds to one bin where data exist and the value 
# tells you which bin the chlorophyll data refers to. 
# To be able to join this data to my sites, I need to extract the coordinates. 
# This uses the bin numbers - stored in bin_num  integer or numeric vector (0-based); numrows default set to 8640

nasa_bin_to_latlon <- function(bin_nums, numrows = 8640) { 
  bin_nums <- as.numeric(bin_nums)
  is_invalid <- is.na(bin_nums) | bin_nums < 0
  
# Calculate latitude
  #  First, compute the "row number" for each bin.
  #  In NASA's triangular grid indexing:
  #      - Row 0 starts at the North Pole
  #      - Row numbers increase toward the equator, then southward
  #    The formula inverts the triangular number sequence:
  #         T(n) = n(n+1)/2
  #    Solving for n gives:
  row_number <- floor((sqrt(1 + 8 * bin_nums) - 1) / 2)
  
  #  Then, convert row number → latitude (in degrees)
  #    A sinusoidal grid has evenly spaced latitude bands.
  #    - Subtract from 90° because bin indexing begins at the North Pole
  #    - Place the coordinate at the center of the row 
  #    - Adjusts for the latitude thickness of each row: Δϕ = 180/numrows
  latitude <- 90.0 - (row_number + 0.5) * 180.0 / numrows 
  
# Calculate longitude
  # First identify the first bin index of this row.
  #    Because rows have triangular-length indexing (1, 3, 5, … bins),
  #    the total number of bins before row r is:
  #       T(r-1) = r(r-1)/2
  first_bin_in_row <- row_number * (row_number - 1) / 2
  
  # Then, compute the bin’s column position within its row.
  #    This is its horizontal position in the sinusoidal grid.
  column_in_row <- bin_nums - first_bin_in_row
  
  # Then, determine number of bins in this row.
  #    In the sinusoidal grid, row r has:
  #         bins = 2 * (r + 1) - i.e., many bins near equator, fewer near poles.
  bins_in_this_row <- 2 * (row_number + 1)
  
  # Then, convert column position → longitude
  # - Multiply by 360 → degrees in [0, 360)
  # - Center of the bin horizontally
  # - Divide by total bins → fraction of full 360° circle
                               
  longitude <- 360.0 * (column_in_row + 0.5) / bins_in_this_row  - 180.0  
  latitude[is_invalid] <- NA_real_
  longitude[is_invalid] <- NA_real_
  list(lat = latitude, lon = longitude) #returns coordinates
}

## A function to read NASA L3b files and extract chla_means----

extract_raw_data <- function(filepath, AOI) {
  # open file
    h5file <- hdf5r::H5File$new(filepath, mode = "r")
    
  # reads the relevant subfolders from the hdf5 files
    # subfolder containing the chla values
    chla_raw <- h5file[["level-3_binned_data/chlor_a"]]$read()
    # subfolder containing the bin_names 
    binlist_raw <- h5file[["level-3_binned_data/BinList"]]$read()
    
  # converts subfolders to data.frame for ease of use
    chla_df  <- as.data.frame(chla_raw, stringsAsFactors = FALSE)
    binlist_df <- as.data.frame(binlist_raw, stringsAsFactors = FALSE)
    
  # calculate chla mean
  # - this will be chla sum / nobs 

    # Extract the values as numbers (for ease of use in the equation)
    chla_sum <- as.numeric(chla_df[["sum"]])
    nobs <- as.numeric(binlist_df[["nobs"]])
    
    # Compute the mean
    chla_mean <- rep(NA_real_, length(nobs))
    valid <- which(!is.na(nobs) & nobs > 0 & !is.na(chla_sum)) # checking that the lists are not empty and that nobs > 0
    if (length(valid)) chla_mean[valid] <- chla_sum[valid] / nobs[valid]

    
  # Convert NASA bins to lat/lon
    bin_nums <- as.numeric(binlist_df[["bin_num"]])
    coords <- nasa_bin_to_latlon(bin_nums)
    
  # Build data.table with relevant variables   
    dt <- data.table::data.table(
      # add extracted variables
      chla_mean = chla_mean,
      nobs = nobs,
      bin_num = bin_nums,
      lat = coords$lat,
      lon = coords$lon
    )[
      # filter for region of interest
      lat >= coords_list$S & lat <= coords_list$N &
        lon >= coords_list$W & lon <= coords_list$E
    ][
      # add metadata
      , `:=`(
        date = parse_date(basename(filepath)),
        source_file = basename(filepath)
      )
    ]
}

## A function to extract lat and lon from the OCNET satellite files ----
  # need to revisit this!! ...
extract_latlon_ocnet <- function(filepath, AOI = NULL,
                                 lon_name = "longitude",
                                 lat_name = "latitude") {
  nc <- nc_open(filepath)
  on.exit(nc_close(nc))
  
  lon <- ncvar_get(nc, lon_name)
  lat <- ncvar_get(nc, lat_name)
  
  if (!is.null(AOI)) {
    lon <- lon[lon >= AOI$W & lon <= AOI$E]
    lat <- lat[lat >= AOI$S & lat <= AOI$N]
  }
  
  data.table(lon = lon, lat = lat)
}


## A function that assigns to each station, the nearest satellite cell (store bin_num and distance to cell) ----

assign_nearest_cell <- function(raw_sat_data, clean_locations){
  
  unique_cells <- raw_sat_data %>% distinct(bin_num, lat, lon)
  
  # convert to sf (assumes columns as specified)
  stations_as_sf  <- st_as_sf(clean_locations, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
  unique_cells_as_sf  <- st_as_sf(unique_cells, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  
  # find index of nearest feature in first_sf for each point in clean_sf
  nearest_indices <- st_nearest_feature(stations_as_sf, unique_cells_as_sf)
  # nearest_idx[i] is the row index in first_sf closest to clean_sf[i,]
  
  # extract corresponding bin id
  stations_as_sf$closest_cell <- unique_cells$bin_num[nearest_indices]
  
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
    closest_cell = clean_locations$closest_cell[i]
    
    # start = first day of the month 11 months before the reference month
    start_month <- floor_date(date, "month") %m-% months(11)
    
    # end = last day of the reference month
    end_month <- ceiling_date(date, "month") - days(1)
    
    relevant_cells <- raw_sat_data %>%
      filter(
        bin_num == closest_cell,
        date >= start_month,
        date <= end_month
      ) 
    
    #find their chlorophyll average
    weighted_mean <- sum(relevant_cells[[metric_name]] * relevant_cells$nobs, na.rm = TRUE) /
      sum(relevant_cells$nobs, na.rm = TRUE)
    
    
    #add it to the table
    clean_locations[[metric_annual_mean]][i] = weighted_mean
  }
  
  clean_locations
}


# End of extracting.r
