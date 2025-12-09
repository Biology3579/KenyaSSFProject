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

# Extracting chl_a data ----

## Libraries ----
library(tidyverse)   # 
library(data.table)  # fast table operations (fread/fwrite, rbindlist)
library(hdf5r)       # read HDF5-style files (NASA L3b)
library(arrow)       # write_parquet
library(sf)          # spatial matching utilities (st_nearest_feature)
library(dplyr)       # light data manipulation helpers (used in a few helpers)

## Processing raw data ----

### A function to parse the date from filename ----
# Extracts first 8-digit sequence (YYYYMMDD) from the HDF5-NetCDF file and returns Date or NA
parse_date <- function(filename) {
  date_str <- stringr::str_extract(filename, "\\d{8}") #  searches the filename and pulls out the date
  if (is.na(date_str)) return(NA) # returns NA if it cannot find the date
  as.Date(date_str, format = "%Y%m%d") #reads the 8 digits as a date
}

### A function to NASA bin numbers to lat/lon (vectorized) ----
# Every MODIS L3b bin has a fixed location on the globe, 
# Each row in BinList corresponds to one bin where data exist and the value 
# tells you which bin the chlorophyll data refers to. 
# To be able to join this data to my sites, I need to extract the coordinates. 
# This uses the bin numbers - stored in bin_num  integer or numeric vector (0-based); numrows default set to 8640
nasa_bin_to_latlon <- function(bin_nums, 
                               numrows = 8640) { # 8640 rows corresponds to resolution used by MODIS L3b global bins
 
  bin_nums <- as.numeric(bin_nums) # Convert vectors to numeric to avoid maths errors (in further equations)
  is_invalid <- is.na(bin_nums) | bin_nums < 0 #deals with potential processing errors
  
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
  # - Divide by total bins → fraction of full 360° circle                              numrows = 8640) {                                numrows = 8640) {                                numrows = 8640) { 
  longitude <- 360.0 * (column_in_row + 0.5) / bins_in_this_row  - 180.0  

  # Set invalid bins to NA
  latitude[is_invalid] <- NA_real_
  longitude[is_invalid] <- NA_real_
  
  list(lat = latitude, lon = longitude) #returns coordinates
}


## A function to read and process NASA L3b files ----
# Reads each NASA L3b bin files (HDF5/NetCDF style HDF5 layout).
# and extracts relevant values to compute chlor_mean, 
# convert bin numbers to lat/lon, 
# and returns combined data.table. 

### Process single file 
process_single_file <- function(filepath, AOI = NULL) {
 
  # open file
    h5file <- hdf5r::H5File$new(filepath, mode = "r")
    
  # reads the relevant subfolders from the hdf5 files
    # subfolder containing the chlor_a values
    chlor_raw <- h5file[["level-3_binned_data/chlor_a"]]$read()
    # subfolder containing the bin_names 
    binlist_raw <- h5file[["level-3_binned_data/BinList"]]$read()
    
  # converts subfolders to data.frame for ease of use
    chlor_df  <- as.data.frame(chlor_raw, stringsAsFactors = FALSE)
    binlist_df <- as.data.frame(binlist_raw, stringsAsFactors = FALSE)
    
  # calculate chlor_a mean
  # - this will be chlor_a sum / nobs 
  
    # 1. Extract the values as numbers (for ease of use in the equation)
    chlor_sum <- as.numeric(chlor_df[["sum"]])
    nobs <- as.numeric(binlist_df[["nobs"]])
    
    # 2. Compute the mean
    chlor_a_mean <- rep(NA_real_, length(nobs))
    valid <- which(!is.na(nobs) & nobs > 0 & !is.na(chlor_sum))
    if (length(valid)) chlor_a_mean[valid] <- chlor_sum[valid] / nobs[valid]

    
  # convert NASA bins to lat/lon
    # 1. extract the bin numbers
    bin_nums <- as.numeric(binlist_df[["bin_num"]])
    # 2. convert numbers to coordinates
    coords <- nasa_bin_to_latlon(bin_nums)
 
  # build data.table with relevant variables   
    dt <- data.table::data.table(
      # add extracted variables
      chlor_a_mean = chlor_a_mean,
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

## Combine with master location

# Utility: match satellite bins to station points by year-month ----
# sta_sf: sf of stations (with date column and geometry)
# sat_sf: sf of satellite bins (with date column and geometry)
# returns integer vector of indices into sat_sf for each station row
match_sat_to_stations_by_month <- function(sta_sf, sat_sf, station_date_col = "date", sat_date_col = "date") {
  if (!inherits(sta_sf, "sf") || !inherits(sat_sf, "sf")) stop("sta_sf and sat_sf must be sf objects")
  # create year-month keys
  sta_ym <- format(as.Date(sta_sf[[station_date_col]]), "%Y-%m")
  sat_ym <- format(as.Date(sat_sf[[sat_date_col]]), "%Y-%m")
  
  nearest_matches <- integer(nrow(sta_sf))
  ym_list <- unique(sta_ym)
  
  for (ym_i in ym_list) {
    sta_idx <- which(sta_ym == ym_i)
    sat_idx <- which(sat_ym == ym_i)
    
    if (length(sat_idx) == 0) {
      # fallback: global nearest
      nearest_matches[sta_idx] <- sf::st_nearest_feature(sta_sf[sta_idx, ], sat_sf)
    } else {
      rel_idx <- sf::st_nearest_feature(sta_sf[sta_idx, ], sat_sf[sat_idx, ])
      # rel_idx gives positions relative to sat_idx; map to global index
      nearest_matches[sta_idx] <- sat_idx[rel_idx]
    }
  }
  
  nearest_matches
}

# Quick QA: count missing grav_tot or other field in a table (invisibly returns count) ----
# dt: data.frame / data.table / sf (if sf, geometry dropped)
count_missing_field <- function(dt, field_name = "grav_tot", quiet = FALSE) {
  df <- if (inherits(dt, "sf")) sf::st_drop_geometry(dt) else as.data.frame(dt)
  if (!field_name %in% names(df)) {
    if (!quiet) message("Field '", field_name, "' not present")
    return(invisible(NA_integer_))
  }
  n_missing <- sum(is.na(df[[field_name]]))
  if (!quiet) {
    if (n_missing > 0) {
      message("⚠ Warning: ", n_missing, " rows missing ", field_name)
    } else {
      message("All rows have a non-missing ", field_name)
    }
  }
  invisible(n_missing)
}

# End of extracting.r
