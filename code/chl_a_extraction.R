# Install/update required packages
install.packages(c("abind", "RNetCDF"))

# Install ncdfCF from GitHub (required for read_l5b to work)
# If you don't have it, install devtools first: install.packages("devtools")
# devtools::install_github("R-CF/ncdfCF") 
library(ncdfCF)
library(abind)
library(RNetCDF)


ncdf4::nc_open(fn)

# --- 3. Loop, Extract, and Bind the Data ---

# Define the file path using the here package
file_path <- here::here("raw_data", "requested_files") 

# List all NetCDF files (.nc) in that directory
# full.names = TRUE ensures you get the complete path for each file
file_list <- list.files(file_path, pattern = "\\.nc$", full.names = TRUE)

# Double-check that you found the correct number of files (3,422)
print(paste("Found", length(file_list), "files to process."))

# Define Area of Interest (W, E, S, N) again, as the custom function expects (Min Lon, Max Lon, Min Lat, Max Lat)
wio_aoi <- c(37, 47, -15, 1) 

# Create an empty list to store the daily Chlor_a arrays
daily_chl_list <- list()

#' Read a MODIS level-3 binned (L3b) file
#'
#' This function will read a MODIS L3b file and return a `CFData` object with
#' the data in physical units and registered to a latitude-longitude coordinate
#' system on an an authalic sphere of radius 6,378,145 m.
#'
#' This is very crude code to test the processing of L3b files. There are no
#' significant checks to capture and correct for edge cases or other potential
#' problems. The code is also not very performant, emphasis is on proof-of-concept.
#' Report any problems in the GitHub issue dedicated to building this
#' functionality into ncdfCF: https://github.com/R-CF/ncdfCF/issues/4.
#'
#' Once tested, this code is likely to be merged into ncdfCF.
#'
#' This code requires version 0.2.1 or higher of ncdfCF and library abind.
#'
#' @param fn Fully qualified name of the MODIS L3b netCDF file.
#' @param aoi Subset of data to extract, vector of four values minimum and
#' maximum longitude, minimum and maximum latitude.
#'
#' @return A `CFData` object.
#' @export
#'
#' @noRd
read_l5b <- function(fn, aoi = c(-180, 180, -90, 90)) {
  require(abind)
  
  ds <- open_ncdf(fn)
  if (length(names(ds)))
    stop("This file has recognized data variables. It is probably not a L3b formatted file. Quitting.")
  
  units <- ds$root$attribute("units")
  if (nzchar(units)) {
    units <- strsplit(units, ":")[[1L]]
    variable <- units[1L]
    units <- units[2L]
  } else stop("Required attribute 'units' not found. Quitting.")
  
  grp <- ds$root$subgroups[["level-3_binned_data"]]
  if (is.null(grp) || length(names(grp$NCvars)) != 3L)
    stop("This file does not have the required netCDF variables for L3b data. Quitting.")
  
  # Read the data
  h <- grp$handle
  binList <- RNetCDF::var.get.nc(h, "BinList")
  binIndex <- RNetCDF::var.get.nc(h, "BinIndex")
  binData <- RNetCDF::var.get.nc(h, variable)
  
  # Calculate global center lat-long and bounds
  numRows <- length(binIndex$start_num)
  lat <- (1:numRows - 0.5) * 180 / numRows - 90
  lat_bnds <- 0:numRows * 180 / numRows - 90
  
  lon_bins <- binIndex$max[numRows * 0.5]
  lon <- (1:lon_bins - 0.5) * 360 / lon_bins - 180
  lon_bnds <- 0:lon_bins * 360 / lon_bins - 180
  
  # Get the binned data in rows
  data_row <- findInterval(binList$bin_num, binIndex$start_num)
  data_row_bin <- as.integer(binList$bin_num - binIndex$start_num[data_row]) + 1L
  
  # Subset the latitude
  startRow <- floor((aoi[3L] + 90) * numRows / 180) + 1L
  endRow <- floor((aoi[4L] + 90) * numRows / 180)
  aoi_rows <- endRow - startRow + 1L
  lat <- lat[startRow:endRow]
  lat_bnds <- lat_bnds[startRow:(endRow+1L)]
  
  # Subset the longitude
  startBin <- floor((aoi[1L] + 180) * lon_bins / 360) + 1L
  endBin <- floor((aoi[2L] + 180) * lon_bins / 360)
  aoi_bins <- endBin - startBin + 1L
  lon <- lon[startBin:endBin]
  lon_bnds <- lon_bnds[startBin:(endBin+1L)]
  
  l <- lapply(startRow:endRow, function(r) {
    out <- rep(NA_real_, lon_bins)
    idx <- which(data_row == r)
    if (!length(idx)) return(out[startBin:endBin])
    
    d <- binData$sum[idx] / binList$weights[idx]   # The data of the physical property
    b <- data_row_bin[idx]                         # The bins to put the data in
    
    # Expand to lon_bins grid cells
    allocate <- ceiling(1:lon_bins * binIndex$max[r] / lon_bins)
    
    # Allocate the data to the expanded grid cells
    for (bin in seq_along(b))
      out[which(allocate == b[bin])] <- d[bin]
    
    out[startBin:endBin]
  })
  data <- abind::abind(l, along = 2)
  
  out_group <- makeMemoryGroup(-1L, "Memory_group", "/Memory_group", paste("Regridding L5b file", ds$name),
                               paste0(format(Sys.time(), "%FT%T%z"), " R package ncdfCF(", packageVersion("ncdfCF"), ")::read_l5b()"))
  
  axis_lon <- makeLongitudeAxis(-1L, "longitude", out_group, aoi_bins, lon, rbind(lon_bnds[-aoi_bins], lon_bnds[-1L]), "degrees_east")
  axis_lat <- makeLatitudeAxis(-1L, "latitude", out_group, aoi_rows, lat, rbind(lat_bnds[-aoi_rows], lat_bnds[-1L]), "degrees_north")
  
  atts <- ds$attributes()
  atts <- atts[!(atts$name %in% c("data_bins", "percent_data_bins", "binning_scheme")), ]
  atts[atts$name == "title", "value"] <- paste0(atts[atts$name == "title", ]$value[[1L]], " - regridded")
  atts[atts$name == "units", "value"] <- units
  atts[atts$name == "processing_level", "value"] <- paste0(atts[atts$name == "processing_level", ]$value[[1L]], "; regridded to lat-long using ncdfCF::read_l5b()")
  
  crs <- 'GEOGCRS["Unknown",DATUM["unknown",ELLIPSOID["unknown",6378145,0,LENGTHUNIT["metre",1]]],PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433]],CS[ellipsoidal,2],AXIS["Geodetic latitude (Lat)",north,ORDER[1]],AXIS["Geodetic longitude (Lon)",east,ORDER[2]],ANGLEUNIT["degree",0.0174532925199433]]'
  CFData$new(variable, out_group, data, list(longitude = axis_lon, latitude = axis_lat), crs, atts)
}

# --- Your original code (re-run this after defining the function) ---

# ... (Libraries and file paths should already be run) ...

# Define Area of Interest (W, E, S, N) again
wio_aoi <- c(37, 47, -15, 1) 
daily_chl_list <- list()

for (i in seq_along(file_list)) {
  fn <- file_list[i]
  
  # 1. Read and regrid the L3b file (This will now work!)
  chl_cfdata <- read_l5b(fn, aoi = wio_aoi)
  
  # 2. Extract the data as a 2D array
  chl_array <- chl_cfdata$array()
  
  # 3. Store the 2D array in the list
  daily_chl_list[[i]] <- chl_array
  
  # Optional: Print progress
  if (i %% 100 == 0) {
    cat(paste("Processed:", i, "of", length(file_list), "files...\n"))
  }
}

# 4. Join all the 2D daily arrays into a single 3D array (X, Y, Time)
chl_3d_array <- abind::abind(daily_chl_list, along = 3)

# Calculate the mean across the time dimension (3rd dimension)
# na.rm = TRUE ensures that the mean is calculated even if some days have no data
mean_chl_array <- apply(chl_3d_array, 1:2, mean, na.rm = TRUE)

# The result is a 2D array representing the 6-year mean Chlorophyll concentration
# for your WIO study area.

# assume master_location is your station table and all_data exists
master_location$Date <- as.Date(master_location$Date)  # or parse with lubridate if needed
summary(master_location$Date)  # check result

# quick bounding-box check:
range(all_data$lat); range(all_data$lon)
range(master_location$Latitude); range(master_location$Longitude)

# count how many stations lie inside satellite bounding box:
inside_idx <- with(master_location,
                   Latitude >= min(all_data$lat) & Latitude <= max(all_data$lat) &
                     Longitude >= min(all_data$lon) & Longitude <= max(all_data$lon))
table(inside = inside_idx)
