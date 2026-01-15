#OCNET data 
#
# to finish...

# Source structure exploration function
source(here::here("functions", "extracting.R"))

#smaller bounding box
coords_list_smaller <- list(
  W = 35.0,   # West longitude
  E = 44.0,   # East longitude
  S = -15.0,  # South latitude
  N = 1.0     # North latitude
)

#ocnet filepath
path <- "C:/KenyaSSFProject/raw_data/ocnet_chla_data/OCNET_chla_2015.nc"


#extracting lat and lon coords from ocnet 
data_list_2018_ocnet <- extract_latlon_ocnet(path, coords_list_smaller)
