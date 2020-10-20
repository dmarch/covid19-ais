#---------------------------------------------------------------------------
# 01_coastal_waters                                                         
#---------------------------------------------------------------------------
# This script generates a vectorial layer that represents the coastal waters
# of the study area, considering the national jurisdictions.
# This layer is then used to extract T-AIS data
#
# Calculates the area covered in the temporal analysis
#
# Input data:
# - Bounding box of the study area
# - EEZ map
#---------------------------------------------------------------------------


## load libraries
library(sf)
library(raster)

## here we store data
output_data <- "data/out/ais-wmed"
if (!dir.exists(output_data)) dir.create(output_data, recursive = TRUE)


# bb: bounding box of the AIS data
# bb_expanded: bounding box to cover the coastline and calculate the buffer
bb <- st_as_sfc(st_bbox(c(xmin = -1.60, xmax = 9.93, ymax = 35.56, ymin = 43.45), crs = st_crs(4326)))
bb_expanded <- st_as_sfc(st_bbox(c(xmin = -3, xmax = 12, ymax = 33, ymin = 45), crs = st_crs(4326)))


# Import EEZ
# intersect with bb
eez <- st_read("data/input/marine_regions/World_EEZ_v11_20191118_gpkg/eez_v11.gpkg")
eez_wmed <- st_intersection(eez, bb)
#plot(st_geometry(eez_wmed))

# Import national boundaries GADM
# intersect with expanded bb to get the coastline
gadm <- st_read("data/input/gadm/gadm36_0.shp")
gadm_wmed <- st_intersection(gadm, bb_expanded)
#plot(st_geometry(gadm_wmed))

# Transform to projected CRS
eez_wmed <- st_transform(eez_wmed, 3035)  # transform to overlap with AIS rasters
gadm_wmed <- st_transform(gadm_wmed, 3035)  # transform to overlap with AIS rasters

# Buffer coastline 24 nautical miles
gadm_wmed_24nm <- st_buffer(gadm_wmed, 44400) %>%  st_union # 24 nautical miles, 44.4 km
plot(st_geometry(gadm_wmed_24nm))

# Intersect with study area
eez_24nm <- st_intersection(eez_wmed, gadm_wmed_24nm)
plot(eez_24nm["GEONAME"])

# Filter for EU countries
eez_24nm_eu <- eez_24nm %>% filter(SOVEREIGN1 %in% c("Italy", "France", "Spain"))

# Calculate area
area <- st_area(eez_24nm_eu)
units(area) = "km^2"
sum(area)  # 164318.2 [km^2]

# Backtransform to lonlat
eez_24nm_eu <- st_transform(eez_24nm_eu, 4326)  # transform to overlap with AIS rasters
plot(eez_24nm_eu["GEONAME"])

# save
outfile <- paste0(output_data, "/eez_24nm_eu.gpkg")
st_write(eez_24nm_eu, outfile)

