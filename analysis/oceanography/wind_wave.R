

# derive frequency of high wind and waves on a monthly basis
# resolution for wave product is 0.5 x 0.5 degrees

library(raster)

nc_file <- "data/input/CDS/cds_hourly_wave.nc"
r <- brick(nc_file)

# extract dates
# get month-year
# reclass all values into 0-1 according to thresholds defined
# sum values per month-year and divide per total observations per month
# calculate difference between equivalent months