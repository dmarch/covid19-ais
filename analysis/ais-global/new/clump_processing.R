# We keep:
# - large clumps (patches > 1000 cells)
# - small clumps wich are adjacent to the coastline (closest distance to coastline < 27,750 m)
#
# We remove:
# - clumps with no variation in SD and count


#---------------------------------------------------------------------------
# clump_processing        Process clumps
#---------------------------------------------------------------------------
# This script analyses clumps from S-AIS data from exact Earth 
# After visual inspection, we detected isolated cells or cells with anomalous
# data. To provide a QC analysis, we detected clumps and characterized them by
# different factors:
# ocean area size: very small patches are likely to be the results of spureous detection
# closest distance to shore: isolated patches near the coast can be the result of local movements
# variability in speed and number of vessels: visual inspection indicates that most the
# clumps in higher latitudes have exactly same speed and number of vessels. Such pattern is
# rare, and we can use the SD to identify those clumps.

library(raster)
library(fasterize)
library(sf)
library(pals)
source("scr/fun_common.R")  # bb()
source("scr/fun_ais.R")


# define input directory with sAIS data
input_dir <- "data/input/exacEarth/"

# create output directory
out_dir <- "data/out/ais-global/clumps/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


#-------------------------------
# 1. Import data
#-------------------------------

# import raster
dist2coast <- raster("data/out/ais-global/dist2coast.tif")


# import ocean mask
r_area <- raster("data/out/ais-global/ocean_area.tif")
#r_area_cell <- raster("data/out/ais-global/areacell.nc")

# select months to process
dates <- c(
  seq.Date(as.Date("2019-01-01"), as.Date("2019-06-01"), by = "month"),
  seq.Date(as.Date("2020-01-01"), as.Date("2020-06-01"), by = "month")
) %>% format("%Y%m%d")

#-----------------------------------------
# 2. Detect and get info for each clump
#-----------------------------------------

clump_list <- list()

for (i in 1:length(dates)){
  
  print(paste("Processing month", i, "from", length(dates)))
  
  # import shapefile
  idate <- dates[i]
  shp_file <- list.files(input_dir, full.names = TRUE, recursive = TRUE, pattern = sprintf("%s.*.shp$", idate))
  
  # import shapefile
  ais_sf <- st_read(shp_file)
  
  # filter sAIS to get mask
  iclumps <- getClumpInfo(pol = ais_sf, r_area = r_area, dist2coast = dist2coast)
  iclumps$date <- idate
  clump_list[[i]] <- iclumps
}

# combine data
clump_data <- data.table::rbindlist(clump_list)

# export file
outfile <- paste0(out_dir, "clumps.csv")
write.csv(clump_data, outfile, row.names=FALSE)


#-----------------------------------------
# 3. Exploratory analysis
#-----------------------------------------

nclumps <- nrow(clump_data)

# get number per different cases
filter(clump_data, count > 1000) %>% nrow()
filter(clump_data, count == 1) %>% nrow()/nrow(clump_data)
filter(clump_data, sd_speed == 0 | sd_count == 0) %>% nrow()/nrow(clump_data)
filter(clump_data, grid_area < 769) %>% nrow()/nrow(clump_data)




# first, we filter out clumps with only one cell and large clumps
clump_data <- filter(clump_data,
                     count > 1,
                     count < 1000,
                     #dist2coast >= 27750,
                     grid_area > 769,
                     sd_speed > 0,
                     sd_count > 0)


# histogram or boxplot per variable
boxplot(clump_data$sd_speed)
boxplot(clump_data$sd_count)
boxplot(clump_data$area)
boxplot(clump_data$dist2coast)

hist(clump_data$sd_speed)
hist(clump_data$sd_count)
hist(clump_data$area)
hist(clump_data$dist2coast)


# set threshold for size
q25 <- quantile(clump_data$area, prob=.25)  # km2
# 1398.438
