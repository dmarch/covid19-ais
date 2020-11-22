#---------------------------------------------------------------------------
# ais2sog          Transform AIS data into average Speed Over Ground
#---------------------------------------------------------------------------
# This script proceses S-AIS data from exact Earth


library(raster)
library(fasterize)
library(sf)
library(pals)
source("https://raw.githubusercontent.com/dmarch/abigoos/master/R/utils.R")  # bb()
source("scr/fun_ais.R")


# define input directory with sAIS data
input_dir <- "data/input/exacEarth/"

# create output directory
out_dir <- "data/out/ais-global/sog/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)



#-------------------------------
# 1. Import data
#-------------------------------

# import ocean mask
#r_area <- raster("data/out/ais-global/oceanmask.nc")
r_area <- raster("data/out/ais-global/ocean_area.tif")

# use ocean area as mask
r_area_grid <- area(r_area)
r_area <- mask(r_area_grid, r_area)

# select variables to compare
vars <- c("SOG_mean")

# select months to process
dates <- c(
  seq.Date(as.Date("2019-01-01"), as.Date("2019-07-01"), by = "month"),
  seq.Date(as.Date("2020-01-01"), as.Date("2020-07-01"), by = "month")
) %>% format("%Y%m%d")



#-------------------------------
# 2. Process SOG data
#-------------------------------

for (i in 1:length(dates)){
  
  print(paste("Processing month", i, "from", length(dates)))
  
  # import shapefile
  idate <- dates[i]
  shp_file <- list.files(input_dir, full.names = TRUE, recursive = TRUE, pattern = sprintf("%s.*.shp$", idate))
  
  # import shapefile
  ais_sf <- st_read(shp_file)
  
  # filter sAIS to get mask
  m <- filterSAIS2(pol = ais_sf, r_area = r_area)
  m_moll <- projectRaster(m, res = 27750, crs = "+proj=moll", method="ngb")

  # transform to raster
  r <- fasterize(ais_sf, m, field = vars)
  
  # mask cells in ocean mask
  r <- mask(r, m)
  
  # export raster
  writeRaster(r, paste0(out_dir, sprintf("SOG_%s_avg.tif", idate)), overwrite=TRUE)
  
  # plot density
  minval <- minValue(r)
  maxval <- maxValue(r)#quantile(r, 0.99)
  pngfile <- paste0(out_dir, sprintf("SOG_%s_avg.png", idate))
  png(pngfile, width=3000, height=1750, res=300)
  plotDensMol(r = r, zlim = c(minval, maxval), mollT = TRUE, logT = FALSE,
              col = rev(brewer.spectral(101)), main = sprintf("SOG (%s)", idate),
              axis_at = c(minval, maxval), axis_labels = c(round(minval,3), round(maxval,3)))
  dev.off()
}
