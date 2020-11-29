#---------------------------------------------------------------------------
# 02_ais2density           Transform AIS data into densities
#---------------------------------------------------------------------------
# This script proceses S-AIS data from exact Earth
# Input data: exact Earth shapefiles (stored in different folder per month)
# Output data: marine traffic density map per ship type. Includes TIF files
# (Mollweide and lonlat) and one PNG plot.

library(raster)
library(fasterize)
library(sf)
library(pals)
source("scr/fun_common.R")  # bb()
source("scr/fun_ais.R")


# define input directory with sAIS data
input_dir <- "data/input/exacEarth/"

# create output directory
out_dir <- "data/out/ais-global/count/"
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
  
#r_area_cell <- raster("data/out/ais-global/areacell.nc")

# select variables to compare
vars <- c("COUNT", "FISHING", "PASSENGER", "CARGO", "TANKER", "OTHER")

# select months to process
dates <- c(
  seq.Date(as.Date("2019-01-01"), as.Date("2019-06-01"), by = "month"),
  seq.Date(as.Date("2020-01-01"), as.Date("2020-06-01"), by = "month")
) %>% format("%Y%m%d")

#-------------------------------
# 2. Process traffic density
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
  #m_moll <- projectRaster(m, res = 27750, crs = "+proj=moll", method="ngb")
  #writeRaster(m, paste0(out_dir, sprintf("%s_mask.tif", idate)), overwrite = TRUE)
  
  # create output folder for selected month
  out_dir_month <- paste0(out_dir, idate, "/")
  if (!dir.exists(out_dir_month)) dir.create(out_dir_month, recursive = TRUE)
  
  # process for each variable
  for (j in 1:length(vars)){
    
    # convert counts in polygon to raster format
    jvar <- vars[j]
    count <- fasterize(ais_sf, r_area, field = jvar)

    # mask with filter layer
    r <- mask(count, m)
    
    # fill gaps with 0
    r[is.na(r)] <- 0
    
    # mask with ocean at 0.25 x 0.25
    # we also store the unmasked version for the pressure analysis
    r <- mask(r, r_area)
    
    # export raster
    # writeRaster(dens, paste0(out_dir, sprintf("%s_%s_dens_unmask.tif", idate, jvar)), overwrite=TRUE)
    writeRaster(r, paste0(out_dir_month, sprintf("%s_%s_cnt.tif", idate, jvar)), overwrite=TRUE)

    # plot density
    minval <- minValue(r)
    maxval <- maxValue(r)#quantile(r, 0.99)
    pngfile <- paste0(out_dir_month, sprintf("%s_%s_cnt.png", idate, jvar))
    png(pngfile, width=3000, height=1750, res=300)
    plotDensMol(r = r, zlim = c(minval, maxval), mollT = TRUE, logT = log1p,
                col = rev(brewer.spectral(101)), main = sprintf("Vessel counts (%s %s)", idate, jvar),
                axis_at = c(minval, maxval), axis_labels = c(round(minval,3), round(maxval,3)))
    dev.off()
  }
}

