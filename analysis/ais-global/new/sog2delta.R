#---------------------------------------------------------------------------
# 03_dens2delta           Calculate absolote difference maps
#---------------------------------------------------------------------------
# This script proceses density maps


library(raster)
library(fasterize)
library(sf)
library(pals)
source("scr/fun_common.R")  # bb()
source("scr/fun_ais.R")


# set input directory
input_dir <- "data/out/ais-global/sog/"
input_dir_dens <- "data/out/ais-global/density/"

# create output directory
out_dir <- "data/out/ais-global/sog/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# select variables to compare
vars <- c("SOG")


# select months to process

### Option 1: select months to process (2020 vs 2019)
dates_post <- c(
  seq.Date(as.Date("2020-01-01"), as.Date("2020-06-01"), by = "month")
) %>% format("%Y%m%d")

dates_pre <- c(
  seq.Date(as.Date("2019-01-01"), as.Date("2019-06-01"), by = "month")
) %>% format("%Y%m%d")




#-------------------------------
# 1. Calculate changes (Mollweide)
#-------------------------------

for (i in 1:length(dates_post)){
  
  print(paste("Processing month", i, "from", length(dates_post)))
  
  # get month
  idate_pre <- dates_pre[i]
  idate_post <- dates_post[i]
  
  # create output folder for selected month
  out_dir_month <- paste0(out_dir, idate_post, "_", idate_pre, "/")
  if (!dir.exists(out_dir_month)) dir.create(out_dir_month, recursive = TRUE)

  for (j in 1:length(vars)){
  
    # select files
    # use projection in original lon/lat to make the subtraction
    jvar <- vars[j]
    tif1 <- list.files(input_dir, full.names = TRUE, recursive = TRUE, pattern = sprintf("SOG_%s_avg.tif", idate_pre))
    tif2 <- list.files(input_dir, full.names = TRUE, recursive = TRUE, pattern = sprintf("SOG_%s_avg.tif", idate_post))
    
    # import raster
    r1 <- raster(tif1)
    r2 <- raster(tif2)
    
    # calculate delta
    delta <- r2 - r1
    
    # percentage change
    per <- (delta / r1) * 100
  
    # export raster
    #writeRaster(cv, paste0(out_dir_month, sprintf("%s_%s_%s_cv_mol.tif", jvar, idate_post, idate_pre)), overwrite=TRUE)
    #writeRaster(delta, paste0(out_dir_month, sprintf("%s_%s_%s_delta.tif", jvar, idate_post, idate_pre)), overwrite=TRUE)
    #writeRaster(per, paste0(out_dir_month, sprintf("%s_%s_%s_per.tif", jvar, idate_post, idate_pre)), overwrite=TRUE)  
    

    # create mask of high traffic
    # import density maps
    tif1 <- list.files(input_dir_dens, full.names = TRUE, recursive = TRUE, pattern = sprintf("%s_COUNT_dens.tif$", idate_pre))
    tif2 <- list.files(input_dir_dens, full.names = TRUE, recursive = TRUE, pattern = sprintf("%s_COUNT_dens.tif$", idate_post))
    # import raster
    r1 <- raster(tif1)
    r2 <- raster(tif2)
    # quantiles
    q <- quantile(r1, prob = 0.8)
    high19 <- reclassify(r1, c(-Inf, q, NA, q, Inf, 1))  # create a ocean mask
    q <- quantile(r2, prob = 0.8)
    high20 <- reclassify(r2, c(-Inf, q, NA, q, Inf, 1))  # create a ocean mask
    high <- high19 * high20

    # mask delta
    delta <- mask(delta, high)

    # transform into Mollweide
    delta_mol <- projectRaster(delta, res = 27750, crs = "+proj=moll", method="bilinear")
    writeRaster(delta_mol, paste0(out_dir_month, sprintf("%s_%s_%s_delta_mol.tif", jvar, idate_post, idate_pre)), overwrite=TRUE)
    
    #### Plot delta
    pngfile <- paste0(out_dir_month, sprintf("%s_%s_%s_delta_mol.png", jvar, idate_post, idate_pre))
    png(pngfile, width=3000, height=1750, res=300)
    plotDelta(delta_mol, main = sprintf("Delta SOG (%s - %s)", idate_post, idate_pre))
    dev.off()
  }
}


