#-------------------------------------------------------------------------
# 04_sector_contribution
#-------------------------------------------------------------------------


## load libraries
library(raster)
library(dplyr)
library(lubridate)
source("https://raw.githubusercontent.com/dmarch/abigoos/master/R/utils.R")  # bb()
source("scr/fun_ais.R")

# set input dir
input_dir <- "data/out/ais-global/delta/"

# create output directory
out_dir <- "data/out/ais-global/contrib/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# select variables to summarize
vars <- c("FISHING", "PASSENGER", "CARGO", "TANKER", "OTHER")

# set dates
dates_post <- c(
  seq.Date(as.Date("2020-01-01"), as.Date("2020-06-01"), by = "month")
) %>% format("%Y%m%d")

dates_pre <- c(
  seq.Date(as.Date("2019-01-01"), as.Date("2019-06-01"), by = "month")
) %>% format("%Y%m%d")



## Calculate percentage on a monhtly basis
for (i in 1:length(dates_post)){

  
  # select dates
  idate_post <- dates_post[i]
  idate_pre <- dates_pre[i]
  
  # monthly stack
  ms <- stack()  
  
  # import data
  for(j in 1:length(vars)){
    
    # select variable
    jvar <- vars[j]

    # import raster
    tif_file <- list.files(input_dir, full.names = TRUE, recursive=TRUE, pattern=sprintf("%s_%s_%s_delta_mol.tif", jvar, idate_post, idate_pre))
    r <- raster(tif_file)
    
    # create stack
    ms <- stack(ms, r)
  }

  # calculate % contribution of each sector (delta values)
  sabs <- abs(ms)  # calculate absolute values
  ssum <- sum(sabs, na.rm=TRUE)  # sum values
  sprop <- 100 * (sabs/ssum)
  names(sprop) <- vars
  sprop <- setZ(sprop, rep(idate_post, nlayers(sprop)))
  
  # for each layer, save raster and plot
  for(k in 1:nlayers(sprop)){
    
    # extract layer
    r <- subset(sprop, k)
    ivar <- names(r)
    
    # save raster
    writeRaster(r, paste0(out_dir, sprintf("contrib_%s_%s.tif", idate_post, ivar)), overwrite=TRUE)
  }
}



## Calculate percentage average during all period

for(j in 1:length(vars)){

  # select variable
  jvar <- vars[j]
  
  # import raster
  tif_file <- list.files(out_dir, full.names = TRUE, recursive=FALSE, pattern=paste0(jvar, ".tif$"))
  r <- stack(tif_file)
  
  # calculate average
  delta_u <- mean(r, na.rm=TRUE)
  delta_sd <- calc(r, fun=sd, na.rm=TRUE)
  delta_cv <- delta_sd/delta_u
  
  # plot mean
  minval <- 0
  maxval <- 100
  pngfile <- paste0(out_dir, sprintf("%s_contrib.png", jvar))
  png(pngfile, width=3000, height=1750, res=300)
  plotDensMol(r = delta_u, zlim = c(minval, maxval), mollT = FALSE, logT = FALSE,
              col = rev(brewer.spectral(101)), main = sprintf("%s Contribution (Percentage)", jvar),
              axis_at = c(minval, maxval), axis_labels = c(minval, maxval))
  dev.off()
  
  # plot sd
  minval <- 0
  maxval <- maxValue(delta_sd)
  pngfile <- paste0(out_dir, sprintf("%s_contrib_sd.png", jvar))
  png(pngfile, width=3000, height=1750, res=300)
  plotDensMol(r = delta_sd, zlim = c(minval, maxval), mollT = FALSE, logT = FALSE,
              col = rev(brewer.spectral(101)), main = sprintf("%s Contribution (Percentage - SD)", jvar),
              axis_at = c(minval, maxval), axis_labels = c(minval, maxval))
  dev.off()
  
  # plot cv
  minval <- 0
  maxval <- maxValue(delta_cv)
  pngfile <- paste0(out_dir, sprintf("%s_contrib_cv.png", jvar))
  png(pngfile, width=3000, height=1750, res=300)
  plotDensMol(r = delta_cv, zlim = c(minval, maxval), mollT = FALSE, logT = FALSE,
              col = rev(brewer.spectral(101)), main = sprintf("%s Contribution (Percentage - cv)", jvar),
              axis_at = c(minval, maxval), axis_labels = c(minval, round(maxval,1)))
  dev.off()
}








