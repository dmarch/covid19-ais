#---------------------------------------------------------------------------
# 03_dens2delta           Calculate absolote difference maps
#---------------------------------------------------------------------------
# This script proceses density maps


library(raster)
library(fasterize)
library(sf)
library(pals)
source("https://raw.githubusercontent.com/dmarch/abigoos/master/R/utils.R")  # bb()
source("scr/fun_ais.R")


# set input directory
input_dir <- "data/out/ais-global/density/"

# create output directory
out_dir <- "data/out/ais-global/delta/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# select variables to compare
vars <- c("COUNT", "FISHING", "PASSENGER", "CARGO", "TANKER", "OTHER")
dates <- c("20190401", "20200401")



#-------------------------------
# 1. Calculate changes (Mollweide)
#-------------------------------

for (j in 1:length(vars)){
  
  # select files
  jvar <- vars[j]
  tif1 <- list.files(input_dir, full.names = TRUE, pattern = sprintf("%s_%s_dens_moll.tif$", dates[1], jvar))
  tif2 <- list.files(input_dir, full.names = TRUE, pattern = sprintf("%s_%s_dens_moll.tif$", dates[2], jvar))
  
  # import raster
  r1 <- raster(tif1)
  r2 <- raster(tif2)
  
  # create common mask between both rasters
  # identifies cells with presence of route density in any of the two periods
  m <- sum(r1, r2, na.rm=TRUE)
  m[m==0] <- NA
  m[m>0] <- 0
  
  # add zero to density maps
  r1z <- sum(r1, m, na.rm=TRUE)
  r2z <- sum(r2, m, na.rm=TRUE)
  
  ## calculate coefficient of variation (CV)
  s <- stack(r1z, r2z)
  u <- mean(s, na.rm=TRUE)
  s <- calc(s, fun=sd, na.rm=TRUE)
  cv <- s/u
  
  # calculate delta
  delta <- sum(r2z, r1z*(-1), na.rm=TRUE)
  delta <- mask(delta, m)
  
  # percentage change
  per <- (delta / r1z) * 100

  # export raster
  writeRaster(cv, paste0(out_dir, sprintf("%s_cv_mol.tif", jvar)), overwrite=TRUE)
  writeRaster(delta, paste0(out_dir, sprintf("%s_delta_mol.tif", jvar)), overwrite=TRUE)
  writeRaster(per, paste0(out_dir, sprintf("%s_per_mol.tif", jvar)), overwrite=TRUE)  
}


#-------------------------------
# 2. Calculate changes (WGS84)
#-------------------------------

for (j in 1:length(vars)){
  
  # select files
  jvar <- vars[j]
  tif1 <- list.files(input_dir, full.names = TRUE, pattern = sprintf("%s_%s_dens.tif$", dates[1], jvar))
  tif2 <- list.files(input_dir, full.names = TRUE, pattern = sprintf("%s_%s_dens.tif$", dates[2], jvar))
  
  # import raster
  r1 <- raster(tif1)
  r2 <- raster(tif2)

  # create common mask between both rasters
  # identifies cells with presence of route density in any of the two periods
  m <- sum(r1, r2, na.rm=TRUE)
  m[m==0] <- NA
  m[m>0] <- 0

  # add zero to density maps
  r1z <- sum(r1, m, na.rm=TRUE)
  r2z <- sum(r2, m, na.rm=TRUE)
  
  ## calculate coefficient of variation (CV)
  s <- stack(r1z, r2z)
  u <- mean(s, na.rm=TRUE)
  s <- calc(s, fun=sd, na.rm=TRUE)
  cv <- s/u
  
  # calculate delta
  delta <- sum(r2z, r1z*(-1), na.rm=TRUE)
  delta <- mask(delta, m)
  
  # percentage change
  per <- (delta / r1z) * 100

  # export raster
  writeRaster(cv, paste0(out_dir, sprintf("%s_cv.tif", jvar)), overwrite=TRUE)
  writeRaster(delta, paste0(out_dir, sprintf("%s_delta.tif", jvar)), overwrite=TRUE)
  writeRaster(per, paste0(out_dir, sprintf("%s_per.tif", jvar)), overwrite=TRUE)  
}


#-------------------------------
# 3. Plot delta
#-------------------------------


for(i in 1:length(vars)){
  
  ivar <- vars[i]
  
  # import delta
  delta <- raster(sprintf("data/out/ais-global/delta/%s_delta_mol.tif", ivar))
  
  # get min and max values
  mindelta <- minValue(delta)
  maxdelta <- maxValue(delta)
  
  # get 99th percentiles
  qmin <- quantile(delta, c(0.01))
  qmax <- quantile(delta, c(0.99))
  #qmin <- -0.05
  #qmax <- 0.05
  
  # define regular interval
  m <- max(abs(qmin), qmax)
  delta[delta>qmax]<-qmax
  delta[delta<qmin]<-qmin
  
  # define regular interval
  intervals <- m/50
  
  # get breaks
  low_breaks <- c(seq(qmin-intervals, 0-intervals, by=intervals))
  high_breaks <- c(seq(0, qmax+intervals, by=intervals))
  breaks <- c(low_breaks, high_breaks)
  
  # diverging assymetric color ramp
  high <- brewer.blues(4)
  low <-  rev(brewer.reds(4))
  low_cols <- colorRampPalette(low)(length(low_breaks)-1)
  high_cols <- colorRampPalette(high)(length(high_breaks)-1)
  cols <- c(low_cols,"#f7f7f7","#f7f7f7", high_cols)

  # figure plot
  pngfile <- sprintf("data/out/ais-global/delta/%s_delta_mol.png", ivar)
  png(pngfile, width=3000, height=1750, res=300)
  plotDensMol(r = delta, zlim = c(qmin, qmax), mollT = FALSE, logT = FALSE, breaks=breaks,
              col = cols, main = sprintf("%s (delta)", ivar),
              axis_at = c(qmin, 0,qmax),
              axis_labels = c(paste("<",round(qmin,2)), 0, paste(">",round(qmax,2))),
              legend_horizontal=TRUE)
  dev.off()
}

