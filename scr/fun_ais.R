#-----------------------------------------------------------------------------------
# fun_ais      Suite of functions to process AIS data
#-----------------------------------------------------------------------------------
# countToDensity
# filterSAIS        Filter S-AIS data
# plotDensMol       plot density map in mollweide
# point_on_land     Check if location overlap with landmask
# trans2            Transform and normalize two rasters
#-----------------------------------------------------------------------------------






#-----------------------------------------------------------------------------------
# countToDensity      
#-----------------------------------------------------------------------------------
# pol:    polygon grid layer with number of vessels
# r_area: raster with area values per cell
# var:    variable from polygon layer to rasterize
countToDensity <- function(pol, r_area, var){
  
  # rasterterize
  r <- fasterize(pol, r_area, field = var)
  
  # add zero values
  #r[is.na(r)] <- 0
  
  # set 0 values as NA
  r[r==0] <- NA
  
  # mask cells in ocean mask
  #r <- mask(r, r_area)
  
  # calculate density
  r_dens <- r / r_area
  
  return(r_dens)
}
#-----------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------
# filterSAIS          Filter S-AIS data
#-----------------------------------------------------------------------------------
filterSAIS <- function(pol = ais_sf, r_area = r_area){
  
  # filter out cells with outlier speeds
  # exploration of data show that cells are found in isolated cells, borders, or anomalies
  rsog_mean <- fasterize(pol, r_area, field = "SOG_mean")
  rsog <- rsog_mean
  q99 <- quantile(rsog, 0.99)
  rsog[rsog >= q99] <- NA
  rsog[rsog < q99] <- 1
  
  # combine with ocean mask
  ocean_mask <- r_area/r_area
  mask1 <- rsog * ocean_mask
  
  # detect clumps of connected cells
  # after visual inspection, we select the two largers clumps
  # the rest are small clump on land, isolated cells at sea
  rclump <- clump(mask1, directions=8) 
  f <- data.frame(freq(rclump))
  f <- f[-nrow(f),]  # remove last row with NA values
  max_clump <- f$value[f$count>150]
  rclump[!rclump %in% max_clump] <- NA  # filter out the rest of the clumps
  
  # detect large clumps without variation in speed
  # we have find one in 202004 in the Artic. Likely to be sporeous data.
  novar <- zonal(rsog_mean, rclump, fun="sd") %>%
    data.frame() %>%
    dplyr::filter(sd == 0) %>%
    dplyr::select(zone) 
  if (length(novar)>0) rclump[rclump %in% novar] <- NA  # filter out the rest of the clumps
 
  # same but considering the variation in the number of total vessels
  rcount <- fasterize(pol, r_area, field = "COUNT")
  novar <- zonal(rcount, rclump, fun="sd") %>%
    data.frame() %>%
    dplyr::filter(sd == 0) %>%
    dplyr::select(zone) 
  if (length(novar)>0) rclump[rclump %in% novar] <- NA  # filter out the rest of the clumps
  
  # set all clump ids into 1 to create the mask
  rclump[rclump %in% max_clump] <- 1  # set remaining clumps to 1
  return(rclump)
}
#-----------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------
# plotDensMol       plot density map in mollweide
#-----------------------------------------------------------------------------------
plotDensMol <- function(r, zlim, logT, mollT, col, main, axis_at, axis_labels, legend_horizontal = TRUE, breaks=NULL){
  
  # import landmask
  data(countriesHigh, package = "rworldxtra", envir = environment())
  countriesHigh <- spTransform(countriesHigh, CRS("+proj=moll +ellps=WGS84"))
  box <- bb(xmin = -180, xmax = 180, ymin = -90, ymax = 90, crs="+proj=moll +ellps=WGS84")
  
  
  # labels for legend
  zlim_lab <- zlim  
  
  # reproject to mollweide
  if(mollT == TRUE){
    r <- projectRaster(r, crs = "+proj=moll +ellps=WGS84", method="bilinear")
    countriesHigh <- spTransform(countriesHigh, CRS("+proj=moll +ellps=WGS84"))
    box <- bb(xmin = -180, xmax = 180, ymin = -90, ymax = 90, crs="+proj=moll +ellps=WGS84")
  }
  
  # log-transform if required
  if(logT == TRUE){
    r = log(r)
    zlim = log(zlim)
    axis_at = log(axis_at)
  }
  
  # plot
  plot(r, maxpixels=1036800, col=col, zlim = zlim, main=main, legend=F, axes=FALSE, box=FALSE, breaks=breaks)  # raster plot
  plot(countriesHigh, col="grey80", border="grey80", add=TRUE)  # land mask
  if(mollT == TRUE) plot(box, border="grey60", add=TRUE)  # box
  plot(box, border="grey60", add=TRUE)  # box
  plot(r, zlim =  zlim, legend.only=TRUE, horizontal = legend_horizontal, col=col, legend.width=1, legend.shrink=0.4,
       axis.args=list(at=axis_at, labels=axis_labels, cex.axis=1.2))
  # legend.args=list(text=expression(Traffic~density~(vessels~km^-2)),
  #                  side=1, font=2, line=2, cex=1.2))
}
#-----------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------
# plotDens       plot density map in mollweide
#-----------------------------------------------------------------------------------
plotDens <- function(r, zlim, logT, mollT, col, main, axis_at, axis_labels, legend_horizontal = TRUE, breaks=NULL){
  
  # import landmask
  data(countriesHigh, package = "rworldxtra", envir = environment())

  # labels for legend
  zlim_lab <- zlim  
  
  # reproject to mollweide
  if(mollT == TRUE){
    r <- projectRaster(r, crs = "+proj=moll +ellps=WGS84", method="bilinear")
    countriesHigh <- spTransform(countriesHigh, CRS("+proj=moll +ellps=WGS84"))
    box <- bb(xmin = -180, xmax = 180, ymin = -90, ymax = 90, crs="+proj=moll +ellps=WGS84")
  }
  
  # log-transform if required
  if(logT == TRUE){
    r = log(r)
    zlim = log(zlim)
    axis_at = log(axis_at)
  }
  
  # plot
  plot(r, maxpixels=1036800, col=col, zlim = zlim, main=main, legend=F, axes=FALSE, box=FALSE, breaks=breaks)  # raster plot
  plot(countriesHigh, col="grey80", border="grey80", add=TRUE)  # land mask
  if(mollT == TRUE) plot(box, border="grey60", add=TRUE)  # box
  plot(r, zlim =  zlim, legend.only=TRUE, horizontal = legend_horizontal, col=col, legend.width=1, legend.shrink=0.4,
       axis.args=list(at=axis_at, labels=axis_labels, cex.axis=1.2))
  # legend.args=list(text=expression(Traffic~density~(vessels~km^-2)),
  #                  side=1, font=2, line=2, cex=1.2))
}
#-----------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------------
# point_on_land    Check if location overlap with landmask
#----------------------------------------------------------------------------------------------
point_on_land <- function(lon, lat, land = NULL){
  
  require(sp)
  require(maptools)
  require(maps)
  
  # If no landmask provided, get continent map
  if(is.null(land)){
    land <- map("world", fill = TRUE, plot = FALSE)
    land <- map2SpatialPolygons(land, IDs=land$names,proj4string=CRS("+proj=longlat +ellps=WGS84"))
  }
  
  # Convert to spatial point
  xy <- cbind(lon,lat)
  pts <- SpatialPoints(xy, proj4string=land@proj4string)
  
  # Overlay points with continents
  ov <- over(pts, land)  # returns NA when no overlap, and poly ID when there is an overlap
  
  # Return output
  return(ov)
}
#----------------------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------
# trans          Transform and normalize two rasters
#-----------------------------------------------------------------------------------
# This function log-transform and normalizes two rasters
# It uses common maximum and minum values to rescale
trans <- function(x, y){
  

  s <- stack(x, y)   # create stack
  slog <- log1p(s)  # log-transform
  maxval <- mean(quantile(slog, probs=c(0.9999)))
  minval <- min(minValue(slog))
  s_norm <- (slog - minval)/(maxval-minval)  # normalize
  s_norm[s_norm > 1] <- 1  # values are capped to a maximum value of 1
  names(s_norm) <- c("x", "y")
  return(s_norm)
  
  # # log transform
  # xlog <- log(x)
  # ylog <- log(y)
  # 
  # # normalize with the 99.99th quantile (Halpern et al. 2019)
  # vals <- c(values(xlog), values(ylog))
  # minval <- quantile(vals, probs=c(0.0001),na.rm=TRUE)
  # maxval <- quantile(vals, probs=c(0.9999),na.rm=TRUE)
  # x_norm <- (xlog - minval)/(maxval-minval)
  # y_norm <- (ylog - minval)/(maxval-minval)
  # 
  # 
  # # prepare output
  # s <- stack(x_norm, y_norm)
  # names(s) <- c("x", "y")
  # return(s)
}
#========================================================#


