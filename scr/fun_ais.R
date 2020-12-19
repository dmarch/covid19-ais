#-----------------------------------------------------------------------------------
# fun_ais      Suite of functions to process AIS data
#-----------------------------------------------------------------------------------
# countToDensity
# filterSAIS        Filter S-AIS data
# getClumpInfo      Get information for cell patches (clumps)
# plotDelta       plot delta map in mollweide
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
# filterSAIS2          Filter S-AIS data
#-----------------------------------------------------------------------------------
filterSAIS2 <- function(pol = ais_sf, r_area = r_area){#, dist2coast = dist2coast){
  # For each density map, we detected the patches of connected cells.
  # For each cell, we consider the closest 8 cells as adjacent (Queen's case).
  # For each clump, we then calculated the ocean area, the closest distance to the coastline,
  # the SD for the average speed and SD of the total number of vessels.
  
  require(dplyr)
  
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
  
  # calculate SD for mean speed
  clump_sd_speed <- zonal(rsog_mean, rclump, fun="sd") %>%
    data.frame() %>%
    dplyr::select(zone, sd) %>%
    dplyr::rename(value = zone, sd_speed = sd)
  
  # calculate SD for number of vessels
  rcount <- fasterize(pol, r_area, field = "COUNT")
  clump_sd_count <- zonal(rcount, rclump, fun="sd") %>%
    data.frame() %>%
    dplyr::select(zone, sd) %>%
    dplyr::rename(value = zone, sd_count = sd)
  
  # calculate total number of vessels
  clump_sum_count <- zonal(rcount, rclump, fun="sum") %>%
    data.frame() %>%
    dplyr::select(zone, sum) %>%
    dplyr::rename(value = zone, sum_count = sum)
  
  # calculate ocean area
  clump_ocean_area <- zonal(r_area, rclump, fun="sum") %>%
    data.frame() %>%
    dplyr::select(zone, sum) %>%
    dplyr::rename(value = zone, ocean_area = sum)
  
  # calculate grid area
  grid_area <- area(r_area)
  clump_grid_area <- zonal(grid_area, rclump, fun="sum") %>%
    data.frame() %>%
    dplyr::select(zone, sum) %>%
    dplyr::rename(value = zone, grid_area = sum)
  
  # combine all variables
  clumps <- f %>%
    left_join(clump_sd_speed, by="value") %>%
    left_join(clump_sd_count, by="value") %>%
    left_join(clump_sum_count, by="value") %>%
    left_join(clump_ocean_area, by="value") %>%
    left_join(clump_grid_area, by="value")
  
  # filter data
  clump_data <- filter(clumps,
                       grid_area > 769,
                       sd_speed > 0,
                       sd_count > 0)
  
  # set all clump ids into 1 to create the mask
  rclump[!rclump %in% clump_data$value] <- NA  # set remaining clumps to 1
  rclump[rclump %in% clump_data$value] <- 1  # set remaining clumps to 1
  return(rclump)
}
#-----------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------
# getClumpInfo         Get information for cell patches (clumps)
#-----------------------------------------------------------------------------------
getClumpInfo <- function(pol = ais_sf, r_area = r_area, dist2coast = dist2coast){
  # For each density map, we detected the patches of connected cells.
  # For each cell, we consider the closest 8 cells as adjacent (Queen's case).
  # For each clump, we then calculated the ocean area, the closest distance to the coastline,
  # the SD for the average speed and SD of the total number of vessels.
  
  require(dplyr)
  
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
  
  # calculate SD for mean speed
  clump_sd_speed <- zonal(rsog_mean, rclump, fun="sd") %>%
    data.frame() %>%
    dplyr::select(zone, sd) %>%
    dplyr::rename(value = zone, sd_speed = sd)
  
  # calculate SD for number of vessels
  rcount <- fasterize(pol, r_area, field = "COUNT")
  clump_sd_count <- zonal(rcount, rclump, fun="sd") %>%
    data.frame() %>%
    dplyr::select(zone, sd) %>%
    dplyr::rename(value = zone, sd_count = sd)
  
  # calculate total number of vessels
  clump_sum_count <- zonal(rcount, rclump, fun="sum") %>%
    data.frame() %>%
    dplyr::select(zone, sum) %>%
    dplyr::rename(value = zone, sum_count = sum)
  
  # calculate ocean area
  clump_ocean_area <- zonal(r_area, rclump, fun="sum") %>%
    data.frame() %>%
    dplyr::select(zone, sum) %>%
    dplyr::rename(value = zone, ocean_area = sum)
  
  # calculate grid area
  grid_area <- area(r_area)
  clump_grid_area <- zonal(grid_area, rclump, fun="sum") %>%
    data.frame() %>%
    dplyr::select(zone, sum) %>%
    dplyr::rename(value = zone, grid_area = sum)
  
  # calculate closest distance to the coastline
  clump_dist2coast <- zonal(dist2coast, rclump, fun="min") %>%
    data.frame() %>%
    dplyr::select(zone, min) %>%
    dplyr::rename(value = zone, dist2coast = min)
  
  # combine all variables
  clumps <- f %>%
    left_join(clump_sd_speed, by="value") %>%
    left_join(clump_sd_count, by="value") %>%
    left_join(clump_sum_count, by="value") %>%
    left_join(clump_ocean_area, by="value") %>%
    left_join(clump_grid_area, by="value") %>%
    left_join(clump_dist2coast, by="value")
  
  return(clumps)
}
#-----------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------
# plotDelta       plot delta map in mollweide
#-----------------------------------------------------------------------------------
plotDelta <- function(delta, main = sprintf("Accumulated Delta Jan-Jun 2020 (%s)", jvar)){
  
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
  qmax <- m
  qmin <- m*(-1)
  delta[delta>qmax]<-qmax
  delta[delta<qmin]<-qmin
  
  # define regular interval
  intervals <- m/50
  
  # get breaks
  # we consider the case where delta values may result in positive values
  # (e.g. 2020-01)
  if (qmin >= 0) breaks <- seq(qmin, qmax, by = intervals)
  if (qmin < 0){
    low_breaks <- c(seq(qmin-intervals, 0-intervals, by=intervals))
    high_breaks <- c(seq(0, qmax+intervals, by=intervals))
    breaks <- c(low_breaks, high_breaks)
  }
  
  # diverging assymetric color ramp
  high <- brewer.blues(4)
  low <-  rev(brewer.reds(4))
  if (qmin >= 0) cols <- high
  if (qmin < 0){
    low_cols <- colorRampPalette(low)(length(low_breaks)-1)
    high_cols <- colorRampPalette(high)(length(high_breaks)-1)
    cols <- c(low_cols,"#f7f7f7","#f7f7f7", high_cols)
  }
  
  # figure plot
  plotDensMol(r = delta, zlim = c(qmin, qmax), mollT = FALSE, logT = FALSE, breaks=breaks,
              col = cols, main = main,
              axis_at = c(qmin, 0,qmax),
              axis_labels = c(paste("<",round(qmin,2)), 0, paste(">",round(qmax,2))),
              legend_horizontal=TRUE)
}
#-----------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------
# plotDelta       plot delta map in mollweide
#-----------------------------------------------------------------------------------
plotDelta2 <- function(delta, percentile = 0.99, main = sprintf("Accumulated Delta Jan-Jun 2020 (%s)", jvar)){
  
  # get min and max values
  mindelta <- minValue(delta)
  maxdelta <- maxValue(delta)
  
  # get 99th percentiles
  if (is.null(percentile)){
    qmin <- mindelta
    qmax <- maxdelta
    
    # define regular interval
    m <- max(abs(qmin), qmax)
    intervals <- m/50
    
  } else {
    qmin <- quantile(delta, 1-percentile)
    qmax <- quantile(delta, percentile)
    
    # define regular interval
    m <- max(abs(qmin), qmax)
    qmax <- m
    qmin <- m*(-1)
    delta[delta>qmax]<-qmax
    delta[delta<qmin]<-qmin
    
    # define regular interval
    intervals <- m/50
    
  }
  
  # get breaks
  # we consider the case where delta values may result in positive values
  # (e.g. 2020-01)
  if (qmin >= 0) breaks <- seq(qmin, qmax, by = intervals)
  if (qmin < 0){
    low_breaks <- c(seq(qmin-intervals, 0-intervals, by=intervals))
    high_breaks <- c(seq(0, qmax+intervals, by=intervals))
    breaks <- c(low_breaks, high_breaks)
  }
  
  # diverging assymetric color ramp
  high <- brewer.blues(4)
  low <-  rev(brewer.reds(4))
  if (qmin >= 0) cols <- high
  if (qmin < 0){
    low_cols <- colorRampPalette(low)(length(low_breaks)-1)
    high_cols <- colorRampPalette(high)(length(high_breaks)-1)
    cols <- c(low_cols,"#f7f7f7","#f7f7f7", high_cols)
  }
  
  # figure plot
  plotDensMol(r = delta, zlim = c(qmin, qmax), mollT = FALSE, logT = FALSE, breaks=breaks,
              col = cols, main = main,
              axis_at = c(qmin, 0,qmax),
              axis_labels = c(paste("<",round(qmin,2)), 0, paste(">",round(qmax,2))),
              legend_horizontal=TRUE)
}
#-----------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------
# plotDensMol       plot density map in mollweide
#-----------------------------------------------------------------------------------
plotDensMol <- function(r, zlim, logT, mollT, col, main, axis_at, axis_labels, legend_horizontal = TRUE, breaks=NULL){
  
  require(pals)
  
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
  # if(logT == TRUE){
  #   r = log(r)
  #   zlim = log(zlim)
  #   axis_at = log(axis_at)
  # }
  
  # log-transform if required
  if(is.function(logT)){
    r = logT(r)
    zlim = logT(zlim)
    axis_at = logT(axis_at)
  }
  
  # plot
  plot(r, maxpixels=1036800, col=col, zlim = zlim, main=main, legend=F, axes=FALSE, box=FALSE, breaks=breaks)  # raster plot
  plot(countriesHigh, col="grey80", border="grey80", add=TRUE)  # land mask
  if(mollT == TRUE) plot(box, border="grey60", add=TRUE)  # box
  plot(box, border="grey60", add=TRUE)  # box
  plot(r, zlim =  zlim, legend.only=TRUE, horizontal = legend_horizontal, col=col, legend.width=1, legend.shrink=0.4,
       axis.args=list(at=axis_at, labels=axis_labels, cex.axis=1.8))
  # legend.args=list(text=expression(Traffic~density~(vessels~km^-2)),
  #                  side=1, font=2, line=2, cex=1.2))
}
#-----------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------
# plotDensMol2       plot density map in mollweide
#-----------------------------------------------------------------------------------
plotDensMol2 <- function(r, col, main, legend_horizontal = TRUE){
  
  # import landmask
  data(countriesHigh, package = "rworldxtra", envir = environment())
  countriesHigh <- spTransform(countriesHigh, CRS("+proj=moll +ellps=WGS84"))
  box <- bb(xmin = -180, xmax = 180, ymin = -90, ymax = 90, crs="+proj=moll +ellps=WGS84")
  
  # log-transform
  r = log10(r)
  zlim = c(minValue(r), maxValue(r))
  axis_at = seq(floor(zlim[1]), ceiling(zlim[2]))
  axis_labels = paste(10, '^', axis_at, sep='') 
  
  # plot
  plot(r, maxpixels=1036800, col=col, zlim = zlim, main=main, legend=F, axes=FALSE, box=FALSE)  # raster plot
  plot(countriesHigh, col="grey80", border="grey80", add=TRUE)  # land mask
  plot(box, border="grey60", add=TRUE)  # box
  plot(r, zlim =  zlim, legend.only=TRUE, horizontal = legend_horizontal, col=col, legend.width=1, legend.shrink=0.4,
       axis.args=list(at=axis_at, labels=parse(text=axis_labels), cex.axis=1.8))
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


