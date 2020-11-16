#-------------------------------------------------------------------------------------------------
# temporal_delta.R        Assess temporal variability of deltas across multiple months
#-------------------------------------------------------------------------------------------------
# Calculate sum across all period
# Calculate average and SD
# Find month with higher decrease

source("scr/fun_common.R")  # bb()
source("scr/fun_ais.R")


# set input directory
input_dir <- "data/out/ais-global/delta/"

# create output directory
out_dir <- "data/out/ais-global/delta_summary/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


# import ocean mask
r_area <- raster("data/out/ais-global/ocean_area.tif")
r_area <- projectRaster(r_area, res = 27750, crs = "+proj=moll +ellps=WGS84", method="ngb")

# select variables to compare
vars <- c("COUNT", "FISHING", "PASSENGER", "CARGO", "TANKER", "OTHER")

# set dates
dates_post <- c(
  seq.Date(as.Date("2020-01-01"), as.Date("2020-06-01"), by = "month")
) %>% format("%Y%m%d")

dates_pre <- c(
  seq.Date(as.Date("2019-01-01"), as.Date("2019-06-01"), by = "month")
) %>% format("%Y%m%d")


# import data



for (j in 1:length(vars)){
  
  # select variable
  jvar <- vars[j]
  s <- stack()
  
  # import data
  for(i in 1:length(dates_post)){
    
    idate_post <- dates_post[i]
    idate_pre <- dates_pre[i]
    
    tif_file <- list.files(input_dir, full.names = TRUE, recursive=TRUE, pattern=sprintf("%s_%s_%s_delta_mol.tif", jvar, idate_post, idate_pre))
    r <- raster(tif_file)
    s <- stack(s, r)
  }
  
  # Accumulate change through the January-June 2020 period
  # Calculate sum
  delta_sum <- sum(s, na.rm=TRUE) %>% mask(r_area)
  #delta_sum1 <- calc(delta_sum, fun=function(x) sign(x)*log1p(abs(x)))
  #delta_sum2 <- calc(delta_sum, fun=function(x) asinh(exp(5)*x)) #exp(2)*
  
  # Calculate average and SD
  delta_u <- mean(s, na.rm=TRUE)
  delta_sd <- calc(s, fun=sd, na.rm=TRUE)
  
  # Export data
  writeRaster(delta_sum, paste0(out_dir, sprintf("%s_delta_sum.tif", jvar)), overwrite=TRUE)
  writeRaster(delta_u, paste0(out_dir, sprintf("%s_delta_u.tif", jvar)), overwrite=TRUE)
  writeRaster(delta_sd, paste0(out_dir, sprintf("%s_delta_sd.tif", jvar)), overwrite=TRUE)
  
  
  #### Plot accumulated
  
  delta <- delta_sum
  
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
  pngfile <- paste0(out_dir, sprintf("%s_delta_accumulated.png", jvar))
  png(pngfile, width=3000, height=1750, res=300)
  plotDensMol(r = delta, zlim = c(qmin, qmax), mollT = FALSE, logT = FALSE, breaks=breaks,
              col = cols, main = sprintf("Accumulated Delta Jan-Jun 2020 (%s)", jvar),
              axis_at = c(qmin, 0,qmax),
              axis_labels = c(paste("<",round(qmin,2)), 0, paste(">",round(qmax,2))),
              legend_horizontal=TRUE)
  dev.off()
  
}




# This transformation is useful to transform skewed data that contain negative values or zeros



delta_trans2 <- calc(delta, fun=function(x) IHS(x,theta=50))
plot(delta_trans2)
hist(delta_trans2)
