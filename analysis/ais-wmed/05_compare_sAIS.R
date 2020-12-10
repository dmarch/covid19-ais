# correlation between T-AIS and S-AIS

library(tidyverse)
library(raster)
library(rasterVis)
library(grid)

# set dates
# for each date, import T-AIS and S-AIS same date and vessel type
# for each pair, calculate pearson at t1
# for each vesel type, calculate average of pearson for all t




# set input directories
input_global <- "data/out/ais-global/density/"
input_wmed <- "data/out/ais-wmed/density/"

# create output directory
out_dir <- "data/out/ais-wmed/compare/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# select months to process
dates <- c(
  seq.Date(as.Date("2019-01-01"), as.Date("2019-07-01"), by = "month"),
  seq.Date(as.Date("2020-01-01"), as.Date("2020-07-01"), by = "month")
)

# select variables to compare
vars <- c("COUNT", "FISHING", "PASSENGER", "CARGO", "TANKER", "OTHER")

# import ocean mask
r_mask <- raster("data/out/ais-global/oceanmask.nc")
r_mask <- r_mask/r_mask

#-------------------------------
# Calculate delta and pearson
#-------------------------------

for (i in 1:length(dates)){
  
  print(paste("Processing month", i, "from", length(dates)))
  
  # select date
  idate <- format(dates[i], "%Y%m%d")
  
  # create output folder for selected month
  out_dir_month <- paste0(out_dir, idate, "/")
  if (!dir.exists(out_dir_month)) dir.create(out_dir_month, recursive = TRUE)
  
  for (j in 1:length(vars)){
    
    # select files
    jvar <- vars[j]
    tif1 <- list.files(input_global, full.names = TRUE, recursive = TRUE, pattern = sprintf("%s_%s_dens.tif$", idate, jvar))
    tif2 <- list.files(input_wmed, full.names = TRUE, recursive = TRUE, pattern = sprintf("%s_%s_dens.tif$", idate, jvar))
    
    # import raster
    r_glob <- raster(tif1)
    r_med <- raster(tif2)

    # crop global to same extent as the Mediterranean
    r_glob <- crop(r_glob,r_med)
    
    # calculate delta
    delta <- r_med - r_glob
    
    # correlation
    x <- corLocal(r_glob, r_med, ngb=5, method=c("pearson"), test=TRUE )
    x <- x * r_mask

    # export rasters
    writeRaster(delta, paste0(out_dir_month, sprintf("%s_%s_delta.tif", idate, jvar)), overwrite=TRUE)
    writeRaster(x[[1]], paste0(out_dir_month, sprintf("%s_%s_cor.tif", idate, jvar)), overwrite=TRUE)
    writeRaster(x[[2]], paste0(out_dir_month, sprintf("%s_%s_pval.tif", idate, jvar)), overwrite=TRUE)
  }
}


#--------------------------------------------------------------------------
# Plot monthly facets (density for T-AIS)
#--------------------------------------------------------------------------

# import land mask
data(countriesHigh, package = "rworldxtra", envir = environment())

# set color
col <- rasterTheme(region = rev(brewer.pal(9, 'Spectral')))

#for (j in 1:length(vars)){

  jvar <- "COUNT"

  # import rasters
  input <- input_wmed  # use "input_wmed" or "input_global"
  tif_files <- list.files(input, full.names = TRUE, recursive = TRUE, pattern = sprintf("%s_dens.tif$", jvar))
  b <- stack(tif_files)
  #b <- crop(b, r_med)
  
  # facet plot
  land <- crop(countriesHigh, b)
  p <- levelplot(b, layout=c(7, 2),
                 par.settings = col, 
                 at=seq(0, 3.5, length.out=101),
                 names.attr=format(dates, "%b %Y"),
                 colorkey=list(height=0.8, width=1)) +
    layer(sp.lines(land, col="black", lwd=0.8))
  
  
  # Save as png file
  p_png <- paste0(out_dir, sprintf("%s_density_tAIS.png", jvar))
  png(p_png, width = 24, height = 8, units = "cm", res = 300)
  print(p)
  trellis.focus("legend", side="right", clipp.off=TRUE, highlight=FALSE)
  dev.off()
  
  ## calculate coefficient of variation (CV)
  u_med <- mean(b, na.rm=TRUE)
  s_med <- calc(b, fun=sd, na.rm=TRUE)
  cv_med <- s_med/u_med
#}



#--------------------------------------------------------------------------
# Plot monthly facets (density for S-AIS)
#--------------------------------------------------------------------------

r_med <- subset(b,1)

# import rasters
input <- input_global  # use "input_wmed" or "input_global"
tif_files <- list.files(input, full.names = TRUE, recursive = TRUE, pattern = sprintf("%s_dens.tif$", jvar))
b <- stack(tif_files)
b <- crop(b, r_med)

# facet plot
p <- levelplot(b, layout=c(7, 2),
               par.settings = col, 
               at=seq(0, 3.5, length.out=101),
               names.attr=format(dates, "%b %Y"),
               colorkey=list(height=0.8, width=1)) +
  layer(sp.lines(land, col="black", lwd=0.8))


# Save as png file
p_png <- paste0(out_dir, sprintf("%s_density_sAIS.png", jvar))
png(p_png, width = 24, height = 8, units = "cm", res = 300)
print(p)
trellis.focus("legend", side="right", clipp.off=TRUE, highlight=FALSE)
dev.off()

## calculate coefficient of variation (CV)
u_glob <- mean(b, na.rm=TRUE)
s_glob <- calc(b, fun=sd, na.rm=TRUE)
cv_glob <- s_glob/u_glob



#--------------------------------------------------------------------------
# Plot mean
#--------------------------------------------------------------------------

# combine CVs
b_u <- brick(u_glob, u_med)

# plot
p <- levelplot(b_u, layout=c(2, 1),
               par.settings = col, 
               at=seq(0,  max(maxValue(b_u)), length.out=101),
               names.attr=c("S-AIS", "T-AIS"),
               colorkey=list(height=0.8, width=1)) +
  layer(sp.lines(land, col="black", lwd=0.8))


# Save as png file
p_png <- paste0(out_dir, "u.png")
png(p_png, width = 12, height = 8, units = "cm", res = 300)
print(p)
trellis.focus("legend", side="right", clipp.off=TRUE, highlight=FALSE)
dev.off()



#--------------------------------------------------------------------------
# Plot CV
#--------------------------------------------------------------------------

# combine CVs
b_cv <- brick(cv_glob, cv_med)

# plot
p <- levelplot(b_cv, layout=c(2, 1),
               par.settings = col, 
               at=seq(0,  max(maxValue(b_cv)), length.out=101),
               names.attr=c("S-AIS", "T-AIS"),
               colorkey=list(height=0.8, width=1)) +
  layer(sp.lines(land, col="black", lwd=0.8))


# Save as png file
p_png <- paste0(out_dir, "cv.png")
png(p_png, width = 12, height = 8, units = "cm", res = 300)
print(p)
trellis.focus("legend", side="right", clipp.off=TRUE, highlight=FALSE)
dev.off()


#--------------------------------------------------------------------------
# Plot monthly facets (delta)
#--------------------------------------------------------------------------

jvar = "COUNT"

# set color
col <- rasterTheme(region = brewer.pal(9, 'RdBu'))

# import rasters
tif_files <- list.files(out_dir, full.names = TRUE, recursive = TRUE, pattern = sprintf("%s_delta.tif$", jvar))
b <- stack(tif_files)
b <- crop(b, r_med)

# facet plot
p <- levelplot(b, layout=c(7, 2),
          par.settings = col,
          at=seq(-2.7, 2.7, length.out=101),
          names.attr=format(dates, "%b %Y"),
          colorkey=list(height=0.8, width=1)) +
  layer(sp.lines(land, col="black", lwd=0.8))


# Save as png file
p_png <- paste0(out_dir, "delta.png")
png(p_png, width = 24, height = 8, units = "cm", res = 300)
print(p)
trellis.focus("legend", side="right", clipp.off=TRUE, highlight=FALSE)
dev.off()

## calculate mean and SD
u <- mean(b, na.rm=TRUE)
s <- calc(b, fun=sd, na.rm=TRUE)
b_delta <- brick(u, s)

# plot mean
p <- levelplot(u,
               par.settings = col, 
               margin = FALSE,
               at=seq(-1.9, 1.9, length.out=101),
               #names.attr=c("Mean", "SD"),
               colorkey=list(height=0.8, width=1)) +
  layer(sp.lines(land, col="black", lwd=0.8))

# Save as png file
p_png <- paste0(out_dir, "mean_delta.png")
png(p_png, width = 8, height = 8, units = "cm", res = 300)
print(p)
trellis.focus("legend", side="right", clipp.off=TRUE, highlight=FALSE)
dev.off()

# plot sd
col <- rasterTheme(region = rev(brewer.pal(9, 'Spectral')))

p <- levelplot(s,
               par.settings = col, 
               margin = FALSE,
               #at=seq(-1.6, 1.6, length.out=101),
               #names.attr=c("Mean", "SD"),
               colorkey=list(height=0.8, width=1)) +
  layer(sp.lines(land, col="black", lwd=0.8))

# Save as png file
p_png <- paste0(out_dir, "sd_delta.png")
png(p_png, width = 8, height = 8, units = "cm", res = 300)
print(p)
trellis.focus("legend", side="right", clipp.off=TRUE, highlight=FALSE)
dev.off()



#--------------------------------------------------------------------------
# Plot monthly facets (pearson)
#--------------------------------------------------------------------------

# import rasters
tif_files <- list.files(out_dir, full.names = TRUE, recursive = TRUE, pattern = sprintf("%s_cor.tif$", jvar))
b <- stack(tif_files)
b <- crop(b, r_med)

# set color
col <- rasterTheme(region = brewer.pal(9, 'RdBu'))

# facet plot
p <- levelplot(b, layout=c(7, 2),
          par.settings = col,
          at=seq(-1, 1, length.out=101),
          names.attr=format(dates, "%b %Y"),
          colorkey=list(height=0.8, width=1)) +
  layer(sp.lines(land, col="black", lwd=0.8))


# Save as png file
p_png <- paste0(out_dir, "pearson.png")
png(p_png, width = 24, height = 8, units = "cm", res = 300)
print(p)
trellis.focus("legend", side="right", clipp.off=TRUE, highlight=FALSE)
dev.off()


## calculate coefficient of variation (CV)
u <- mean(b, na.rm=TRUE)
s <- calc(b, fun=sd, na.rm=TRUE)
cv <- s/u


