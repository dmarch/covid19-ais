#---------------------------------------------------------------------------
# plot zooms delta          Plot zooms for delta map at regional level
#---------------------------------------------------------------------------
# This script proceses density maps


library(raster)
library(fasterize)
library(sf)
library(pals)
library(rasterVis)
library(pals)
library(RColorBrewer)
source("scr/fun_common.R")  # bb()
source("scr/fun_ais.R")


# set input directory
input_dir <- "data/out/ais-global/delta/"

# create output directory
out_dir <- "results/ais-global/zoom_regions/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# select variables to compare
vars <- c("COUNT")


# select months to process

### Option 1: select months to process (2020 vs 2019)
dates_post <- c(
  seq.Date(as.Date("2020-01-01"), as.Date("2020-06-01"), by = "month")
)

dates_pre <- c(
  seq.Date(as.Date("2019-01-01"), as.Date("2019-06-01"), by = "month")
)


# subsets
ext <- list(eu = extent(-11, 14, 35, 53), #extent(-15, 25, 35, 60),
            gal = extent(-93, -87, -2.8, 1.8)+1,#extent(-92, -78, -3, 11),#extent(-93, -87, -2.8, 1.8),#extent(-100, -80, -8, 8),
            van = extent(-130, -120, 46, 52.5),#extent(-140, -120, 43, 53),
            ch =  extent(116, 135, 26, 41),#extent(110, 145, 20, 45),
            sue = extent(30, 61, 8, 35),
            mal = extent(93, 115, -8, 13),
            ind = extent(110, 145, -22, 8),
            ara = extent(40, 83, -3, 32))


#-------------------------------
# 3. Prepare stack
#-------------------------------


s <- stack()

for (i in 1:length(dates_post)){
  for(j in 1:length(vars)){
  
  print(paste("Processing month", i, "from", length(dates_post)))
  
  # get month
  idate_post <- dates_post[i] %>% format("%Y%m%d")
  idate_pre <- dates_pre[i] %>% format("%Y%m%d")
  
  # set input and output folder for selected month
  input_dir_month <- paste0(input_dir, idate_post, "_", idate_pre, "/")
  
  # get var
  jvar <- vars[j]
    
  # import delta
  delta_file <- paste0(input_dir_month, sprintf("%s_%s_%s_delta.tif", jvar, idate_post, idate_pre))
  delta <- raster(delta_file)
    
  # append to stack
  s <- stack(s, delta)
  }  
}
    


#-------------------------------
# 4. Plot
#-------------------------------

# import land mask
world <- st_read("data/input/gshhg-shp-2.3.7/GSHHS_shp/l/GSHHS_l_L1.shp") 
world <- as(world, Class = "Spatial")
#gadm <- st_read("data/input/gadm/gadm36_0.shp")
#world <- as(gadm, Class = "Spatial")

# set color
#col <- rasterTheme(region = brewer.pal(9, 'RdBu'))


for(i in 1:length(ext)){
  
  # get extent
  e_name <- names(ext)[i]
  e <- ext[[i]]
  
  # crop land
  land <- crop(world, e)
  
  # crop map
  delta <- crop(s, e)

  # get min and max values
  # mindelta <- min(minValue(delta))
  # maxdelta <- max(maxValue(delta))
  # maxval <- max(abs(mindelta), maxdelta)
  # minval <- maxval * (-1)
  
  # reclass values below or above give threshold
  maxval <- 0.2
  minval <- -0.2
  delta[delta > maxval] <- maxval
  delta[delta < minval] <- minval
  
  # color key
  my.at <- seq(minval, maxval, length.out=101)
  my.labels <- c(minval, 0, maxval)
  myColorkey <- list(space="bottom", height=0.4, width=1.2,
                     at=my.at, ## where the colors change
                     labels=list(
                       at=my.labels ## where to print labels
                     ))
  
  # cols
  cols <- colorRampPalette(brewer.pal(9,"RdBu"))(length(my.at) - 1)
  
  # facet plot
  p <- levelplot(delta, layout=c(6, 1),
                 xlab=NULL, ylab=NULL,
                 maxpixels = 1e6,
                 #par.settings = col,
                 col.regions = cols,
                 at = my.at,
                 names.attr=format(dates_post, "%b"),
                 colorkey=myColorkey) + # list(height=0.8, width=1)
    layer(sp.lines(land, col="grey10", lwd=1, alpha=1))
  
  # Save as png file
  p_png <- paste0(out_dir, sprintf("%s_delta_zoom.png", names(ext)[i]))
  png(p_png, width = 25, height = 8, units = "cm", res = 300)
  print(p)
  trellis.focus("legend", side="bottom", clipp.off=TRUE, highlight=FALSE)
  dev.off()
}


# additional plot for legend

myColorkey <- list(space="bottom", height=0.4, width=1.2,
                   at=my.at, ## where the colors change
                   labels=list(
                     at=my.labels ## where to print labels
                   ))

p <- levelplot(delta, xlab=NULL, layout=c(6, 1),
               maxpixels = 1e6,
               par.settings = col,
               names.attr=format(dates_post, "%b"),
               colorkey=myColorkey)

# Save as png file
p_png <- paste0(out_dir, "legend_delta_zoom.png")
png(p_png, width = 26, height = 8, units = "cm", res = 300)
print(p)
trellis.focus("legend", side="bottom", clipp.off=TRUE, highlight=FALSE)
dev.off()





#--------------------------------------------------------------
# Plot same style with average of delta
#--------------------------------------------------------------

# import delta
jvar <- vars[1]
s <- raster(paste0("data/out/ais-global/delta_summary/", sprintf("%s_delta_u.tif", jvar)))
s <- projectRaster(s, crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# subsets
ext <- list(eu = extent(-11, 14, 35, 53), #extent(-15, 25, 35, 60),
            gal = extent(-93, -87, -2.8, 1.8)+2,#extent(-92, -78, -3, 11),#extent(-93, -87, -2.8, 1.8),#extent(-100, -80, -8, 8),
            van = extent(-130, -120, 46, 52.5),#extent(-140, -120, 43, 53),
            ch =  extent(116, 135, 26, 41),#extent(110, 145, 20, 45),
            sue = extent(30, 61, 8, 35),
            mal = extent(93, 115, -8, 13),
            ind = extent(110, 145, -22, 8),
            ara = extent(40, 83, -3, 32))


for(i in 1:length(ext)){
  
  # get extent
  e_name <- names(ext)[i]
  e <- ext[[i]]
  
  # crop land
  land <- crop(world, e)
  
  # crop map
  delta <- crop(s, e)
  
  # get min and max values
  # mindelta <- min(minValue(delta))
  # maxdelta <- max(maxValue(delta))
  # maxval <- max(abs(mindelta), maxdelta)
  # minval <- maxval * (-1)
  
  # reclass values below or above give threshold
  maxval <- 0.2
  minval <- -0.2
  delta[delta > maxval] <- maxval
  delta[delta < minval] <- minval
  
  # color key
  my.at <- seq(minval, maxval, length.out=101)
  my.labels <- c(minval, 0, maxval)
  myColorkey <- list(space="bottom", height=0.4, width=1.2,
                     at=my.at, ## where the colors change
                     labels=list(
                       at=my.labels ## where to print labels
                     ))
  
  
  # cols
  cols <- colorRampPalette(brewer.pal(9,"RdBu"))(length(my.at) - 1)
  
  # facet plot
  p <- levelplot(delta, #layout=c(6, 1),
                 xlab=NULL, ylab=NULL,
                 maxpixels = 1e6,
                 #par.settings = col,
                 col.regions = cols,
                 at = my.at,
                 #names.attr=format(dates_post, "%b"),
                 colorkey=myColorkey) + # list(height=0.8, width=1)
    layer(sp.lines(land, col="grey10", lwd=1, alpha=1))
  
  # Save as png file
  p_png <- paste0(out_dir, sprintf("%s_delta_zoom.png", names(ext)[i]))
  png(p_png, width = 25, height = 8, units = "cm", res = 300)
  print(p)
  trellis.focus("legend", side="bottom", clipp.off=TRUE, highlight=FALSE)
  dev.off()
}






