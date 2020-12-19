# summarize density maps
# calculate average and coefficient of variation of deltas during 2020

library(stringr)
library(raster)
library(lubridate)
source("scr/fun_common.R")  # bb()
source("scr/fun_ais.R")

# set input folder
input_dir <- "data/out/ais-global/density/"

# create output directory
out_dir <- "data/out/ais-global/density/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# select variables to compare
vars <- c("COUNT", "FISHING", "PASSENGER", "CARGO", "TANKER", "OTHER")

# select months to process
dates <- c(
  seq.Date(as.Date("2020-01-01"), as.Date("2020-06-01"), by = "month")
) 

# get year and format dates
year <- year(dates)[1]
dates <- dates %>% format("%Y%m%d")



for(j in 1:length(vars)){
  
  # select variable
  jvar <- vars[j]

  # create empty stack
  s <- stack()

  # loop by month date
  for(i in 1:length(dates)){
    
    idate <- dates[i]

    tif_file <- list.files(input_dir, full.names = TRUE, recursive=TRUE, pattern=sprintf("%s_%s_dens_moll.tif$", idate, jvar))
    r <- raster(tif_file)
    s <- stack(s, r)
  }
  
  # calculate average and CV
  s_u <- mean(s, na.rm=TRUE)
  s_sd <- calc(s, fun=sd, na.rm=TRUE)
  s_cv <- s_sd/s_u
  writeRaster(s_cv, paste0(out_dir, sprintf("%s_%s_dens_cv.tif", year, jvar)), overwrite=TRUE)
  

  # plot average density
  pngfile <- paste0(out_dir, sprintf("%s_%s_dens_avg.png", jvar, year))
  png(pngfile, width=3000, height=1750, res=300)
  plotDensMol2(r = s_u, col = rev(brewer.spectral(101)), main = sprintf("Average density %s (Jan-Jun %s)", jvar, year))
  dev.off()
  
  # plot coefficient of variation
  minval <- minValue(s_cv)
  maxval <- maxValue(s_cv)
  pngfile <- paste0(out_dir, sprintf("%s_%s_dens_cv.png", jvar, year))
  png(pngfile, width=3000, height=1750, res=300)
  plotDensMol(r = s_cv, zlim = c(minval, maxval), mollT = FALSE, logT = FALSE,
              col = rev(brewer.spectral(101)), main = sprintf("Coefficient of variation %s (Jan-Jun %s)", jvar, year),
              axis_at = c(minval, maxval), axis_labels = c(round(minval,3), round(maxval,3)))
  dev.off()
}  



# ### Compare CV between 2019 and 2020
# 
# for(j in 1:length(vars)){
# 
#   jvar <- vars[j]
#   
#   # import both CV
#   r1 <- raster(paste0(out_dir, sprintf("%s_%s_dens_cv.tif", 2019, jvar)))
#   r2 <- raster(paste0(out_dir, sprintf("%s_%s_dens_cv.tif", 2020, jvar)))
#   
#   # calculate difference
#   delta_cv <- r2 - r1
#   #per <- delta_cv/r1
#   
#   # plot
#   pngfile <- paste0(out_dir, sprintf("%s_delta_cv.png", jvar))
#   png(pngfile, width=3000, height=1750, res=300)
#   plotDelta2(delta_cv, percentile = 0.99, main = sprintf("CV 2020-2019 (%s)", jvar))
#   dev.off()
# }
