# summarize density maps
# calculate average and coefficient of variation of deltas during 2020

library(stringr)
library(raster)
library(lubridate)
source("scr/fun_common.R")  # bb()
source("scr/fun_ais.R")

# set input folder
input_dir <- "data/out/ais-global/sog/"

# create output directory
out_dir <- "data/out/ais-global/sog/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# select months to process
dates <- c(
  seq.Date(as.Date("2020-01-01"), as.Date("2020-06-01"), by = "month")
) 

# get year and format dates
year <- year(dates)[1]
dates <- dates %>% format("%Y%m%d")


# create empty stack
s <- stack()

# loop by month date
for(i in 1:length(dates)){
  
  idate <- dates[i]
  
  tif_file <- list.files(input_dir, full.names = TRUE, recursive=TRUE, pattern=sprintf("SOG_%s_avg.tif$", idate))
  r <- raster(tif_file)
  s <- stack(s, r)
}

# calculate average and CV
s_u <- mean(s, na.rm=TRUE)
s_sd <- calc(s, fun=sd, na.rm=TRUE)
s_cv <- s_sd/s_u

# plot average
minval <- minValue(s_u)
maxval <- maxValue(s_u)
pngfile <- paste0(out_dir, sprintf("SOG_%s_avg.png", year))
png(pngfile, width=3000, height=1750, res=300)
plotDensMol(r = s_u, zlim = c(minval, maxval), mollT = TRUE, logT = FALSE,
            col = rev(brewer.spectral(101)), main = sprintf("Average SOG (Jan-Jun %s)", year),
            axis_at = c(minval, maxval), axis_labels = c(round(minval,3), round(maxval,3)))
dev.off()

# plot coefficient of variation
minval <- minValue(s_cv)
maxval <- maxValue(s_cv)
pngfile <- paste0(out_dir, sprintf("SOG_%s_cv.png", year))
png(pngfile, width=3000, height=1750, res=300)
plotDensMol(r = s_cv, zlim = c(minval, maxval), mollT = TRUE, logT = FALSE,
            col = rev(brewer.spectral(101)), main = sprintf("Coefficient of variation SOG (Jan-Jun %s)", year),
            axis_at = c(minval, maxval), axis_labels = c(round(minval,3), round(maxval,3)))
dev.off()
 

