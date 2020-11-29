#-------------------------------------------------------------------------
# 04_sector_contribution
#-------------------------------------------------------------------------


## load libraries
library(raster)
source("https://raw.githubusercontent.com/dmarch/abigoos/master/R/utils.R")  # bb()
source("scr/fun_ais.R")

# set input dir
input_dir <- "data/out/ais-global/delta/"

# create output directory
out_dir <- "data/out/ais-global/contrib/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# select variables to summarize
vars <- c("FISHING", "PASSENGER", "CARGO", "TANKER", "OTHER")



# import files
tif_files <- paste0(input_dir, sprintf("%s_delta_mol.tif", vars))
s <- stack(tif_files)

# calculate % contribution of each sector (delta values)
sabs <- abs(s)  # calculate absolute values
ssum <- sum(sabs, na.rm=TRUE)  # sum values
sprop <- 100 * (sabs/ssum)
names(sprop) <- vars

# for each layer, save raster and plot
for(i in 1:nlayers(sprop)){
  
  # extract layer
  r <- subset(sprop, i)
  ivar <- names(r)
  writeRaster(r, paste0(out_dir, sprintf("%s_contrib.tif", ivar)), overwrite=TRUE)
  
  # plot
  minval <- 0
  maxval <- 100
  pngfile <- paste0(out_dir, sprintf("%s_contrib.png", ivar))
  png(pngfile, width=3000, height=1750, res=300)
  plotDensMol(r = r, zlim = c(minval, maxval), mollT = FALSE, logT = FALSE,
              col = rev(brewer.spectral(101)), main = sprintf("%s Contribution (Percentage)", ivar),
              axis_at = c(minval, maxval), axis_labels = c(minval, maxval))
  dev.off()
}