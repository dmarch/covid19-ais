


library(raster)
library(ncdf4)
library(RColorBrewer)

cols = rev(colorRampPalette(brewer.pal(11, 'Spectral'))(255)) # rainbow color scheme


ncin <- nc_open("data/input/tournadre/shiptarfile1992-2020/prod_gridded_1992-2020.nc")
ncin <- nc_open("data/input/tournadre/shiptarfile1992-2020/ships.nc")

#ncin <- nc_open(file.path(dir_M, "git-annex/impact_acceleration/stressors/shipping/ships.nc"))
print(ncin)
attributes(ncin$var)$names
nc_close(ncin)




## get the data and create a raster stack
raw <- stack("data/input/tournadre/shiptarfile1992-2020/prod_gridded_1992-2020.nc",
             varname="nships_smoothed") #nships_smoothed


plot_data <- rotate(flip(flip(t((raw[[25]])), direction = 'y'), direction = 'x'))
plot(plot_data, col=cols)
#click(tmp)
maps::map('world', col='gray95', fill=T, border='gray80', add=T)



# convert names to year (equal to days since 1992)
names(raw) <- round(as.numeric(sub("X", "", names(raw)))/365.2422)
## convert 0 values to NA (since not clear for each raster whether it is NA or zero)
for (year in 1994:2016){ #year = 2015
  rast_year <- grep(year, names(raw))
  tmp <- raw[[rast_year]]
  tmp <- rotate(flip(flip(t((tmp)), direction = 'y'), direction = 'x'))
  tmp[tmp==0] <- NA 
  writeRaster(tmp, file.path(dir_M, sprintf("git-annex/impact_acceleration/stressors/shipping/int/shipping_raw_%s.tif", year)),
              overwrite=TRUE)
}