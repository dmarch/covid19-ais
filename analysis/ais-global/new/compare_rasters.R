# compare rasters
# We use maps in Mollweide projection to use same window size


library(spatialEco)


# raster.modified.ttest
# implements a moving window version of the Dutileul modified t-test,
# accounting for spatial autocorrelation (Dutilleul 1993; Clifford et al., 1989). 

# create output directory
out_dir <- "data/out/ais-global/raster_correlations/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)



# import month i densities
r1 <- raster("data/out/ais-global/density/20190401/20190401_COUNT_dens_moll.tif")
r2 <- raster("data/out/ais-global/density/20200401/20200401_COUNT_dens_moll.tif")


r1 <- crop(r1, extent(0,4e6,4e6,6e6))
r2 <- crop(r2, extent(0,4e6,4e6,6e6))

# pearson
r.cor <- rasterCorrelation(r1, r2, s=3, type="pearson")
writeRaster(r.cor, paste0(out_dir, "20200401_pearson.tif"), overwrite=TRUE)
#writeRaster(r.cor, paste0(out_dir, sprintf("%s_pearson.tif", jvar)), overwrite=TRUE)



# raster.modified.ttest
# implements a moving window version of the Dutileul modified t-test,
# accounting for spatial autocorrelation (Dutilleul 1993; Clifford et al., 1989). 

r1 <- as(r1, "SpatialPixelsDataFrame")
r2 <- as(r2, "SpatialPixelsDataFrame")

corr <- raster.modified.ttest(r1, r2, d = 27750*2)
b<- brick(corr)

writeRaster(r.cor, paste0(out_dir, "20200401_pearson.tif"), overwrite=TRUE)
