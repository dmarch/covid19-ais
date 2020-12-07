# compare rasters
# We use maps in Mollweide projection to use same window size


library(spatialEco)
library(elsa)

# raster.modified.ttest
# implements a moving window version of the Dutileul modified t-test,
# accounting for spatial autocorrelation (Dutilleul 1993; Clifford et al., 1989). 

# create output directory
out_dir <- "data/out/ais-global/raster_correlations/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)



# import month i densities
r1 <- raster("data/out/ais-global/density/20190401/20190401_COUNT_dens_moll.tif")
r2 <- raster("data/out/ais-global/density/20200401/20200401_COUNT_dens_moll.tif")


r1 <- crop(r1, extent(-1e6,4e6,3e6,6e6))
r2 <- crop(r2, extent(-1e6,4e6,3e6,6e6))
delta <- r2-r1
per <- delta/r1

# pearson
r.cor <- rasterCorrelation(r1, r2, s=3, type="pearson")
writeRaster(r.cor, paste0(out_dir, "20200401_pearson.tif"), overwrite=TRUE)
#writeRaster(r.cor, paste0(out_dir, sprintf("%s_pearson.tif", jvar)), overwrite=TRUE)



# elsa
cl <- colorRampPalette(c('darkblue','yellow','red','black'))(100) # specifying a color scheme
e <- elsa(delta,d=27750*10,categorical=FALSE) 
plot(e, col=cl)


et <- elsa.test(delta,d=27750*10,n=99, categorical=FALSE)
en <- entrogram(delta, width = 2000, cutoff = 27750*10)
#> the input is considered as a continuous variable...

plot(en)

co <- correlogram(delta, width=2000,cutoff=27750*100)
plot(co)

v <- Variogram(delta, width = 27750, cutoff = 27750*15)
plot(v)


lisa.g <- lisa(delta, d1=0,d2=27750*5,statistic='g*')
plot(lisa.g, col=cl)

lisa.i <- lisa(r2, d1=0,d2=27750*5,statistic='I')
plot(lisa.i,col=cl,main="Local Moran's I")


# raster.modified.ttest
# implements a moving window version of the Dutileul modified t-test,
# accounting for spatial autocorrelation (Dutilleul 1993; Clifford et al., 1989). 

# create common mask
r1mask <- r1/r1
r2mask <- r2/r2
rmask <- r1mask * r2mask
r1 <- mask(r1, rmask)
r2 <- mask(r2, rmask)

r1 <- as(r1, "SpatialPixelsDataFrame")
r2 <- as(r2, "SpatialPixelsDataFrame")

corr <- raster.modified.ttest(r1, r2, d = 27750*2)
b<- brick(corr)


#---------------------------------------------------------------------
# Mann-Kendall trend
#---------------------------------------------------------------------
require(raster)
require(Kendall)
library(devtools)
source_gist("https://gist.github.com/hakimabdi/8ca1d2783b0ab8aec5cac993573b08ee")

# smooth
pre_mol_s <- smooth.time.series(pre_mol, f = 0.8, smooth.data = FALSE)

s<-stack(r1, r2)
trend <- MKraster(s, type = "both")
plot(trend)


x <- corLocal(r2, r1, ngb=5, method=c("spearman"), test=TRUE )



#---------------------------------------------------------------------
# EFDR
#---------------------------------------------------------------------
# import delta
library(EFDR)

delta <- raster("data/out/ais-global/delta/20200401_20190401/COUNT_20200401_20190401_delta_mol.tif")

delta <- crop(delta, extent(0,4e6,4e6,6e6))



# resample into a raster with number of cols and rows as power of two
r <- raster(nrows=64, ncols=128, ext=extent(delta), crs = crs(delta))  # 512 and 1024, 64 vs 128
delta <- resample(delta, r)


# convert to matrix
Z <- as.matrix(delta)
Z[is.na(Z)] <- 0  # NA are not accepted


# Set parameters
alpha =     0.05        # 5% significant level
wf =        "la8"       # wavelet name
J =         3           # 3 resolutions
b =         11          # 11 neighbours for EFDR tests
n.hyp =     c(5,10,20,30,40,60,80,100,400,800,1200,2400,4096)  # find optimal number of tests in 
# EFDR from this list
iteration = 100         # number of iterations in MC run when finding optimal n in n.hyp
idp =       0.5         # inverse distance weighting power
nmax =      15          # maximum number of points to consider when carrying out inverse 
# distance interpolation
parallel =  parallel::detectCores()/2 # use half the number of available cores
set.seed(20)            # same random seed for reproducibility


# Conduct tests
m1 <- test.bonferroni(Z, wf=wf,J=J, alpha = alpha)
m2 <- test.fdr(Z, wf=wf,J=J, alpha = alpha)
m3 <- test.los(Z, wf=wf,J=J, alpha = alpha)
m4 <- test.efdr(Z, wf="la8",J=J, alpha = alpha, n.hyp = n.hyp, 
                b=b,iteration=iteration,parallel = parallel)


# Convert to raster
r1 <- raster(m1$Z, crs=crs(delta))
extent(r1) <- extent(delta)
r1 <- mask(r1, delta)

r2 <- raster(m2$Z, crs=crs(delta))
extent(r2) <- extent(delta)
r2 <- mask(r2, delta)

r3 <- raster(m3$Z, crs=crs(delta))
extent(r3) <- extent(delta)
r3 <- mask(r3, delta)

r4 <- raster(m4$Z, crs=crs(delta))
extent(r4) <- extent(delta)
r4 <- mask(r4, delta)


#----------------------------
library(diffeR)
d <- differenceMR(comp = r2, ref = r1)





# NOT RUN {
old.par <- par(no.readonly = TRUE)
grid1 <- raster(system.file("external/GRID1_INT.rst", package="diffeR"))
grid2 <- raster(system.file("external/GRID2_INT.rst", package="diffeR"))
strata <- raster(system.file("external/strata_int.rst", package="diffeR"))
MAD(grid1, grid2, strata, eval="original")
MAD(grid1, grid2, strata, eval="multiple")

# }
# NOT RUN {
veg_obs1 <- raster(system.file("external/veg_obs1.rst", package="diffeR")) 
veg_pre1 <- raster(system.file("external/veg_pre1.rst", package="diffeR")) 
veg_mask1 <- raster(system.file("external/veg_mask1.rst", package="diffeR")) 
MADscatterplot(veg_obs1, veg_pre1, veg_mask1) 
MAD(veg_obs1, veg_pre1, veg_mask1, eval="multiple") 
# }
# NOT RUN {
par(old.par)
# }



