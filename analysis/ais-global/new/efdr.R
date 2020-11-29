library(EFDR)
source("scr/fun_common.R")  # bb()
source("scr/fun_ais.R")




# import delta
delta <- raster("data/out/ais-global/delta/20200401_20190401/OTHER_20200401_20190401_delta_mol.tif")

plotDelta(delta, main = "prova")



# aggregate for faster
delta <- aggregate(delta, fun=mean, fact=4)

# resample into a raster with number of cols and rows as power of two
r <- raster(nrows=128, ncols=256, ext=extent(delta), crs = crs(delta))  # 512 and 1024, 64 vs 128
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


s <- stack(r1, r2, r3)
names(s) <- c("Bonferroni","FDR" ,"LOS", "EFDR")

plot(s, col=rev(brewer.spectral(101)))


writeRaster(s, "data/test/efdr.grd", overwrite=TRUE)

s <- stack("data/test/efdr.grd")


efdr <- s$EFDR
plotDelta(s$EFDR, main = "prova")
