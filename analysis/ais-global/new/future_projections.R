# Before continue, check similar approaches.


# projections no covid
# We assume that percentage variation in 2020 
# Use January as baseline
# Calculate percentage variation in 2019 Feb-June from January 2019
# Project such variation in 2020 Feb-June using January 2020. This will create a projected density.
# Calculate delta between observed density and projected density
# Use zonal 3 x 3


# For February 202
# import January 2020 density
# import percentage Feb2019-Jan2019
# multiply


# import January densities
jan19 <- raster("data/out/ais-global/density/20190101/20190101_COUNT_dens_moll.tif")
jan20 <- raster("data/out/ais-global/density/20200101/20200101_COUNT_dens_moll.tif")

# import month i densities
m19 <- raster("data/out/ais-global/density/20190401/20190401_COUNT_dens_moll.tif")
m20 <- raster("data/out/ais-global/density/20200401/20200401_COUNT_dens_moll.tif")

# focal
jan19f <- focal(jan19, w=matrix(1,5,5), fun=mean, na.rm=TRUE)
jan20f <- focal(jan20, w=matrix(1,5,5), fun=mean, na.rm=TRUE)
m19f <- focal(m19, w=matrix(1,5,5), fun=mean, na.rm=TRUE)
m20f <- focal(m20, w=matrix(1,5,5), fun=mean, na.rm=TRUE)


# calculate projection for month i in 2020
mproj <- (jan20 * m19)/jan19

# create common mask between both rasters
# identifies cells with presence of route density in any of the two periods
m <- sum(m20, mproj, na.rm=TRUE)
m[m==0] <- NA
m[m>0] <- 0

# add zero to density maps
r1z <- sum(mproj, m, na.rm=TRUE)
r2z <- sum(m20, m, na.rm=TRUE)

# calculate delta
delta <- sum(r2z, r1z*(-1), na.rm=TRUE)
delta <- mask(delta, m)


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
pngfile <- paste0(out_dir_month, sprintf("%s_%s_%s_delta_mol.png", jvar, idate_post, idate_pre))
png(pngfile, width=3000, height=1750, res=300)
plotDensMol(r = delta, zlim = c(qmin, qmax), mollT = FALSE, logT = FALSE, breaks=breaks,
            col = cols, #main = sprintf("Delta %s (%s vs %s)", jvar, idate_post, idate_pre),
            axis_at = c(qmin, 0,qmax),
            axis_labels = c(paste("<",round(qmin,2)), 0, paste(">",round(qmax,2))),
            legend_horizontal=TRUE)
dev.off()













asinh_trans <- function(){
  trans_new(name = 'asinh', transform = function(x) asinh(x), 
            inverse = function(x) sinh(x))
}

asinh(0.05)
            

delta_trans3 <- calc(delta, fun=function(x) asinh(exp(5)*x))

delta_trans1 <- calc(delta, fun=function(x) sign(x)*log10(abs(x)))
delta_trans2 <- delta_trans1/maxValue(delta_trans1)
plot(delta_trans2)

trans <- function(r){
  l = sign(r)*log1p(abs(r))
  l = l/max(l)
  return(l)
}
https://stats.stackexchange.com/questions/1444/how-should-i-transform-non-negative-data-including-zeros




log10(10) - log10(5)
10^0.30103


r = -1000:1000

l = sign(r)*log1p(abs(r))
l = l/max(l)
plot(r, l, type = "l", xlab = "Original", ylab = "Transformed", col = adjustcolor("red", alpha = 0.5), lwd = 3)





#We scale both to fit (-1,1)
for(i in exp(seq(-10, 100, 10))){
  s = asinh(i*r)
  
  s = s / max(s)
  lines(r, s, col = adjustcolor("blue", alpha = 0.2), lwd = 3)
}
legend("topleft", c("asinh(x)", "sign(x) log(abs(x)+1)"), col = c("blue", "red"), lty = 1)