


library(elsa)
#> Loading required package: sp
#> Loading required package: raster

file <- system.file('external/dem_example.grd',package='elsa') 

r <- raster(file) # reading a raster map (Dogital Elevation Model: DEM)

plot(r, main='DEM: a continuous raster map')




# ELSA statistic for the local distance of 2Km:
e <- elsa(r,d=2000,categorical=FALSE) 

cl <- colorRampPalette(c('darkblue','yellow','red','black'))(100) # specifying a color scheme

# Following is the map of ELSA, lower values represent higher local spatial autocorrelation
plot(e,col=cl, main='ELSA for the local distance of 2 Km')


# ELSA statistic for the local distance of 4Km:
e2 <- elsa(r,d=4000,categorical=FALSE) 

plot(e2,col=cl, main='ELSA for the local distance of 4 Km')





en <- entrogram(r, width = 2000, cutoff = 15000)
#> the input is considered as a continuous variable...

plot(en)



v <- Variogram(r, width = 2000, cutoff = 1400000)

plot(v)

moran(r, d1=0, d2=2000)


co <- correlogram(r, width=2000,cutoff=80000)

plot(co)




et <- elsa.test(r,d=2000,n=99, categorical=FALSE)

plot(et)




### Marine Traffic
r <- raster("data/out/ais-global/density/20190101/20190101_COUNT_dens_moll.tif")
r <- log10(r)
e <- elsa(r,d=555000,categorical=FALSE) 
plot(e,col=cl)

e2 <- elsa(r,d=277500,categorical=FALSE) 
plot(e2,col=cl)


et <- elsa.test(r,d=277500,n=99, categorical=FALSE)


v <- Variogram(r, width = 55500, cutoff = 1000000)

r2 <- log10(r)
co2 <- correlogram(r2, width = 55500, cutoff = 10000000)
plot(co)


en <- entrogram(r, width = 55500, cutoff = 1000000)
plot(en)





## Not run: 
library(gstat)                                         
library(sp)                                            

data(meuse)                                            
data(meuse.grid)                                       
coordinates(meuse) <- ~x + y                           
coordinates(meuse.grid) <- ~x + y                      

# GRID-1 log(copper):                                              
v1 <- variogram(log(copper) ~ 1, meuse)                  
x1 <- fit.variogram(v1, vgm(1, "Sph", 800, 1))           
G1 <- krige(zinc ~ 1, meuse, meuse.grid, x1, nmax = 30)
gridded(G1) <- TRUE                                      
G1@data = as.data.frame(G1@data[,-2])

# GRID-2 log(elev):                                              
v2 <- variogram(log(elev) ~ 1, meuse)                  
x2 <- fit.variogram(v2, vgm(.1, "Sph", 1000, .6))        
G2 <- krige(elev ~ 1, meuse, meuse.grid, x2, nmax = 30)
gridded(G2) <- TRUE    
G2@data <- as.data.frame(G2@data[,-2])
G2@data[,1] <- G2@data[,1]

corr <- raster.modified.ttest(G1, G2)
plot(raster::raster(corr,1))

corr.rand <- raster.modifed.ttest(G1, G2, sub.sample = TRUE, type = "random")	 
corr.hex <- raster.modifed.ttest(G1, G2, sub.sample = TRUE, d = 500, size = 1000)	
head(corr.hex@data)
bubble(corr.hex, "corr") 



s.pix1 <- as(raster(xmn=0,xmx=10,ymn=0,ymx=10,res=1,vals=sample(1:5,100,replace=T)), "SpatialGridDataFrame")
s.pix1[sample(1:100,25)]<-NA
s.pix2 <- as(raster(xmn=0,xmx=10,ymn=0,ymx=10,res=1,vals=sample(1:5,100,replace=T)), "SpatialGridDataFrame")
s.pix2[sample(1:100,25)]<-NA

corr <- spatialEco::raster.modifed.ttest(x=s.pix1, y=s.pix2)

plot(raster(corr), main="spatially adjusted raster correlation")





r1 <- raster("data/out/ais-global/density/20190401/20190401_COUNT_dens_moll.tif")
r2 <- raster("data/out/ais-global/density/20200401/20200401_COUNT_dens_moll.tif")

r1 <- as(r1, "SpatialPixelsDataFrame")
r2 <- as(r2, "SpatialPixelsDataFrame")


corr <- raster.modified.ttest(r1, r2, d=27750*2)
b<- brick(corr)

corr.rand <- raster.modified.ttest(r1, r2, sub.sample = TRUE, d = 27750*2, size = 1000)
