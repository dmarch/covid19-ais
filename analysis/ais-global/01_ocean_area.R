#-------------------------------------------------------------------------------
# ocean_area        Estimates ocean area per cell
#-------------------------------------------------------------------------------
# AIS density maps were provided in EPSG 4326.
# We calculate ocean area per cell taking into account:
# - Difference in cell size by latitude
# - Land coverage in coastal cells



library(sf)
library(raster)
library(fasterize)

#-----------------------------------------------------------------
# Set paths
#-----------------------------------------------------------------

output_data <- "data/out/ais-global"
if (!dir.exists(output_data)) dir.create(output_data, recursive = TRUE)



#-----------------------------------------------
# 1. Import data
#-----------------------------------------------

# Import land mask
# We used GADM
shp_gadm <- "data/input/gadm/gadm36_0.shp"
sf_gadm <- st_read(shp_gadm)


#-----------------------------------------------
# 2. Calculate area accounting to latitude
#-----------------------------------------------

# create empty raster at 0.25 x 0.25 degrees (resolution of AIS density maps)
r <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90,
            crs=CRS("+proj=longlat +datum=WGS84"),
            resolution=c(0.25,0.25), vals=NULL)

# calculate area for each cell
r_area <- area(r)  # km2


#-----------------------------------------------
# 3. Calculate ocean coverage
#-----------------------------------------------
# rasterize() from raster package takes too much time
# We rather use the fasterize package. However, this function
# does not allow the calculation of coverage. Therefore, we use a two-step
# processing:
# (1) rasterize at higher resolution (x10),
# (2) aggregate at original resolution and count the number of land cells.

# create higher resolution raster
r_high <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90,
                 crs=CRS("+proj=longlat +datum=WGS84"),
                 resolution=c(0.025,0.025), vals=NULL)

# rasterize
r_land <- fasterize(sf_gadm, r_high) 

# return to original resolution and calculate land cover
r_land_cover <- aggregate(r_land, fact=10, fun=sum)
r_land_cover[is.na(r_land_cover)] <- 0  # NA cells are ocean, so change to 0



#-----------------------------------------------
# 4. Filter ocean cells based on coverage threshold
#-----------------------------------------------

# create mask based on threshold
cover_thrs <- 90
r_land_cover[r_land_cover >= cover_thrs] <- NA
r_land_cover[r_land_cover < cover_thrs] <- 1

# calculate area in ocean cells
#r_ocean_area <- r_area * ((100 - r_land_cover)/100)
#r_ocean_area[r_ocean_area == 0] <- NA

# filter cells
r_ocean_area <- r_area * r_land_cover


# transform to mollweide
r_ocean_area_moll <- projectRaster(r_ocean_area, res = 27750, crs = "+proj=moll +ellps=WGS84", method="ngb")
r_ocean_area_moll <- r_ocean_area_moll/r_ocean_area_moll

# write raster
writeRaster(r_ocean_area, "data/out/ais-global/oceanmask.nc", format="CDF", overwrite=TRUE)
writeRaster(r_area, "data/out/ais-global/areacell.nc", format="CDF", overwrite=TRUE)
writeRaster(r_ocean_area_moll, "data/out/ais-global/oceanmask_mol.nc", format="CDF", overwrite=TRUE)

