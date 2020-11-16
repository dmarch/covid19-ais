#-------------------------------------------------------------------------------
# ocean_area        Estimates ocean area per cell
#-------------------------------------------------------------------------------
# AIS density maps were provided in EPSG 4326.
# We calculate ocean area per cell taking into account:
# - Difference in cell size by latitude
# - Land coverage in coastal cells
#
# Output:
# "ocean_area.tif": raster layer in lon/lat representing the ocean area for each cell.



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
#shp_gadm <- "data/input/gadm/gadm36_0.shp"
#sf_gadm <- st_read(shp_gadm)

# Import GSHHS
sf_l1 <- st_read("data/input/gshhg-shp-2.3.7/GSHHS_shp/f/GSHHS_f_L1.shp")
sf_l6 <- st_read("data/input/gshhg-shp-2.3.7/GSHHS_shp/f/GSHHS_f_L6.shp")


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
# 3. Calculate ocean coverage on each cell
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
#r_land <- fasterize(sf_gadm, r_high) 

r_land_l1 <- fasterize(sf_l1, r_high) 
r_land_l6 <- fasterize(sf_l6, r_high)
r_land <- sum(r_land_l1, r_land_l6, na.rm=TRUE)
#r_land[r_land==0] <- NA

# return to original resolution and calculate land cover
r_land_cover <- aggregate(r_land, fact=10, fun=sum)
#r_land_cover[is.na(r_land_cover)] <- 0  # NA cells are ocean, so change to 0
writeRaster(r_land_cover, "data/out/ais-global/land_cover.tif", overwrite=TRUE)



#-----------------------------------------------
# 4. Filter ocean cells based on coverage threshold
#-----------------------------------------------

# create mask based on threshold (v1 use this part)
cover_thrs <- 95
r_land_cover[r_land_cover > cover_thrs] <- NA
#r_land_cover[r_land_cover < cover_thrs] <- 1

# calculate area in ocean cells (v1 commented)
r_ocean_area <- r_area * ((100 - r_land_cover)/100)
#r_ocean_area[r_ocean_area == 0] <- NA

# filter cells  (v1 use this part)
#r_ocean_area <- r_area * r_land_cover


# transform to mollweide
#r_ocean_area_moll <- projectRaster(r_ocean_area, res = 27750, crs = "+proj=moll +ellps=WGS84", method="ngb")
#r_ocean_area_moll <- r_ocean_area_moll/r_ocean_area_moll  (v1 use this part)

# write raster
writeRaster(r_ocean_area, "data/out/ais-global/ocean_area.tif", overwrite=TRUE)
#writeRaster(r_area, "data/out/ais-global/areacell.nc", format="CDF", overwrite=TRUE)
#writeRaster(r_ocean_area_moll, "data/out/ais-global/oceanmask_mol.nc", format="CDF", overwrite=TRUE)



#------------------------------------------------
# 5. Create distance to coast raster
#------------------------------------------------
# Note this step is time consuming

# # distance to shore
# r_ocean_area[is.na(r_ocean_area)] <- 999999
# r_ocean_area[r_ocean_area < 1000] <- NA
# dist2coast <- distance(r_ocean_area)
# 
# # save raster
# writeRaster(dist2coast, "data/out/ais-global/dist2coast.tif")

