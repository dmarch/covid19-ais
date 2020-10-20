# choke_points


# import chokepoints
# generate buffer
# overlap with cells
# extract density and delta values at each time step.



library(sf)


# import chockepoints
choke_file <- "data/input/chokepoints/chokepoints.csv"
choke <- read.csv(choke_file)

# create spatial buffer
radius <- 20000  # in meters
choke_sf <- st_as_sf(choke, coords = c("lon", "lat"), crs = 4326)
choke_proj <- st_transform(choke_sf, crs = st_crs('+proj=moll'))
choke_buff <- st_buffer(choke_proj, radius)
choke_buff <- st_buffer(choke_sf, 0.5)

# export geopackage
st_write(choke_buff, "data/input/chokepoints/chokepoints_buff_geo.gpkg")
