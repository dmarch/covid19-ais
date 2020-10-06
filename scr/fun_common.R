#----------------------------------------------------------------------------------
# Suite of common function
#----------------------------------------------------------------------------------
# clip_to_globe


#----------------------------------------------------------------------
# bb             Create world bounding box at custom CRS
#----------------------------------------------------------------------
bb <- function(xmin, xmax, ymin, ymax, crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"){
  # Arguments
  # xmin        Minimum longitude
  # xmax        Maximum longitude
  # ymin        Minimum latitude
  # ymax        Maximum latitude
  # crs         CRS definition
  # 
  # Details:
  # By default, the bounding box is returned in longlat.
  #
  # Description:
  # This function has been designed to create a bounding box defined by multiple points.
  # Using extent() from raster package only provides 4 points, which is not enough to plot
  # a box using other projections (eg. Mollweide).
  # We then need to create a bounding box with more density of points in order to
  # change its shape with transforming to other projections.
  #
  # Source: "https://raw.githubusercontent.com/dmarch/abigoos/master/R/utils.R"
  
  
  library(raster)
  library(rgeos)
  library(sp)
  
  # create bounding box with four points and project to mercator
  e <- as(raster::extent(xmin, xmax, ymin, ymax), "SpatialPolygons")
  proj4string(e) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  
  # transform to line and sample new points
  el <- as(e, "SpatialLines")
  len <- gLength(el)
  i <- gInterpolate(el, d=seq(0, len, by=0.5), normalized = FALSE)  # sample points at 0.5 degrees
  
  # reconvert to spatial polygon
  P1 = Polygon(i)
  Ps1 = SpatialPolygons(list(Polygons(list(P1), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  
  # reproject
  Ps1.prj <- spTransform(Ps1, crs)
  
  return(Ps1.prj)
}
#----------------------------------------------------------------------



clip_to_globe <- function(x) {
  ### Credit for this function to the OHI-Science team.
  ### https://github.com/OHI-Science/ohiprep_v2020/blob/18880e6ab7eb1607f754660b67d4b3b46f2a6043/globalprep/spp/v2020/_setup/common_fxns.R
  ### for SF features, transform to wgs84, clip to +-180 and +-90
  epsg <- st_crs(x)$epsg
  if(epsg != 4326 | is.na(epsg)) {
    message('Original EPSG = ', epsg, '; Proj4 = ', st_crs(x)$proj4string,
            '\n...converting to EPSG:4326 WGS84 for clipping')
    x <- st_transform(x, 4326)
  }
  x_bbox <- st_bbox(x)
  if(x_bbox$xmin < -180 | x_bbox$xmax > +180 |
     x_bbox$ymin <  -90 | x_bbox$ymax >  +90) {
    message('Some bounds outside +-180 and +-90 - clipping')
    z <- st_crop(x, y = c('xmin' = -180,
                          'ymin' =  -90,
                          'xmax' = +180,
                          'ymax' =  +90)) %>%
      st_cast('MULTIPOLYGON')
  } else {
    message('All bounds OK, no clipping necessary')
    z <- x
  }
  return(z)
}