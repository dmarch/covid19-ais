#----------------------------------------------------------------------------------
# Suite of common function
#----------------------------------------------------------------------------------
# clip_to_globe



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