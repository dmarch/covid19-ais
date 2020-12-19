#-------------------------------------------------------------------------------------------------
# temporal_delta.R        Assess temporal variability of deltas across multiple months
#-------------------------------------------------------------------------------------------------
# Calculate sum across all period
# Calculate average and SD
# Find month with higher decrease

source("scr/fun_common.R")  # bb()
source("scr/fun_ais.R")
library(ggplot2)
library(ggridges)
library(egg)
library(dplyr)
library(raster)
library(pals)
library(reshape2)
library(exactextractr)
library(tmap)
library(doParallel)
library(foreach)
library(parallel)
library(exactextractr)
library(sf)
library(raster)
library(tmap)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(xlsx)
library(egg)
library(reshape2)
library(flextable)
library(officer)
library(fasterize)

# create output directory
input_dir <- "data/out/ais-global/delta_summary/"

# select variables to compare
vars <- c("COUNT", "FISHING", "PASSENGER", "CARGO", "TANKER", "OTHER")


#----------------------------------------------------
# Part 1. Prepare EEZ and combination table
#----------------------------------------------------

# import EEZ
# source: marineregions
eez <- st_read("data/input/marine_regions/World_EEZ_v11_20191118_gpkg/eez_v11_covid.gpkg")
poly_moll <- eez %>%
  filter(!POL_TYPE %in% c("Joint regime")) %>% # remove joint regimes
  clip_to_globe() %>%  # ensures no coordinates outside +-180, +-90
  st_transform("+proj=moll")


# Prepare cluster for parallel computing
cl <- makeCluster(6)
registerDoParallel(cl)

# Extract data for each combination  
data <- bind_rows(foreach(j = 1:length(vars), .packages=c("dplyr", "raster", "exactextractr")) %dopar% { 
  
#for (j in 1:length(vars)){
  
  # select variable
  jvar <- vars[j]

  # import data
  rast <- raster(paste0(input_dir, sprintf("%s_delta_u.tif", jvar)))
 
  # summarize data per polygon
  # weighted by the fraction of each cell that is covered by the polygon
  extr <- exact_extract(rast, poly_moll, "mean")
  
  # generate data.frame
  df <- data.frame(MRGID = poly_moll$MRGID, var = jvar, measure = "delta", value = extr)
  df
})

stopCluster(cl)  # Stop cluster





#----------------------------------------------------
# Part 5. Plot map
#----------------------------------------------------

# combine with polygon map
poly_moll$MRGID <- as.character(poly_moll$MRGID)
data$MRGID <- as.character(data$MRGID)
poly_ais <- left_join(poly_moll, data, by = c("MRGID" = "MRGID"))

# transform to data.frame
df <- data.frame(poly_ais)
df <- dplyr::select(df, -geom)

# summarize
df %>%
  group_by(var) %>%
  summarize(negative_eez = length(which(value < 0)),
            positive_eez = length(which(value > 0)),
            no_change = length(which(value == 0)),
            total = length(value),
            proportion_neg = 100 * (negative_eez / total),
            proportion_pos = 100 * (positive_eez / total),
            proportion_nochange = 100 * (no_change / total))

length(unique(df$MRGID))  # number of EEZ: 255





# Export table
#outfile <- "data/out/ais-global/eez_change.csv"
#write.csv(df, outfile, row.names = FALSE)

# import landmask
data(countriesHigh, package = "rworldxtra", envir = environment())
countriesHigh <- spTransform(countriesHigh, CRS("+proj=moll +ellps=WGS84"))
land <- st_as_sf(countriesHigh)

# create box
box <- bb(xmin = -180, xmax = 180, ymin = -90, ymax = 90, crs="+proj=moll +ellps=WGS84")

# common scale to compare with other regions
breaks <- c(-0.8, -0.4, -0.1, -0.05, -0.02, -0.01, -0.005, 0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.4)
labels <- c("-0.8 to -0.4", "-0.4 to -0.1", "-0.1 to -0.05", "-0.05 to -0.02", "-0.02 to -0.01", "-0.01 to -0.005", "-0.005 to 0.00",
            "0.00 to 0.005", "0.005 to 0.01", "0.01 to 0.02", "0.02 to 0.05", "0.05 to 0.1", "0.1 to 0.4")

for(i in 1:length(vars)){
  
  # select variable
  ivar <- vars[i]
  
  # plot variable
  p1 <- tm_shape(filter(poly_ais, var == ivar), projection="+proj=moll +ellps=WGS84") +
    tm_polygons("delta",
                style = "fixed",
                breaks = breaks,
                labels = labels,
                midpoint = 0,
                palette = "RdBu",
                border.col = "grey60",
                border.alpha = 0.3,
                textNA = "No data",
                colorNA = "white") +
    tm_shape(land, projection="+proj=moll +ellps=WGS84") +
    tm_polygons(col="grey60", border.alpha = 0) +
    tm_shape(box) +
    tm_polygons(alpha=0, border.col = "grey60") +
    tm_layout(frame = FALSE, legend.title.size=1, legend.outside = TRUE)
  
  
  # export plot
  out_file <- paste0(input_dir, sprintf("%s_eez_delta.png", ivar))
  tmap_save(tm = p1, filename = out_file, width=25, height=12, units = "cm")
}

