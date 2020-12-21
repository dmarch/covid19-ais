#------------------------------------------------------------------------------
# extract_dens_by_EEZ
#------------------------------------------------------------------------------
# This script extract vessel density information per EEZ for each variable and timing
# The output is a data.frame where each row represent a grid cell assigned to a EEZ


## libraries
library(raster)
library(dplyr)
library(sf)
library(fasterize)
library(data.table)
library(foreach)
library(parallel)
library(doParallel)
source("scr/fun_common.R")
source("scr/fun_table.R")
source("scr/fun_ais.R")


## Prepare clusters
cores <- 10 # detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)


#----------------------------------------------------
# 1. Set input parameters
#----------------------------------------------------

# set input directory
input_dir <- "data/out/ais-global/delta/"

# output directory
out_dir <- "data/out/ais-global/eez/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# select variables to summarize
vars <- c("COUNT", "FISHING", "PASSENGER", "CARGO", "TANKER", "OTHER")

# set dates
dates_post <- c(
  seq.Date(as.Date("2020-01-01"), as.Date("2020-06-01"), by = "month")
) %>% format("%Y%m%d")

dates_pre <- c(
  seq.Date(as.Date("2019-01-01"), as.Date("2019-06-01"), by = "month")
) %>% format("%Y%m%d")



#-----------------------------------------------
# 2. Rasterize EEZ
#-----------------------------------------------

#create empty raster at 0.25 x 0.25 degrees (resolution of AIS density maps)
r <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90,
            crs=CRS("+proj=longlat +datum=WGS84"),
            resolution=c(0.25,0.25), vals=NULL)


# create higher resolution raster
r_high <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90,
                 crs=CRS("+proj=longlat +datum=WGS84"),
                 resolution=c(0.025,0.025), vals=NULL)

# import EEZ
# source: marineregions
eez <- st_read("data/input/marine_regions/World_EEZ_v11_20191118_gpkg/eez_v11_covid.gpkg")
eez <- eez %>%
  mutate(TERRITORY1 = as.character(TERRITORY1),
         SOVEREIGN1 = as.character(SOVEREIGN1)) %>%
  filter(!POL_TYPE %in% c("Joint regime", "Overlapping claim"), # remove joint regimes
         TERRITORY1!="Antarctica",  # remove Antarctica
         TERRITORY1 == SOVEREIGN1,  # select main territories (excludes overseas)
         AREA_KM2>(769*3)) # remove EEZ smaller than 3 grid sizes at equator

# rasterize at high resolution
mrgid <- fasterize(eez, r_high, "MRGID")
#writeRaster(mrgid, "data/out/ais-global/eez/mrgid_high.tif")

# return to original resolution and calculate surface
r_eez <- resample(mrgid, r, method="ngb")
#writeRaster(r_eez, "data/out/ais-global/eez/mrgid_resamp.tif")


#----------------------------------------------------
# 4. Process data
#----------------------------------------------------

data_list <- list()
counter <- 1

for (i in 1:length(dates_post)){  ## process monthly data
  for (j in 1:length(vars)){  ## process variable 
    
    
    # # import raster
    # input_dir_month <- paste0(input_dir, format(idate, "%Y%m%d"))
    # tif <- list.files(input_dir_month, full.names = TRUE, pattern = sprintf("%s_%s_dens.tif$", format(idate, "%Y%m%d"), jvar))
    # cnt <- raster(tif)
    
    # set data to import
    idate_post <- dates_post[i]
    idate_pre <- dates_pre[i]
    jvar <- vars[j]
    print(paste("processing", idate_post, jvar))
    
    # import raster
    tif_file <- list.files(input_dir, full.names = TRUE, recursive=TRUE, pattern=sprintf("%s_%s_%s_delta.tif", jvar, idate_post, idate_pre))
    cnt <- raster(tif_file)
    
    
    # combine with EEZ
    b <- brick(r_eez, cnt)
    names(b) <- c("mrgid", "delta")
    
    # convert to data.frame
    df <- as.data.frame(b, xy=TRUE, na.rm=TRUE) %>%
      mutate(idcell = cellFromXY(r_eez, cbind(x, y)), # Add cell ID
             date = idate,
             var = jvar)  
    
    data_list[[counter]] <- df
    counter <- counter + 1
  }
}



#----------------------------------------------------
# 5. Prepare data.frame for later analysis
#----------------------------------------------------

## combine data
data <- rbindlist(data_list)

## rename and reorder
data <- data %>% 
  rename(lon = x, lat = y, MRGID = mrgid) %>%
  dplyr::select(idcell, lon, lat, date, var, delta, MRGID)

## combine with EEZ
eez <- eez %>% st_drop_geometry() # remove geometry to convert to data.frame
data <- data %>%
  left_join(dplyr::select(eez, MRGID, TERRITORY1, SOVEREIGN1, ISO_SOV1), by = c("MRGID" = "MRGID"))

# select combinations of EEZ and variable with data for all months (n=128)
select_countries <- data %>%
  #filter(vessels_km2 > 0) %>%
  group_by(idcell, var) %>%
  summarize(n_months = n()) %>%
  dplyr::filter(n_months == 6)


# filter data
data <- data %>%
  right_join(select_countries, by = c("idcell" = "idcell", "var" = "var"))

# export data.frame
write.csv(data, paste0(out_dir, "eez_ais_delta.csv"), row.names = FALSE)
