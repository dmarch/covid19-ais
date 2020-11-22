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
input_dir <- "data/out/ais-global/density/"

# output directory
out_dir <- "data/out/ais-global/eez/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# select variables to summarize
vars <- c("COUNT", "FISHING", "PASSENGER", "CARGO", "TANKER", "OTHER")

# select months to process
dates <- c(
  seq.Date(as.Date("2019-01-01"), as.Date("2019-06-01"), by = "month"),
  seq.Date(as.Date("2020-01-01"), as.Date("2020-06-01"), by = "month")
) %>% format("%Y%m%d")

#-----------------------------------------------
# 2. Calculate area accounting to latitude
#-----------------------------------------------

#create empty raster at 0.25 x 0.25 degrees (resolution of AIS density maps)
r <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90,
            crs=CRS("+proj=longlat +datum=WGS84"),
            resolution=c(0.25,0.25), vals=NULL)

# calculate area for each cell
r_area <- area(r)  # km2

# create higher resolution raster
r_high <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90,
                 crs=CRS("+proj=longlat +datum=WGS84"),
                 resolution=c(0.025,0.025), vals=NULL)


#----------------------------------------------------
# 3. Prepare EEZ
#----------------------------------------------------

# import EEZ
# source: marineregions
eez <- st_read("data/input/marine_regions/World_EEZ_v11_20191118_gpkg/eez_v11_covid.gpkg")
eez <- eez %>%
  filter(!POL_TYPE %in% c("Joint regime")) # remove joint regimes


#----------------------------------------------------
# 4. Process data
#----------------------------------------------------

data_list <- list()
counter <- 1

for (i in 1:length(dates)){  ## process monthly data
  for (j in 1:length(vars)){  ## process variable 
    
    
    # set data to import
    idate <- dates[i]
    jvar <- vars[j]
    print(paste("processing", idate, jvar))
    
    # import raster
    input_dir_month <- paste0(input_dir, idate)
    tif <- list.files(input_dir_month, full.names = TRUE, pattern = sprintf("%s_%s_dens.tif$", idate, jvar))
    cnt <- raster(tif)
    
    var_list <- rbindlist(foreach(k=1:length(eez$MRGID), .packages = c("dplyr", "raster", "fasterize", "sf")) %dopar% { # length(eez$MRGID)
      
      #print(paste("eez", k, "out", length(eez$MRGID)))
      
      # select eez
      imrgid <- eez$MRGID[k]
      #ieez <- dplyr::filter(eez, MRGID == imrgid)
      ieez <- eez[eez$MRGID == imrgid,]
      
      # rasterize into high resolution raster
      r_eez <- fasterize(ieez, r_high)
      
      # crop raster to rasterized extension
      r_eez_crop <- crop(r_eez, extent(ieez)+1)
      
      # return to original resolution and calculate surface
      r_eez_cover <- aggregate(r_eez_crop, fact=10, fun=sum)
      
      # align with area
      r_area_crop <- crop(r_area, r_eez_cover)
      r_eez_align <- resample(r_eez_cover, r_area_crop)
      
      # calculate area in EEZ cells
      r_eez_area <- r_area_crop * ((r_eez_align)/100)
      
      # create mask per eez
      r_eez_mask <- (r_eez_area/r_eez_area)
      
      # get ship counts
      r_eez_cnt <- cnt * r_eez_mask
      
      # set zero values
      r_eez_cnt[is.na(r_eez_cnt)]<-0
      r_eez_cnt <- r_eez_cnt * r_eez_mask
      
      # create stack with surface and count of vessels
      s <- stack(r_eez_area, r_eez_cnt)
      names(s) <- c("eez_km2", "vessels_km2")
      
      # convert to data.frame
      df <- as.data.frame(s, xy=TRUE, na.rm=TRUE)
      df <- mutate(df, date = idate, var = jvar, MRGID = imrgid)
      df
      #data_list[[counter]] <- df
      #counter <- counter + 1
    })
    
    data_list[[counter]] <- var_list
    counter <- counter + 1
  }
}

## Stop cluster
stopCluster(cl)


#----------------------------------------------------
# 5. Prepare data.frame for later analysis
#----------------------------------------------------

## combine data
data <- rbindlist(data_list)

## Add cell ID
data$idcell <- cellFromXY(cnt, cbind(data$x, data$y))

## estimate number of vessels
data$vessels <- round(data$vessels_km2 * data$eez_km2)

## rename and reorder
data <- data %>% 
  rename(lon = x, lat = y) %>%
  dplyr::select(idcell, lon, lat, date, var, vessels, vessels_km2, eez_km2, MRGID)

## combine with EEZ
eez <- eez %>% st_drop_geometry() # remove geometry to convert to data.frame
data <- data %>%
  left_join(dplyr::select(eez, MRGID, TERRITORY1, SOVEREIGN1, ISO_SOV1), by = c("MRGID" = "MRGID"))

# filter cells with NA
data <- filter(data, !is.na(vessels))

# export data.frame
write.csv(data, paste0(out_dir, "eez_ais_density.csv"), row.names = FALSE)
