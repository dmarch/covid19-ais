#---------------------------------------------------------------------------
# 04_extract_ais_grid                                                        
#---------------------------------------------------------------------------
#
# This script counts the number of unique vessels per ship type per grid cell and month
#
# This script extracts AIS data from SOCIB database. It cannot be
# reproduced without access to the raw data.
# The output is stored for further reproducibility of next steps.
#
# Note:
# Based on 02_extract_ais.R. See required previous steps in that script
#---------------------------------------------------------------------------


library(foreach)
library(parallel)
library(doParallel)
library(data.table)
library(dplyr)
library(lubridate)
library(stringr)
library(raster)
library(sf)
source("../socib-ais/scr/processing_tools.R")  # sources from external repo
source("scr/fun_ais.R")


#-----------------------------------------------------------------
# Set paths
#-----------------------------------------------------------------

output_data <- "data/out/ais-wmed/density"
if (!dir.exists(output_data)) dir.create(output_data, recursive = TRUE)


#-----------------------------------------------------------------
# Create raster template matching S-AIS data
#-----------------------------------------------------------------

# bb: bounding box of the AIS data
bb <- st_as_sfc(st_bbox(c(xmin = -1.60, xmax = 9.93, ymax = 35.56, ymin = 43.45), crs = st_crs(4326)))
bb_sp <- as(bb, "Spatial")

# get raster at same resolution as S-AIS
r_area <- raster("data/out/ais-global/areacell.nc")  # see ais-global/01_ocean_area.R

r <- raster("data/out/ais-global/oceanmask.nc")

r <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90,
            crs=CRS("+proj=longlat +datum=WGS84"),
            resolution=c(0.25,0.25), vals=NULL)

# crop raster
r <- crop(r, bb_sp)


#-----------------------------------------------------------------
# Prepare data to import
#-----------------------------------------------------------------

# Set start and end date of data inspection
sDate <- as.Date("2019-01-01")
eDate <- as.Date("2020-07-31") # "2020-06-04"
dates <- seq.Date(sDate, eDate, "1 day")

# Check dates with data
calendar_file <- "../socib-ais/out/ais_daily_L0.csv"
dates <- dates[check_dates(dates, file = calendar_file)] # filter dates by AIS availability

# create table with dates and months
dates_df <- data.frame(date = dates, month = month(dates), ym = floor_date(dates, "month"))
dates_df <- filter(dates_df, month <= 7)  # filter months of interest


#-----------------------------------------------------------------
# Get number of vessels per day and ship type
#-----------------------------------------------------------------

## Prepare clusters
cores <- 12 # detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)


unique_months <- unique(dates_df$ym)
for(j in 1:length(unique_months)){  # 14 months to analyse
  
  # select dates for month i
  jmonth <- unique_months[j]
  jdates <- dates_df %>% dplyr::filter(ym == jmonth) %>% dplyr::select(date)
  dates <- jdates$date

  # process AIS data
  cnt <- rbindlist(foreach(i=1:length(dates), .packages = c("lubridate", "dplyr", "janitor"))  %dopar% {
    
    # Import AIS
    year <- year(dates[i])
    data <- import_ais(date = dates[i], path_dyn = "D:/Data/AIS/socib/L0/dynamic/", path_sta = paste0("D:/Data/AIS/socib/L1/static/", year, "/"))
    
    # Exclude codes outside the correct numerical range
    # (i.e. keep MMSI codes with first digits between 2 and 7)
    data <- filter(data, mmsi >= 200000000, mmsi < 800000000)
    
    # exclude anchored or moored vessels
    # data <- filter(data, !status %in% c(1,5), speed >  1.852)
    
    # exclude some types
    # exclude_type_summary <- c("", "NULL", "Wing in Grnd")
    # exclude_type_name <- c("SAR Aircraft")
    # data <- filter(data, !ais_type_summary %in% exclude_type_summary, !type_name %in% exclude_type_name)
    
    # reclass categories
    data$type <- NA
    data$type[data$ais_type_summary %in% c("High Speed Craft", "Passenger")] <- "PASSENGER"
    data$type[data$ais_type_summary %in% c("Other", "Unspecified", "Search and Rescue", "Special Craft", "Tug", "Sailing Vessel", "Pleasure Craft")] <- "OTHER"
    #data$type[data$ais_type_summary %in% c("Sailing Vessel", "Pleasure Craft")] <- "Recreational"
    data$type[data$ais_type_summary %in% c("Cargo")] <- "CARGO"
    data$type[data$ais_type_summary %in% c("Tanker")] <- "TANKER"
    data$type[data$ais_type_summary %in% c("Fishing")] <- "FISHING"
    
    # filter vessels without static information
    data <- filter(data, !is.na(data$type))
    
    # filter coastal vessels and add EEZ
    #eez <- point_on_land(lon = data$lon, lat =data$lat, land = eez_24nm)
    #data <- cbind(data, eez)
    
    # filter vessels outside the coastal waters
    #data <- filter(data, !is.na(ISO_SOV1))
    
    # count number of vessel by type
    # cnt <- data %>%
    #   group_by(type) %>%  # add ISO_SOV1 to split between countries
    #   summarize(vessels = length(unique(mmsi)),
    #             messages = sum(n)) %>% 
    #   adorn_totals("row", name = "All") 
    
    # Add date
    #cnt$date <- dates[i]
    #cnt
    
    # select required data for further steps
    data <- dplyr::select(data, mmsi, timestamp, lon, lat, type)
    data
  })
  
  # rasterize the number of vessels
  # get cell id for each anchoring event
  # then get unique combination of mmsi and cell id
  cnt$cellID <- cellFromXY(r, cbind(cnt$lon, cnt$lat)) 
  cnt$cellID <- as.character(cnt$cellID)
  cell_sum <- cnt %>%
    group_by(cellID) %>%
    summarize(n = length(unique(mmsi)))
  cell_sum$cellID <- as.numeric(cell_sum$cellID)
  
  # rasterize
  rcount <- r
  rcount[cell_sum$cellID] <- cell_sum$n
  writeRaster(ranchor_vessels , filename=paste0(path_out_balearics,"/grid/anchor_bal_vessels.nc"), format="CDF", overwrite=TRUE) 
  
  # rasterize by vessel type
  
  cell_sum <- cnt %>%
    group_by(cellID, type) %>%
    summarize(n = length(unique(mmsi)))
  cell_sum$cellID <- as.numeric(cell_sum$cellID)
  
  types <- unique(cell_sum$type)
  
  for (k in 1:length(types)){
    
    tdata <- dplyr::filter(cell_sum, type == types[k])
    
  }
  
  
  
}




## Stop cluster
stopCluster(cl)

## write csv
outfile <- paste0(output_data, "/vessels_day.csv")
write.csv(cnt, outfile, row.names = FALSE)