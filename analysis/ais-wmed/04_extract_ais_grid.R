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
library(janitor)
library(pals)
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
r_mask <- raster("data/out/ais-global/oceanmask.nc")


# r <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90,
#             crs=CRS("+proj=longlat +datum=WGS84"),
#             resolution=c(0.25,0.25), vals=NULL)

# crop raster
r_area <- crop(r_area, bb_sp)
r <- r_area
r[] <- 0

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

# reformat month
dates_df$ym <- format(dates_df$ym, "%Y%m%d")


#-----------------------------------------------------------------
# Get number of vessels per day and ship type
#-----------------------------------------------------------------

## Prepare clusters
cores <- 12 # detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)

data(countriesHigh, package = "rworldxtra", envir = environment())

unique_months <- unique(dates_df$ym)
for(j in 1:length(unique_months)){  # 14 months to analyse
  
  print(paste("Processing month", j, "from", length(unique_months)))
  
  # select dates for month i
  jmonth <- unique_months[j]
  jdates <- dates_df %>% dplyr::filter(ym == jmonth) %>% dplyr::select(date)
  dates <- jdates$date

  # create output folder for selected month
  out_dir_month <- paste0(output_data, "/", jmonth, "/")
  if (!dir.exists(out_dir_month)) dir.create(out_dir_month, recursive = TRUE)
  
  # process AIS data
  cnt <- rbindlist(foreach(i=1:length(dates), .packages = c("lubridate", "dplyr", "janitor"))  %dopar% {
    
    # Import AIS
    year <- year(dates[i])
    data <- import_ais(date = dates[i], path_dyn = "D:/Data/AIS/socib/L0/dynamic/", path_sta = paste0("D:/Data/AIS/socib/L1/static/", year, "/"))
    
    # Exclude codes outside the correct numerical range
    # (i.e. keep MMSI codes with first digits between 2 and 7)
    data <- filter(data, mmsi >= 200000000, mmsi < 800000000)
    
    # exclude anchored or moored vessels
    data <- filter(data, !status %in% c(1,5), speed >  1.852)
    
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
  
  # get grid cell id for each vessel observation
  cnt$cellID <- cellFromXY(r, cbind(cnt$lon, cnt$lat)) 
  cnt$cellID <- as.character(cnt$cellID)
  
  # filter out observations that are outside the grid
  cnt <- filter(cnt, !is.na(cellID))
  
  # calculate number of unique vessels per ship type per grid cell
  cell_type <- cnt %>%
    group_by(cellID, type) %>%
    summarize(n = length(unique(mmsi)))
  
  # calculate number of ALL unique vessels
  cell_all <- cnt %>%
    group_by(cellID) %>%
    summarize(n = length(unique(mmsi))) %>%
    mutate(type = "COUNT") %>%
    relocate(type, .after = cellID)

  # combine global summariesand summaries per ship type
  cell_sum <- bind_rows(cell_all, cell_type)
  cell_sum$cellID <- as.numeric(cell_sum$cellID)

  # list unique categories  
  types <- unique(cell_sum$type)
  
  for (k in 1:length(types)){
    
    ktype <- types[k]
    tdata <- dplyr::filter(cell_sum, type == ktype)
    
    # rasterize
    rcount <- r
    rcount[tdata$cellID] <- tdata$n
    
    # convert to densities
    rcount[rcount==0] <- NA
    rdens <- rcount / r_mask

    # export raster
    writeRaster(rdens, paste0(out_dir_month, sprintf("%s_%s_dens.tif", jmonth, ktype)), overwrite=TRUE)

    # plot density
    minval <- minValue(rdens)
    maxval <- maxValue(rdens)
    pngfile <- paste0(out_dir_month, sprintf("%s_%s_dens.png", jmonth, ktype))
    png(pngfile, width=3000, height=1750, res=300)
    plot(rdens, col=rev(brewer.spectral(101)), main = sprintf("Density %s (%s)", jmonth, ktype))
    plot(countriesHigh, col=NA, border="black", add=TRUE)  # land mask
    # plotDens(r = rdens, zlim = c(minval, maxval), mollT = FALSE, logT = FALSE,
    #             col = rev(brewer.spectral(101)), main = sprintf("Density %s (%s)", jmonth, ktype),
    #             axis_at = c(minval, maxval), axis_labels = c(round(minval,3), round(maxval,3)))
    dev.off()
  }
}


## Stop cluster
stopCluster(cl)

