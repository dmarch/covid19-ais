#---------------------------------------------------------------------------
# 02_extract_ais                                                        
#---------------------------------------------------------------------------
# This script extracts AIS data from SOCIB database. It cannot be
# reproduced without access to the raw data.
# The output is stored for further reproducibility of next steps.
#
# Note:
# Adding new dates, requires a previous pre-processing of raw AIS data
# from 'socib-ais' repo, following this order:
# 1. "01_check_daily_data.R": inspects raw data 
# 2. "01b_check_daily_data_L0.R": inspects previous pre-processed data
# 3. "02_process_raw_mt.R": process either dynamic and static data.
# 4. "03_process_static_L1.R": process static data to combine daily files into calendar year data
# 5. "01b_check_daily_data_L0.R": run again to incorporate new data
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

output_data <- "data/out/ais-wmed"
if (!dir.exists(output_data)) dir.create(output_data, recursive = TRUE)


#-----------------------------------------------------------------
# Import auxilary data
#-----------------------------------------------------------------

# subset of EEZ for EU zone (see 01_coastal_waters.R)
eez_24nm <- st_read("data/out/ais-wmed/eez_24nm_eu.gpkg") %>%
  dplyr::select(SOVEREIGN1, ISO_SOV1) %>% as_Spatial()

eez_24nm_sf <- st_read("data/out/ais-wmed/eez_24nm_eu.gpkg") %>%
  dplyr::select(SOVEREIGN1, ISO_SOV1) 

#-----------------------------------------------------------------
# Prepare data to import
#-----------------------------------------------------------------

# Set start and end date of data inspection
sDate <- as.Date("2016-01-01")
eDate <- as.Date("2020-12-12") # "2020-06-04"
dates <- seq.Date(sDate, eDate, "1 day")

# Check dates with data
calendar_file <- "../socib-ais/out/ais_daily_L0.csv"
dates <- dates[check_dates(dates, file = calendar_file)] # filter dates by AIS availability



#-----------------------------------------------------------------
# Get number of vessels per day and ship type
#-----------------------------------------------------------------

## Prepare clusters
cores <- 12 # detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)


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
  exclude_type_summary <- c("", "NULL", "Wing in Grnd")
  exclude_type_name <- c("SAR Aircraft")
  data <- filter(data, !ais_type_summary %in% exclude_type_summary, !type_name %in% exclude_type_name)
  
  # reclass categories
  data$type <- NA
  data$type[data$ais_type_summary %in% c("High Speed Craft", "Passenger")] <- "Passenger"
  data$type[data$ais_type_summary %in% c("Other", "Unspecified", "Search and Rescue", "Special Craft", "Tug")] <- "Other"
  data$type[data$ais_type_summary %in% c("Sailing Vessel", "Pleasure Craft")] <- "Recreational"
  data$type[data$ais_type_summary %in% c("Cargo")] <- "Cargo"
  data$type[data$ais_type_summary %in% c("Tanker")] <- "Tanker"
  data$type[data$ais_type_summary %in% c("Fishing")] <- "Fishing"
  
  # filter vessels without static information
  data <- filter(data, !is.na(data$type))
  
  # filter coastal vessels and add EEZ
  eez <- point_on_land(lon = data$lon, lat =data$lat, land = eez_24nm)
  data <- cbind(data, eez)
  
  # filter vessels outside the coastal waters
  data <- filter(data, !is.na(ISO_SOV1))
  
  # count number of vessel by type
  cnt <- data %>%
    group_by(type) %>%  # add ISO_SOV1 to split between countries
    summarize(vessels = length(unique(mmsi)),
              messages = sum(n)) %>% 
    adorn_totals("row", name = "All") 
  
  # Add date
  cnt$date <- dates[i]
  cnt
})


## Stop cluster
stopCluster(cl)

## write csv
outfile <- paste0(output_data, "/vessels_day.csv")
write.csv(cnt, outfile, row.names = FALSE)