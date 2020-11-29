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
library(ggplot2)
library(egg)
source("../socib-ais/scr/processing_tools.R")  # sources from external repo
source("scr/fun_ais.R")


#-----------------------------------------------------------------
# Set paths
#-----------------------------------------------------------------

output_data <- "data/out/ais-wmed/flag"
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
# Import MIDs codes
#-----------------------------------------------------------------
mids <- read.csv("https://raw.githubusercontent.com/michaeljfazio/MIDs/master/mids.csv", header=FALSE)
names(mids) <- c("MID", "ISO2", "ISO3", "unknow", "Country")
mids$MID <- as.character(mids$MID)

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


unique_months <- unique(dates_df$ym)
for(j in 1:length(unique_months)){  # 14 months to analyse
  
  print(paste("Processing month", j, "from", length(unique_months)))
  
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
    
    # select required data for further steps
    data <- dplyr::select(data, mmsi, timestamp, lon, lat, type, flag)
    data
  })
  
  # get country iso from standardize table
  cnt$MID <- substr(cnt$mmsi, start = 1, stop = 3)
  cnt <- cnt %>% left_join(mids, by="MID")
  
  # calculate number of unique vessels per flag type and ship type
  flag_type <- cnt %>%
    filter(!is.na(Country)) %>%
    group_by(type, Country, ISO3) %>%
    summarize(n = length(unique(mmsi))) %>%
    mutate(month = jmonth)
  
  # export table
  outfile <- paste0(output_data, "/", jmonth, "_flag.csv")
  write.csv(flag_type, outfile, row.names = FALSE)

}


## Stop cluster
stopCluster(cl)



#-----------------------------------------------------------------
# Analyse data
#-----------------------------------------------------------------

# import data
loc_files <- list.files(output_data, full.names = TRUE, pattern="flag.csv")
data <- rbindlist(lapply(loc_files, fread))

# calculate totals
data_totals <- data %>%
  group_by(Country, ISO3, month) %>%
  summarize(n = sum(n)) %>%
  mutate(type = "All")
data <- bind_rows(data, data_totals)


# select top 10 countries
top_countries <- data_totals %>%
  group_by(Country) %>%
  summarize(n = sum(n)) %>%
  arrange(desc(n)) %>%
  slice(1:10)


# format date
data$date <- ymd(data$month)



# Plot time series all countries in the same plot
p1 <- ggplot(filter(data, Country %in% top_countries$Country, date >= as.Date("2020-01-01"), type=="All"), aes(x = date, group = Country)) +
  geom_line(aes(y = n, color = Country), size = 1) +
  #geom_vline(xintercept = as.Date("2020-03-11"), linetype="dotted") +
  scale_colour_brewer(palette="Set1") +
  #scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  #ylim(c(0,100))+
  #ylab("Stringency Index") + xlab("") +
  facet_wrap(Country ~ ., ncol = 3) +
  theme_article() +
  theme(legend.position = "none")
  #theme(legend.position = c(0.85, 0.2), legend.title = element_blank())

# Faceted by country
p2 <- ggplot(filter(data, CountryName %in% WMedCountries), aes(x = Date, group = CountryName)) +
  geom_line(aes(y = StringencyIndex, color = CountryName), size = 1) +
  geom_vline(xintercept = as.Date("2020-03-11"), linetype="dotted") +
  scale_colour_brewer(palette="Set2") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  ylab("Stringency Index") + xlab("") +
  facet_wrap(CountryName ~ ., ncol = 3) +
  theme_article() +
  theme(legend.position = "none")

