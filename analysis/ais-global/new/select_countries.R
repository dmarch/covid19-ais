#------------------------------------------------------------------
# select_countries
#------------------------------------------------------------------
# Select countries for analysis with following criteria:
# - EEZ type and size
# - Stringency Index data availability



library(sf)
library(dplyr)
library(lubridate)

#----------------------------------------------
# Criteria 1. EEZ type and size
#----------------------------------------------

# import EEZ
# source: marineregions
eez <- st_read("data/input/marine_regions/World_EEZ_v11_20191118_gpkg/eez_v11_covid.gpkg")
nrow(eez) # 278

eez <- eez %>%
  mutate(TERRITORY1 = as.character(TERRITORY1),
         SOVEREIGN1 = as.character(SOVEREIGN1)) %>%
  filter(!POL_TYPE %in% c("Joint regime", "Overlapping claim"), # remove joint regimes
         TERRITORY1!="Antarctica",  # remove Antarctica
         TERRITORY1 == SOVEREIGN1,  # select main territories (excludes overseas)
         AREA_KM2>(769*3)) # remove EEZ smaller than 3 grid sizes at equator
nrow(eez) # 143 EEZ selected

#----------------------------------------------
# Criteria 2. Stringency index
#----------------------------------------------
output_data <- "data/out/stringency"
data <- read.csv(paste0(output_data, "/stringency.csv"))
data$Date <- ymd(data$Date)
length(unique(data$CountryCode)) # 184 countries with Stringecy Index

data <- data %>%
  filter(CountryCode %in% eez$ISO_SOV1,
         !is.na(StringencyIndex))
length(unique(data$CountryCode)) # 124 countries with EEZ




#----------------------------------------------
# Export list of countries
#----------------------------------------------

selected_iso <- unique(data$CountryCode)
write.csv(selected_iso, "data/out/stringency/countries.csv", row.names=FALSE)


