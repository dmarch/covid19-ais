




# Load libraries
library(tidyverse)
library(rgee)
library(sf)
library(data.table)
library(reshape2)

# Initialize Earth Engine
ee_Initialize()


## import maritime chokepoints
choke_buff <- st_read("data/input/chokepoints/chokepoints_buff_geo.gpkg")

## Cloud mask function
## Select pixels with less than 10% cloud cover
cloudMask <- function(image) {
  image$updateMask(image$select("cloud_fraction")$lt(0.1))
}

## Ocean mask
datamask <- ee$Image('UMD/hansen/global_forest_change_2015')$ # Hansen et al. forest change dataset.
  select('datamask')  # Select the land/water mask.

## Create a binary mask.
## combine both ocean (0) and permanent water bodies (2).
mask <- datamask$eq(0)$Or(datamask$eq(2))

## select months to process
dates <- c(
  seq.Date(as.Date("2019-01-01"), as.Date("2019-06-01"), by = "month"),
  seq.Date(as.Date("2020-01-01"), as.Date("2020-06-01"), by = "month")
) 

data_list <- list()
for (i in 1:length(dates)){
  
  print(paste("Processing date", dates[i]))
  
  # define date range
  start_date <- floor_date(dates[i], 'month') %>% as.character()
  end_date <- (ceiling_date(dates[i], 'month')-1) %>% as.character()
  
  ## Load Tropomi NO2 for month in 2019.
  l3_no2 <- ee$ImageCollection("COPERNICUS/S5P/OFFL/L3_NO2")$
    select(c("tropospheric_NO2_column_number_density","cloud_fraction"))$
    filterDate(start_date, end_date)$
    map(cloudMask)
  
  # calculate monthly mean
  no2_mean <- l3_no2$
    select("tropospheric_NO2_column_number_density")$
    mean()$
    updateMask(mask)
  
  # Extract values
  # calculate average around buffer chokepoint
  choke_no2 <- ee_extract(x = no2_mean, y = choke_buff, sf = FALSE)
  names(choke_no2) <- "no2"
  choke_no2$id <- choke_buff$id
  choke_no2$name <- choke_buff$name
  choke_no2$date <- dates[i]
  data_list[[i]] <- choke_no2
}


# combine data
data <- rbindlist(data_list)
data$year <- year(data$date)
data$month <- month(data$date)

# transform to wide format
wide <- dcast(data, id + name + month ~ year, value.var = "no2") %>%
  rename(dens2019 = `2019`, dens2020 = `2020`) %>%
  mutate(delta = dens2020 - dens2019,
         per = 100*(delta/dens2019),
         change_positive = per > 0)


# plot
p2 <- ggplot(wide, mapping=aes(x = month, y = delta, fill = change_positive)) +
  geom_col(alpha=1, width=0.8) +
  ylab("Relative change (%) in traffic density") + xlab("Month") +
  #scale_x_continuous( breaks = c(1,3,5), labels = month.abb[c(1, 3,5)]) +  # , labels = month.abb[1:6]
  scale_x_continuous(breaks = 1:6) +
  scale_fill_manual(values=c("#e34a33", "#9ecae1")) +
  facet_wrap(name ~ ., ncol = 4) +
  theme_article() +
  guides(fill = FALSE)






## Load Tropomi NO2 for month in 2019.
# no2 <- ee$ImageCollection("COPERNICUS/S5P/OFFL/L3_NO2")$
#   select(c("tropospheric_NO2_column_number_density","cloud_fraction"))$
#   filterDate("2019-04-01", "2019-04-30")$
#   map(cloudMask)$
#   mean()
# 
# ## Load Tropomi NO2 for month in 2020.
# collection2020 <- ee$ImageCollection("COPERNICUS/S5P/OFFL/L3_NO2")$
#   select(c("tropospheric_NO2_column_number_density","cloud_fraction"))$
#   filterDate("2020-04-01", "2020-04-30")$
#   map(cloudMask)









# Extract values
# calculate average
ee_chk_no2 <- ee_extract(x = maskedDiff, y = choke_buff, sf = FALSE)


# # plot map
# Map$setCenter(9.08203, 47.39835, 3)
# Map$addLayer(
#   eeObject = no2_mean,
#   visParams = list(
#     bands = c("tropospheric_NO2_column_number_density"),
#     min = -0.00002,
#     max = 0.00002,
#     palette = c('blue', 'white', 'red')
#   ),
#   name = "NO2 dif"
# )

