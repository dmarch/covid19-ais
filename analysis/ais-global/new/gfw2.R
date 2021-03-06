#------------------------------------------------------------
# GFW
#------------------------------------------------------------


# Create average map of fishing effort 2012-2016
# Correlation with shipping maps



# Load packages
library(tidyverse) # for general data wrangling and plotting
library(furrr) # for parallel operations on lists
library(lubridate) # for working with dates
library(sf) # for vector data 
library(raster) # for working with rasters
library(maps) # additional helpful mapping packages
library(maptools)
library(rgeos)


# create output directory
out_dir <- "data/out/ais-global/gfw/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)




#---------------------------------------------------------------------
# 1. Create average map of fishing effort 2012-2016
#---------------------------------------------------------------------

data_dir <- "data/input/gfw/daily-csvs-100-v1/fishing_effort/daily_csvs/"

# Create dataframe of filenames dates and filter to date range of interest
effort_files <- tibble(
  file = list.files(data_dir, 
                    pattern = '.csv', recursive = T, full.names = T),
  date = ymd(str_extract(file, 
                         pattern = '[[:digit:]]{4}-[[:digit:]]{2}-[[:digit:]]{2}')))

# Generate a vector of dates of interest using ymd from lubridate
effort_dates <- seq(ymd('2016-01-01'), ymd('2016-12-31'), by='days')

# Filter to files within our date range of interest
effort_files <- filter(effort_files, date %in% effort_dates)



# Read in data (uncomment to read in parallel)
plan(multisession) # Windows users should change this to plan(multisession)
effort_df <- furrr::future_map_dfr(effort_files$file, .f = read_csv)

# Add date information
effort_df <- effort_df %>% 
  mutate(year  = year(date),
         month = month(date))

# Specify new (lower) resolution in degrees for aggregating data
res <- 0.25

# Transform data across all fleets and geartypes
effort_df <- effort_df %>% 
  mutate(
    # convert from hundreths of a degree to degrees
    lat_bin = lat_bin / 100, 
    lon_bin = lon_bin / 100,
    # calculate new lat lon bins with desired resolution
    lat_bin = floor(lat_bin/res) * res + 0.5 * res, 
    lon_bin = floor(lon_bin/res) * res + 0.5 * res)

# Re-aggregate the data to 0.25 degrees
effort_df <- effort_df %>% 
  group_by(date, year, month, lon_bin, lat_bin, flag, geartype) %>% 
  summarize(vessel_hours = sum(vessel_hours, na.rm = T),
            fishing_hours = sum(fishing_hours, na.rm = T),
            mmsi_present  = sum(mmsi_present, na.rm = T))


# Aggregate data across all fleets and geartypes
effort_all <- effort_df %>% 
  group_by(lon_bin,lat_bin) %>% 
  summarize(mmsi_present = mean(mmsi_present, na.rm = T),
            vessel_hours = sum(vessel_hours, na.rm = T),
            fishing_hours = sum(fishing_hours, na.rm = T),
            log_fishing_hours = log10(sum(fishing_hours, na.rm = T))) %>% 
  ungroup() %>% 
  mutate(log_fishing_hours = ifelse(log_fishing_hours <= 1, 1, log_fishing_hours),
         log_fishing_hours = ifelse(log_fishing_hours >= 5, 5, log_fishing_hours))# %>% 
  #filter(fishing_hours >= 24)



# Linear green color palette function
effort_pal <- colorRampPalette(c('#0C276C', '#3B9088', '#EEFF00', '#ffffff'), 
                               interpolate = 'linear')

# Map fishing effort
p1 <- effort_all %>%
  ggplot() +
  # geom_sf(data = world_shp, 
  #         fill = '#374a6d', 
  #         color = '#0A1738',
  #         size = 0.1) +
  # geom_sf(data = eezs,
  #         color = '#374a6d',
  #         alpha = 0.2,
  #         fill = NA,
  #         size = 0.1) +
  geom_raster(aes(x = lon_bin, y = lat_bin, fill = log_fishing_hours)) +
  scale_fill_gradientn(
    "Fishing Hours",
    na.value = NA,
    limits = c(1, 5),
    colours = effort_pal(5), # Linear Green
    labels = c("10", "100", "1,000", "10,000", "100,000+"),
    values = scales::rescale(c(0, 1))) +
  labs(fill  = 'Fishing hours (log scale)',
       title = 'Global fishing effort in 2016') +
  guides(fill = guide_colourbar(barwidth = 10)) #+
#gfw_theme



# Select the coordinates and the variable we want to rasterize
effort_all_raster <- effort_all %>% 
  dplyr::select(lon_bin, lat_bin, mmsi_present)

effort_all_raster <- rasterFromXYZ(effort_all_raster, 
                                   crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# transform to molleide
gfw <- projectRaster(effort_all_raster, res = 27750, crs = "+proj=moll", method="ngb")

#---------------------------------------------------------------------
# 2. Compare with ship types
#---------------------------------------------------------------------

# set input folder
input_dir <- "data/out/ais-global/density/"

# select variables to compare
vars <- c("FISHING", "PASSENGER", "CARGO", "TANKER", "OTHER")

# select months to process
dates <- c(
  seq.Date(as.Date("2019-01-01"), as.Date("2019-06-01"), by = "month")
) 

# get year and format dates
year <- year(dates)[1]
dates <- dates %>% format("%Y%m%d")



for(j in 1:length(vars)){
  
  # select variable
  jvar <- vars[j]
  
  # create empty stack
  s <- stack()
  
  # loop by month date
  for(i in 1:length(dates)){
    
    idate <- dates[i]
    
    tif_file <- list.files(input_dir, full.names = TRUE, recursive=TRUE, pattern=sprintf("%s_%s_dens_moll.tif$", idate, jvar))
    r <- raster(tif_file)
    s <- stack(s, r)
  }
  
  # calculate average
  u <- mean(s, na.rm=TRUE)
  
  # log-transform
  ulog <- log10(u)
  
  # compare with gfw
  gfwlog <- log10(gfw)
  gfwlog <- resample(gfwlog, ulog, method="ngb")
  gfwlog[gfwlog==0]<-NA
  
  # aggregate datasets at higher resolution

  
  # pearson
  x <- corLocal(ulog, gfwlog, test=TRUE )

  # scatterplot
  comb <- stack(ulog, gfwlog)
  df <- as.data.frame(comb, na.rm=TRUE)
  plot(df$layer.1, df$layer.2)
 
  library(viridis) 
  library(MASS)
  library(ggplot2)
  get_density <- function(x, y, ...) {
    dens <- MASS::kde2d(x, y, ...)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    return(dens$z[ii])
  }
  df$density <- get_density(df$layer.1, df$layer.2, n = 100)
  ggplot(df) +
    geom_point(aes(x=layer.1, y=layer.2, color = density)) + scale_color_viridis()
  