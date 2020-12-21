#-----------------------------------------------------------------------------
# summary_delta_by_eez       Summarize data per EEZ
#-----------------------------------------------------------------------------

library(dplyr)
library(lubridate)
library(reshape2)
library(lmerTest)
library(lme4)
library(sjPlot)
library(sjmisc)
library(ggplot2)
require(lattice)  
library(ggeffects)
library(pdp)
library(xlsx)
library(merTools)
library(HLMdiag)
library(DHARMa)
library(sf)
#library(mixedup) # remotes::install_github('m-clark/mixedup')


## Set output plots
output_data <- "results/ais-global/eez"
if (!dir.exists(output_data)) dir.create(output_data, recursive = TRUE)



#---------------------------------------------
# Import and prepare density data 
#---------------------------------------------
# import data
# delta data at grid cell level
data <- read.csv("data/out/ais-global/eez/eez_ais_delta.csv")
data$date <- ymd(data$date)


# filter by ship type
vars <- c("COUNT",  "CARGO", "TANKER","PASSENGER", "FISHING","OTHER")


#---------------------------------------------
# Summarize data per EEZ
#---------------------------------------------

sum_eez <- data %>%
  group_by(MRGID, TERRITORY1, SOVEREIGN1, ISO_SOV1, var) %>%
  summarize(mean = mean(delta),
            sd = sd(delta),
            median = median(delta),
            iqr = quantile(delta, 0.75)-quantile(delta, 0.25))



#---------------------------------------------
# Plot maps
#---------------------------------------------

# import country map
library(tmap)
data("World")
#World <- World %>% filter(name != "Antarctica")


# create box
box <- bb(xmin = -180, xmax = 180, ymin = -90, ymax = 90, crs="+proj=moll +ellps=WGS84")


# # import landmask
# data(countriesHigh, package = "rworldxtra", envir = environment())
# countriesHigh <- spTransform(countriesHigh, CRS("+proj=moll +ellps=WGS84"))
# land <- st_as_sf(countriesHigh)


# combine with countries with EEZ by ISO code
map_data <- left_join(World, sum_eez, by = c("iso_a3" = "ISO_SOV1"))
map_data <- filter(map_data, var =="COUNT")



# set breaks
# breaks <- c(-0.4, -0.1, -0.05, -0.02, -0.01, -0.005, 0, 0.005, 0.01, 0.02, 0.05, 0.1)
# labels <- c("-0.4 to -0.1", "-0.1 to -0.05", "-0.05 to -0.02", "-0.02 to -0.01", "-0.01 to -0.005", "-0.005 to 0.00",
#             "0.00 to 0.005", "0.005 to 0.01", "0.01 to 0.02", "0.02 to 0.05", "0.05 to 0.1")


# plot data (average)
library(pals)
p1 <- tm_shape(World, projection="+proj=moll +ellps=WGS84") +
  tm_polygons(col="grey60", border.alpha = 0) +
  tm_shape(map_data, projection="+proj=moll +ellps=WGS84") +
  tm_polygons("median",
              style="quantile",
              midpoint = 0,
              #breaks = breaks,
              #labels = labels,
              n = 7,
              title=("Delta (median)"),
              palette = "RdBu", #brewer.ylgnbu(10) # YlGnBu, -RdYlBu
              border.col = "grey60",
              border.alpha = 0.3) +
  tm_shape(box) +
  tm_polygons(alpha=0, border.col = "grey60") +
  tm_layout(frame = FALSE, legend.title.size=1, legend.outside = TRUE)

# export plot
out_file <- paste0(output_data, "/", sprintf("%s_eez_delta_median_jan_jun.png", ivar))
tmap_save(tm = p1, filename = out_file, width=25, height=12, units = "cm")



# plot data (average)
library(pals)
p2 <-   tm_shape(World, projection="+proj=moll +ellps=WGS84") +
  tm_polygons(col="grey60", border.alpha = 0) +
  tm_shape(map_data, projection="+proj=moll +ellps=WGS84") +
  tm_polygons("iqr",
              style="quantile",
              #midpoint = 0,
              #breaks = breaks,
              #labels = labels,
              title=("Delta (IQR)"),
              palette = "YlGnBu", #brewer.ylgnbu(10) # YlGnBu, -RdYlBu
              border.col = "grey60",
              border.alpha = 0.3) +
  tm_shape(box) +
  tm_polygons(alpha=0, border.col = "grey60") +
  tm_layout(frame = FALSE, legend.title.size=1, legend.outside = TRUE)

# export plot
out_file <- paste0(output_data, "/", sprintf("%s_eez_delta_iqr_jan_jun.png", ivar))
tmap_save(tm = p2, filename = out_file, width=25, height=12, units = "cm")



