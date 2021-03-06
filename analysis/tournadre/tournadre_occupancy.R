# calculate annual occupancy for the time series
# compare with our estimate in 2019
# can we predict expected in 2020?
# we cannot compare by sector


library(raster)
library(ncdf4)
library(RColorBrewer)
library(spatialEco)
source("https://raw.githubusercontent.com/dmarch/abigoos/master/R/utils.R")
source("scr/fun_common.R")  # bb()
source("scr/fun_ais.R")




## get the data and create a raster stack
gridded_nc <- "data/input/tournadre/shiptarfile1992-2020/prod_gridded_1992-2020.nc"
s <- stack(gridded_nc, varname="nships_smoothed") #nships_smoothed

# pre-process data
pre <- raster::rotate(flip(flip(t((s)), direction = 'y'), direction = 'x'))
pre[pre==0] <- NA 
crs(pre) <- "+proj=longlat +datum=WGS84"  # set projection

# calculate density data
r_area <- raster::area(pre)
dens <- pre/r_area

# transform to Mollweide
pre_mol <- projectRaster(dens, res = 27750, crs = "+proj=moll +ellps=WGS84", method="ngb")  # transform to mollweide

# convert names to year (equal to days since 1992)
names(pre_mol) <- round(as.numeric(sub("X", "", names(s)))/365.2422)


#---------------------------------------------------------------------
# Mann-Kendall trend (2010-2019)
#---------------------------------------------------------------------
require(raster)
require(Kendall)
library(devtools)

# select years
ships <- subset(pre_mol, 19:28)
#ships <- crop(ships, extent(-1e6,4e6,3e6,6e6))

# find locations without missing data
s_count <- ships/ships
s_count <- sum(s_count, na.rm=TRUE)
s_count[s_count<5] <- NA
s_count <- s_count/s_count
ships <- mask(ships, s_count)

# imputes missing NA values
#ships_s <- smooth.time.series(ships, f = 0.8, smooth.data = FALSE)

# kendall
k <- raster.kendall(ships, p.value=TRUE, z.value=TRUE, 
                    intercept=TRUE, confidence=TRUE, 
                    tau=TRUE)

# keep brick
writeRaster(k, "data/out/ais-global/occupancy/tournadre_trend.grd")
k <- brick("data/out/ais-global/occupancy/tournadre_trend.grd")

#---------------------------------------------------------------------
# Calculate difference 2019 from average 2016-2018
#---------------------------------------------------------------------

# select years
ships1418 <- subset(pre_mol, 23:27)
ships2019 <- pre_mol$X2019

# calculate average previous years
ships_avg <- mean(ships1418, na.rm=TRUE)
ships_sd <- calc(ships1418, fun=sd, na.rm=TRUE)
ships_sd[ships_sd==0]<-NA

# calculate anomaly
delta <- ships2019 - ships_avg
delta_st <- delta/ships_sd
per <- 100*(delta/ships_avg)

#---------------------------------------------------------------------
# Calculate difference 2019 from average 2016-2018
#---------------------------------------------------------------------

#xmin = -1.60, xmax = 9.93, ymax = 35.56, ymin = 43.45



#---------------------------------------------------------------------
# calculate occupancy for each year
#---------------------------------------------------------------------
data_list <- list()
for(i in 1:nlayers(pre_mol)){
  
  # subset raster
  r <- subset(pre_mol, i)
  
  # get year
  iyear <- format(as.Date(names(pre_mol)[i], format = "X%Y"), format = "%Y")

  # calculate occupancy (any ship density)
  r1q <- quantile(r, prob = 0)
  r1ext <- reclassify(r, c(-Inf, r1q, NA, r1q, Inf, 1))  # create a ocean mask
  occ_km2 <- binsurf(r1ext)  # km2
  
  # create data.frame
  df <- data.frame(year = as.numeric(iyear), occ_km2)
  data_list[[i]] <- df
}

# combine data
data <- data.table::rbindlist(data_list)


# export raster for 2019
r2019 <- pre_mol[[28]]
writeRaster(r2019, "data/input/tournadre/density_2019.tif", overwrite=TRUE)


# plot coefficient of variation
minval <- minValue(r2019)
maxval <- maxValue(r2019)
pngfile <- "data/input/tournadre/density_2019.png"
png(pngfile, width=3000, height=1750, res=300)
plotDensMol(r = r2019, zlim = c(minval, maxval), mollT = FALSE, logT = FALSE,
            col = rev(brewer.spectral(101)), main = "Tournadre 2019",
            axis_at = c(minval, maxval), axis_labels = c(round(minval,3), round(maxval,3)))
dev.off()




plot(data$year, data$occ_km2)


# 90,362,984 km2 in 2019

# 291,037,421 km2 from our data

# 58,218,265
#https://otexts.com/fpp2/stochastic-and-deterministic-trends.html

data <- filter(data, year >= 2000, year < 2020) #year >= 2000, 


library(forecast)
library(tsbox)


# this function plots a forcast evaluation
# forecast_eval = function(forecast, model_name = "", eval_set = test){
#   acc = accuracy(forecast, eval_set)
#   acc = round(acc[6], 2) 
#   model_name = model_name
#   forecast %>% 
#     autoplot() + 
#     autolayer(test, series = "Test set")+
#     theme_classic()+
#     labs(title = paste(model_name, " model of by ", acc, " on average across test set", sep = ""))
# }



# time series
ggplot(data = data, aes(year, occ_km2))+
  geom_line()+ 
  labs(x = "Year", y = "Occupancy", title = "Marine traffic occupancy")+
  theme_minimal()

# convert to time-series object
data_ts <- ts_ts(ts_long(data))

# ACF
ggAcf(data_ts) + theme_light() + labs(title = "ACF plot of Marine traffic occupancy")



aarim = data_ts %>% 
  auto.arima() %>% 
  forecast(h = 5, level = c(80, 95))


p <- aarim %>% 
  autoplot() + 
  xlab("")+
  ylab(expression(Occupancy~(km^2))) +
  theme_article() +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))  # top, right, bottom, left

out_file <- "data/out/ais-global/occupancy/arima_plot.png"
ggsave(out_file, p, width=12, height=10, units = "cm")


# Calculate % change for predicted 2020 in relation to 2019

x1 <- data$occ_km2[data$year==2019]
xmean <- aarim$mean[1]
percentage_change <- 100*(xmean-x1)/x1  #2.981492# 2.53204
