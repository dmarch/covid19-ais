# calculate annual occupancy for the time series
# compare with our estimate in 2019
# can we predict expected in 2020?
# we cannot compare by sector


library(raster)
library(ncdf4)
library(RColorBrewer)
source("https://raw.githubusercontent.com/dmarch/abigoos/master/R/utils.R")





## get the data and create a raster stack
gridded_nc <- "data/input/tournadre/shiptarfile1992-2020/prod_gridded_1992-2020.nc"
s <- stack(gridded_nc, varname="nships_smoothed") #nships_smoothed

# pre-process data
pre <- raster::rotate(flip(flip(t((s)), direction = 'y'), direction = 'x'))
pre[pre==0] <- NA 

# transform to Mollweide
crs(pre) <- "+proj=longlat +datum=WGS84"  # set projection
pre_mol <- projectRaster(pre, res = 27750, crs = "+proj=moll +ellps=WGS84", method="ngb")  # transform to mollweide

# convert names to year (equal to days since 1992)
names(pre_mol) <- round(as.numeric(sub("X", "", names(s)))/365.2422)


# calculate occupancy for each year
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

plot(data$year, data$occ_km2)


# 90,362,984 km2 in 2019

# 291,037,421 km2 from our data

# 58,218,265
https://otexts.com/fpp2/stochastic-and-deterministic-trends.html

data <- filter(data, year < 2020)


library(forecast)
library(tsbox)


# this function plots a forcast evaluation
forecast_eval = function(forecast, model_name = "", eval_set = test){
  acc = accuracy(forecast, eval_set)
  acc = round(acc[6], 2) 
  model_name = model_name
  forecast %>% 
    autoplot() + 
    autolayer(test, series = "Test set")+
    theme_classic()+
    labs(title = paste(model_name, " model of by ", acc, " on average across test set", sep = ""))
}



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
  forecast(h = 5, level = 0.85)
aarim %>% 
  autoplot() + 
  theme_article()

plot(aarim)


# Calculate % change for predicted 2020 in relation to 2019

x1 <- data$occ_km2[data$year==2019]
x2 <- aarim$mean[1]
percentage_change <- 100*(x2-x1)/x1  # 2.53204