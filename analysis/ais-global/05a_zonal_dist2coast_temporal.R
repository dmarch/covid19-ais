#---------------------------------------------------------------------------
# 05a_zonal_dist2coast          Zonal summary by distance to coast gradient
#---------------------------------------------------------------------------
# This script average absolute difference maps by distance to coast


## libraries
library(doParallel)
library(foreach)
library(parallel)
library(dplyr)
library(ggplot2)
library(egg)
library(raster)



# set input directory
input_dir <- "data/out/ais-global/delta/"


# create output directory
out_dir <- "results/ais-global/zonal/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


# select variables to summarize
vars <- c("COUNT", "FISHING", "PASSENGER", "CARGO", "TANKER", "OTHER")


# set dates
dates_post <- c(
  seq.Date(as.Date("2020-01-01"), as.Date("2020-06-01"), by = "month")
) %>% format("%Y%m%d")

dates_pre <- c(
  seq.Date(as.Date("2019-01-01"), as.Date("2019-06-01"), by = "month")
) %>% format("%Y%m%d")


#------------------------------------------------
# Create distance to coast raster
#------------------------------------------------
# Note this step is time consuming

# import ocean mask
r_area <- raster("data/out/ais-global/oceanmask.nc")

# distance to shore
r_area[is.na(r_area)] <- 999999
r_area[r_area < 1000] <- NA
dist2coast <- distance(r_area)

# save raster
writeRaster(dist2coast, "data/out/ais-global/dist2coast.tif")

#------------------------------------------------
# Reproject and reclass distance to coast
#------------------------------------------------

# import raster
dist2coast <- raster("data/out/ais-global/dist2coast.tif")

# transform to mollweide
dcoast_moll <- projectRaster(dist2coast, res = 27750, crs = "+proj=moll +ellps=WGS84", method="bilinear")

## get range values
dcoast_moll[dcoast_moll<0] <- 0
dcoast_moll <- setMinMax(dcoast_moll)
dcoast_max <- signif(maxValue(dcoast_moll), digits=2)

## Reclassify dist2coast at larger intervals
int <- 50000  # set distance interval, in meters
from <- seq(0, dcoast_max, by=int)
to <- seq(int, dcoast_max+int, by=int)
becomes <- seq(int, dcoast_max+int, by=int)
rclmat <- matrix(cbind(from, to, becomes), ncol=3, byrow=FALSE)  # switched 'to' and 'from' due to negative values
rclmat <- rbind(c(-Inf,0,NA), rclmat)  # include NA for surface values
rc <- reclassify(dcoast_moll, rclmat)  # reclassify raster



#------------------------------------------------
# Zona statistics to average change by distance to coast
#------------------------------------------------

# Prepare cluster for parallel computing
#cl <- makeCluster(6)
#registerDoParallel(cl)

#data <- bind_rows(foreach(i = 1:length(vars), .packages=c("dplyr", "raster")) %dopar% { 


data_list <- list()
cnt <- 1  

for (j in 1:length(vars)){
  
  # select variable
  jvar <- vars[j]
  
  # import data
  for(i in 1:length(dates_post)){
    
    # select dates
    idate_post <- dates_post[i]
    idate_pre <- dates_pre[i]
    
    # import raster
    tif_file <- list.files(input_dir, full.names = TRUE, recursive=TRUE, pattern=sprintf("%s_%s_%s_delta_mol.tif", jvar, idate_post, idate_pre))
    delta <- raster(tif_file)

    ## zonal statistics
    delta_mean <- zonal(delta, rc, mean, digits=1, na.rm=TRUE)
    delta_sd <- zonal(delta, rc, sd, digits=1, na.rm=TRUE)

  
  ## combine data
  df <- data.frame(var = jvar, date = idate_post, delta_mean) %>% cbind(delta_sd[,2])
  names(df) <- c("var", "date", "dcoast", "delta_mean", "delta_sd")
  data_list[[cnt]] <- df
  cnt <- cnt+1
  }
}


data <- data.table::rbindlist(data_list)


#------------------------------------------------
# Summarize per months
#------------------------------------------------

data <- data %>%
  group_by(var, dcoast) %>%
  summarize(mean = mean(delta_mean, na.rm=TRUE),
            sd = sd(delta_mean, na.rm=TRUE))


#------------------------------------------------
# Plots
#------------------------------------------------

# transform m to km
data$dcoast <- data$dcoast/1000

# set order
data$var <- factor(data$var, levels =c("COUNT", "CARGO", "TANKER", "PASSENGER", "FISHING", "OTHER"))

p <- ggplot(filter(data, dcoast <= 4800), aes(x = dcoast, group = var)) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha=.2, linetype=0, fill = "steelblue") +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_line(aes(y = mean), colour="#3182bd", size = 0.5) +
  scale_colour_brewer(palette="Set2") +
  ylab(expression(Traffic~density~change~(Delta~vessel~transits/km^2))) + xlab("Distance to coast (km)") +
  facet_wrap(var ~ ., ncol = 1) +
  ylim(-0.03, 0.01) +
  theme_article() +
  theme(legend.position = "none",
        plot.margin = unit(c(10,10,10,10), "points"))

p <- tag_facet(p, open = "", close = "")

# export multi-panel plot
out_file <- paste0(out_dir, "delta_distcoast.png")
ggsave(out_file, p, width=10, height=20, units = "cm")

