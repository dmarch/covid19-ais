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



# create output directory
out_dir <- "results/ais-global/zonal/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


# select variables to summarize
vars <- c("COUNT", "FISHING", "PASSENGER", "CARGO", "TANKER", "OTHER")
dates <- c("20200401")



#------------------------------------------------
# Zona statistics to average change by latitude
#------------------------------------------------

# Prepare cluster for parallel computing
cl <- makeCluster(6)
registerDoParallel(cl)

#for (i in 1:length(vars)){
data <- bind_rows(foreach(i = 1:length(vars), .packages=c("dplyr", "raster")) %dopar% { 
  
  ## ser variable
  ivar <- vars[i]
  
  ## import delta
  delta_file <- paste0("data/out/ais-global/delta/", sprintf("%s_delta.tif", ivar))
  delta <- raster(delta_file)
  
  ## import density 2020
  dens_file <- paste0("data/out/ais-global/density/20200401_", sprintf("%s_dens.tif", ivar))
  dens <- raster(dens_file)
  
  ## create latitude raster
  yzone <- init(delta, v="y")
  
  ## zonal statistics
  delta_mean <- zonal(delta, yzone, mean, digits=1, na.rm=TRUE)
  delta_sd <- zonal(delta, yzone, sd, digits=1, na.rm=TRUE)
  dens_mean <- zonal(dens, yzone, mean, digits=1, na.rm=TRUE)
  dens_sd <- zonal(dens, yzone, sd, digits=1, na.rm=TRUE)
  
  ## combine data
  df <- data.frame(var = ivar, delta_mean) %>% cbind(delta_sd[,2]) %>% cbind(dens_mean[,2]) %>% cbind(dens_sd[,2])
  names(df) <- c("var", "latitude", "delta_mean", "delta_sd", "dens_mean", "dens_sd")
  df
})

stopCluster(cl)  # Stop cluster


#------------------------------------------------
# Plots
#------------------------------------------------

# set order
data$var <- factor(data$var, levels =c("COUNT", "CARGO", "TANKER", "PASSENGER", "FISHING", "OTHER"))

p <- ggplot(filter(data, latitude > -70 & latitude < 70), aes(x = latitude, group = var)) + #var!="COUNT", 
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_line(aes(y = delta_mean), colour="#3182bd", size = 0.5) +
  ylab(expression(Traffic~density~change~(Delta~vessels/km^2))) + xlab("Latitude") +
  ylim(-0.03, 0.01) +
  facet_wrap(var ~ ., ncol = 1) +
  theme_article() +
  theme(legend.position = "right")

p <- tag_facet(p, open = "", close = "", tag_pool = letters[-c(1:6)])

# export multi-panel plot
out_file <- paste0(out_dir, "delta_latitude.png")
ggsave(out_file, p, width=10, height=20, units = "cm")


