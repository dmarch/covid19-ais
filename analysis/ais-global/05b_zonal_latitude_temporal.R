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


# set input directory
input_dir <- "data/out/ais-global/delta/"


# create output directory
out_dir <- "results/ais-global/zonal/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


# select variables to summarize
vars <- c("COUNT", "FISHING", "PASSENGER", "CARGO", "TANKER", "OTHER")


# set dates
dates_post <- c(
  seq.Date(as.Date("2020-01-01"), as.Date("2020-07-01"), by = "month")
) %>% format("%Y%m%d")

dates_pre <- c(
  seq.Date(as.Date("2019-01-01"), as.Date("2019-07-01"), by = "month")
) %>% format("%Y%m%d")



#------------------------------------------------
# Zona statistics to average change by latitude
#------------------------------------------------

#for (i in 1:length(vars)){
#data <- bind_rows(foreach(j = 1:length(vars) & i = 1:length(dates_post), .packages=c("dplyr", "raster")) %dopar% { 

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
    tif_file <- list.files(input_dir, full.names = TRUE, recursive=TRUE, pattern=sprintf("%s_%s_%s_delta.tif", jvar, idate_post, idate_pre))
    delta <- raster(tif_file)
    
    ## create latitude raster
    yzone <- init(delta, v="y")
    
    ## zonal statistics
    delta_mean <- zonal(delta, yzone, mean, digits=1, na.rm=TRUE)
    delta_sd <- zonal(delta, yzone, sd, digits=1, na.rm=TRUE)
  
    ## combine data
    df <- data.frame(var = jvar, date = idate_post, delta_mean) %>% cbind(delta_sd[,2])
    names(df) <- c("var", "date", "latitude", "delta_mean", "delta_sd")
    data_list[[cnt]] <- df
    cnt <- cnt+1
  }
}

data <- rbindlist(data_list)


#------------------------------------------------
# Summarize per months
#------------------------------------------------

data <- data %>%
  group_by(var, latitude) %>%
  summarize(mean = mean(delta_mean, na.rm=TRUE),
            sd = sd(delta_mean, na.rm=TRUE))


#------------------------------------------------
# Plots
#------------------------------------------------

# set order
data$var <- factor(data$var, levels =c("COUNT", "CARGO", "TANKER", "PASSENGER", "FISHING", "OTHER"))

p <- ggplot(filter(data, latitude > -70 & latitude < 70), aes(x = latitude)) + #var!="COUNT", 
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha=.2, linetype=0, fill = "steelblue") +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_line(aes(y = mean), colour="#3182bd", size = 0.5) + # colour="#3182bd", 
  scale_colour_brewer(palette = "Blues") +
  ylab(expression(Traffic~density~change~(Delta~vessels/km^2))) + xlab("Latitude") +
  ylim(-0.03, 0.01) +
  facet_wrap(var ~ ., ncol = 1) +
  theme_article() +
  theme(legend.position = "right",
        plot.margin = unit(c(10,10,10,10), "points"))

p <- tag_facet(p, open = "", close = "", tag_pool = letters[-c(1:6)])

# export multi-panel plot
out_file <- paste0(out_dir, "delta_latitude.png")
ggsave(out_file, p, width=10, height=20, units = "cm")


