#-------------------------------------------------------------------------------------------------
# temporal_delta.R        Assess temporal variability of deltas across multiple months
#-------------------------------------------------------------------------------------------------
# Calculate sum across all period
# Calculate average and SD
# Find month with higher decrease

source("scr/fun_common.R")  # bb()
source("scr/fun_ais.R")
library(ggplot2)
library(ggridges)
library(egg)
library(raster)
library(dplyr)

# set input directory
input_dir <- "data/out/ais-global/sog/"

# create output directory
out_dir <- "data/out/ais-global/sog/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


# import ocean mask
r_area <- raster("data/out/ais-global/ocean_area.tif")
r_area <- projectRaster(r_area, res = 27750, crs = "+proj=moll +ellps=WGS84", method="ngb")

# set dates
dates_post <- c(
  seq.Date(as.Date("2020-01-01"), as.Date("2020-06-01"), by = "month")
) %>% format("%Y%m%d")

dates_pre <- c(
  seq.Date(as.Date("2019-01-01"), as.Date("2019-06-01"), by = "month")
) %>% format("%Y%m%d")


histo_list <- list()
cnt <- 1

s <- stack()

# import data
for(i in 1:length(dates_post)){
  
  # select dates
  idate_post <- dates_post[i]
  idate_pre <- dates_pre[i]
  
  # import raster
  tif_file <- list.files(input_dir, full.names = TRUE, recursive=TRUE, pattern=sprintf("SOG_%s_%s_delta_mol.tif", idate_post, idate_pre))
  r <- raster(tif_file)
  
  # create stack
  s <- stack(s, r)
  
  # append decrease values into histo
  vals <- na.omit(values(r))
  #sign_vals <- sign(vals)
  #vals_decrease <- abs(vals[which(sign_vals == -1)])
  hist <- data.frame(var = "SOG", date = idate_post, vals = vals)
  histo_list[[cnt]] <- hist
  cnt <- cnt+1
}

# Calculate average and SD
delta_u <- mean(s, na.rm=TRUE)

# Export data
writeRaster(delta_u, paste0(out_dir, sprintf("SOG_delta_u.tif")), overwrite=TRUE)
#writeRaster(delta_sd, paste0(out_dir, sprintf("%s_delta_sd.tif", jvar)), overwrite=TRUE)


#### Plot average
pngfile <- paste0(out_dir, sprintf("SOG_delta_avg.png"))
png(pngfile, width=3000, height=1750, res=300)
plotDelta(delta_u, main = sprintf("Average Delta Jan-Jun 2020 (SOG)"))
dev.off()


## Store values into list
# vals <- values(delta_u)
# sign_vals <- sign(vals)
# vals_decrease <- abs(vals[which(sign_vals == -1)])
# hist <- data.frame(var = jvar, decrease = vals_decrease)
# histo_list[[j]] <- hist


# combine data
data <- data.table::rbindlist(histo_list)



# reorder vessel categories for plots
data$var <- factor(data$var)

# rename months
data$month <- as.character(data$date)
data$month[which(data$month == "20200101")] <- "Jan"
data$month[which(data$month == "20200201")] <- "Feb"
data$month[which(data$month == "20200301")] <- "Mar"
data$month[which(data$month == "20200401")] <- "Apr"
data$month[which(data$month == "20200501")] <- "May"
data$month[which(data$month == "20200601")] <- "Jun"
data$month <- factor(data$month, levels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun"))



## plot
p <- ggplot(data, aes(x = vals, y = month, fill=var)) +
  geom_density_ridges(aes(fill = var, colour=var), scale = 1.2,#, rel_min_height = 0.01,
                      quantile_lines = T, quantiles = 2, alpha=0.3) +
  #scale_fill_brewer(palette = 4) +
  scale_x_continuous(limits = c(-5, 5)) +
  # scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
  #               labels = scales::trans_format("log10", scales::math_format(10^.x)),
  #               limits = c(1e-4, 1e0)) +
  #facet_wrap(var~., ncol=2) +
  xlab("") + ylab("") +
  theme_article() +
  theme(legend.position = "none", panel.grid.major = element_line(colour = "grey90"), panel.ontop = FALSE)


# export multi-panel plot
out_file <- paste0(out_dir, "monthly_histos_sog.png")
ggsave(out_file, p, width=12, height=10, units = "cm")



# ## plot
# p <- ggplot(data, aes(x = decrease, y = month, fill=var)) +
#   geom_density_ridges(aes(fill = var, colour=var), scale = 1.5,#, rel_min_height = 0.01,
#                       quantile_lines = T, quantiles = 2, alpha=.3) +
#   #scale_fill_brewer(palette = 4) +
#   scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
#                 labels = scales::trans_format("log10", scales::math_format(10^.x)),
#                 limits = c(1e-4, 1e0)) +
#   facet_wrap(var~., ncol=2) +
#   xlab("") + ylab("") +
#   theme_article() +
#   theme(legend.position = "none", panel.grid.major = element_line(colour = "grey90"), panel.ontop = FALSE)
# 
# # add letters
# p <- tag_facet(p, open = "", close = "")
# 
# # export multi-panel plot
# out_file <- paste0(out_dir, "monthly_histos.png")
# ggsave(out_file, p, width=20, height=15, units = "cm")










