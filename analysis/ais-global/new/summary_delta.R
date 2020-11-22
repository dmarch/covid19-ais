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
library(dplyr)
library(raster)
library(pals)

# set input directory
input_dir <- "data/out/ais-global/delta/"

# create output directory
out_dir <- "data/out/ais-global/delta_summary/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


# import ocean mask
r_area <- raster("data/out/ais-global/ocean_area.tif")
r_area <- projectRaster(r_area, res = 27750, crs = "+proj=moll +ellps=WGS84", method="ngb")

# select variables to compare
vars <- c("COUNT", "FISHING", "PASSENGER", "CARGO", "TANKER", "OTHER")

# set dates
dates_post <- c(
  seq.Date(as.Date("2020-01-01"), as.Date("2020-06-01"), by = "month")
) %>% format("%Y%m%d")

dates_pre <- c(
  seq.Date(as.Date("2019-01-01"), as.Date("2019-06-01"), by = "month")
) %>% format("%Y%m%d")


histo_list <- list()
cnt <- 1

for (j in 1:length(vars)){
  
  # select variable
  jvar <- vars[j]
  s <- stack()
  
  # import data
  for(i in 1:length(dates_post)){
    
    # select dates
    idate_post <- dates_post[i]
    idate_pre <- dates_pre[i]
    
    # import raster
    tif_file <- list.files(input_dir, full.names = TRUE, recursive=TRUE, pattern=sprintf("%s_%s_%s_delta_mol.tif", jvar, idate_post, idate_pre))
    r <- raster(tif_file)
    
    # create stack
    s <- stack(s, r)
    
    # append decrease values into histo
    vals <- na.omit(values(r))
    #sign_vals <- sign(vals)
    #vals_decrease <- abs(vals[which(sign_vals == -1)])
    hist <- data.frame(var = jvar, date = idate_post, vals = vals)
    histo_list[[cnt]] <- hist
    cnt <- cnt+1
  }
  
  # Accumulate change through the January-June 2020 period
  # Calculate sum
  delta_sum <- sum(s, na.rm=TRUE) %>% mask(r_area)
  #delta_sum1 <- calc(delta_sum, fun=function(x) sign(x)*log1p(abs(x)))
  #delta_sum2 <- calc(delta_sum, fun=function(x) asinh(exp(5)*x)) #exp(2)*
  
  # Calculate average and SD
  delta_u <- mean(s, na.rm=TRUE)
  #delta_sd <- calc(s, fun=sd, na.rm=TRUE)
  
  
  # Count percentage of negative values
  neg <- s
  neg[neg >= 0] <- NA
  neg[neg < 0] <- 1
  total_neg <- 100 * (sum(neg, na.rm=TRUE) / nlayers(neg))
  total_neg <- mask(total_neg, r_area)
  
  # Export data
  writeRaster(delta_sum, paste0(out_dir, sprintf("%s_delta_sum.tif", jvar)), overwrite=TRUE)
  writeRaster(delta_u, paste0(out_dir, sprintf("%s_delta_u.tif", jvar)), overwrite=TRUE)
  #writeRaster(delta_sd, paste0(out_dir, sprintf("%s_delta_sd.tif", jvar)), overwrite=TRUE)
  
  #### Plot accumulated
  pngfile <- paste0(out_dir, sprintf("%s_delta_sum.png", jvar))
  png(pngfile, width=3000, height=1750, res=300)
  plotDelta(delta_sum, main = sprintf("Accumulated Delta Jan-Jun 2020 (%s)", jvar))
  dev.off()
  
  #### Plot average
  pngfile <- paste0(out_dir, sprintf("%s_delta_avg.png", jvar))
  png(pngfile, width=3000, height=1750, res=300)
  plotDelta(delta_u, main = sprintf("Average Delta Jan-Jun 2020 (%s)", jvar))
  dev.off()
  
  #### Plot percentage of negative values
  minval <- 0
  maxval <- 100
  pngfile <- paste0(out_dir, sprintf("%s_percent_negative.png", jvar))
  png(pngfile, width=3000, height=1750, res=300)
  plotDensMol(r = total_neg, zlim = c(minval, maxval), mollT = FALSE, logT = FALSE,
              col = (viridis(101)), main = sprintf("Persistence of decreases (percentage, %s)", jvar),
              axis_at = c(minval, maxval), axis_labels = c(round(minval,3), round(maxval,3)))
  dev.off()
  
  ## Store values into list
  # vals <- values(delta_u)
  # sign_vals <- sign(vals)
  # vals_decrease <- abs(vals[which(sign_vals == -1)])
  # hist <- data.frame(var = jvar, decrease = vals_decrease)
  # histo_list[[j]] <- hist
}

# combine data
data <- data.table::rbindlist(histo_list)
write.csv(data, )


# reorder vessel categories for plots
data$var <- factor(data$var, levels=c("COUNT", "CARGO", "TANKER", "PASSENGER", "FISHING", "OTHER"))

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
  geom_density_ridges(aes(fill = var, colour=var), scale = 1.5,#, rel_min_height = 0.01,
                      quantile_lines = T, quantiles = 2, alpha=.3) +
  #scale_fill_brewer(palette = 4) +
  scale_x_continuous(#"Modulus transform scale",
                     trans = scales::modulus_trans(-1),
                     limits = c(-0.01, 0.01)) +
  # scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
  #               labels = scales::trans_format("log10", scales::math_format(10^.x)),
  #               limits = c(1e-4, 1e0)) +
  facet_wrap(var~., ncol=2) +
  xlab("") + ylab("") +
  theme_article() +
  theme(legend.position = "none", panel.grid.major = element_line(colour = "grey90"), panel.ontop = FALSE)


# add letters
p <- tag_facet(p, open = "", close = "")

# export multi-panel plot
out_file <- paste0(out_dir, "monthly_histos_posneg.png")
ggsave(out_file, p, width=20, height=15, units = "cm")


#### Boxplots

pd = position_dodge(width = 0.8)
p <- ggplot(filter(data), aes(x = month, y = vals)) +
  geom_hline(aes(yintercept = 0), color = "red") +
  stat_boxplot(geom="errorbar", position=pd, width=0.2, lwd=0.3) +
  geom_boxplot(position=pd, outlier.shape = NA) +
  scale_y_continuous(#"Modulus transform scale",
    #trans = scales::modulus_trans(-1),
    limits = c(-0.012, 0.012)) +
  facet_wrap(var ~ .,  ncol = 2) +
  xlab("") + ylab("") +
    theme_bw(base_size = 12, base_family = "") +
    theme_article()+
    theme(legend.position =  "none")  


# add letters
p <- tag_facet(p, open = "", close = "")

# export multi-panel plot
out_file <- paste0(out_dir, "monthly_boxplots.png")
ggsave(out_file, p, width=20, height=20, units = "cm")

# dat <- ggplot_build(p)$data[[3]]
# dat$col <- ifelse(dat$middle < 0, "red", "grey20")
# p + geom_segment(data=dat, aes(x=xmin, xend=xmax,
#                                y=middle, yend=middle, colour=col), size=1) +
#   scale_color_manual(values=c("grey20", "red"))

  





