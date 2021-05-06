# summarize density maps
# calculate average and coefficient of variation of deltas during 2020

library(stringr)
library(raster)
library(lubridate)
library(rgdal)
library(RColorBrewer)
library(pals)
source("scr/fun_common.R")  # bb()
source("scr/fun_ais.R")

# set input folder
input_dir <- "data/out/ais-global/density/"

# create output directory
out_dir <- "data/out/ais-global/density/BP"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# select variables to compare
vars <- c("COUNT", "FISHING", "PASSENGER", "CARGO", "TANKER", "OTHER")
vars <- c("COUNT")

# select months to process
dates <- c(
  seq.Date(as.Date("2019-01-01"), as.Date("2019-06-01"), by = "month")
) 

# get year and format dates
year <- year(dates)[1]
#dates <- dates %>% format("%Y%m%d")


library("rnaturalearthdata")
library("rnaturalearth")
world <- ne_countries(scale = "medium", returnclass = "sp")

# import landmask
land <- crop(world, extent(-30, 15, 10, 45))

# box
box <- bb(xmin = -30, xmax = 15, ymin = 10, ymax = 45, crs="+init=epsg:4326")



for(j in 1:length(vars)){
  
  # select variable
  jvar <- vars[j]
  
  # create empty stack
  s <- stack()
  
  # loop by month date
  for(i in 1:length(dates)){
    
    # set date
    idate <- dates[i]
    
    # import file
    tif_file <- list.files(input_dir, full.names = TRUE, recursive=TRUE, pattern=sprintf("%s_%s_dens.tif$", format(idate, "%Y%m%d"), jvar))
    r <- raster(tif_file)
    s <- stack(s, r) # append to stack
  }
  
  # crop basedmap to cover the southern ocean (between 90S and 40S latitude)
  s <- crop(s, extent(-40, 5, 0, 45))
  
  # reproject raster to polar stereographic projection (EPSG:3031)
  #s <- projectRaster(s, crs = crs("+init=epsg:3031"), method="ngb")
  
  # get common min and max values
  min <- min(minValue(s))
  max <- max(maxValue(s))
  
  
  for (i in 1:nlayers(s)){
    
    # subset raster
    r <- subset(s,i)
    # set date
    idate <- dates[i]
    # plot title
    main = sprintf("Traffic density - %s - (%s)", jvar, idate)
    
    # plot file
    pngfile <- paste0(out_dir, "/", sprintf("%s_%s_dens.png", jvar, format(idate, "%Y%m%d")))
    png(pngfile, width=3000, height=1750, res=300)
    
    # log-transform
    r = log10(r)
    zlim = c(log10(min), log10(max))
    axis_at = seq(floor(zlim[1]), ceiling(zlim[2]))
    axis_labels = paste(10, '^', axis_at, sep='') 
    
    # plot
    col = rev(brewer.spectral(101))
    plot(box, border="grey20", lwd=5, main=main)
    plot(box, col="white", border=NA, add=TRUE)
    plot(r, maxpixels=1036800, col=col, zlim = zlim, main=main, legend=F, axes=FALSE, box=FALSE, add=TRUE)  # raster plot
    plot(land, col="grey80", border="grey80", add=TRUE)  # land mask
    plot(r, zlim =  zlim, legend.only=TRUE, horizontal = TRUE, col=col, legend.width=1, legend.shrink=0.3,
         axis.args=list(at=axis_at, labels=parse(text=axis_labels), cex.axis=1.2))
    
    dev.off()
  }
  
  
  # calculate average and CV
  s <- mean(s, na.rm=TRUE)
  
  # get common min and max values
  min <- min(minValue(s))
  max <- max(maxValue(s))
  
  # plot title
  main = sprintf("Traffic density - %s - (Average Jan - Jun 2019)", jvar)
  
  # plot file
  pngfile <- paste0(out_dir, "/", sprintf("%s_dens.png", jvar))
  png(pngfile, width=3000, height=1750, res=300)
  
  # log-transform
  r = log10(s)
  zlim = c(log10(min), log10(max))
  axis_at = seq(floor(zlim[1]), ceiling(zlim[2]))
  axis_labels = paste(10, '^', axis_at, sep='') 
  
  # plot
  col = rev(brewer.spectral(101))
  plot(box, border="grey20", lwd=5, main=main)
  plot(box, col="white", border=NA, add=TRUE)
  plot(r, maxpixels=1036800, col=col, zlim = zlim, main=main, legend=F, axes=FALSE, box=FALSE, add=TRUE)  # raster plot
  plot(land, col="grey80", border="grey80", add=TRUE)  # land mask
  plot(r, zlim =  zlim, legend.only=TRUE, horizontal = TRUE, col=col, legend.width=1, legend.shrink=0.3,
       axis.args=list(at=axis_at, labels=parse(text=axis_labels), cex.axis=1.2))
  
  dev.off()
  
}




## NEW plots

library(ggplot2)
library(raster)
library(tidyr)
library(sf)
library(scales)
library(egg)

rasdf <- as.data.frame(s,xy=TRUE)%>%drop_na()


world <- ne_countries(scale = "medium", returnclass = "sf")
box = c(xmin = -40, ymin = 0, xmax = 5, ymax = 45)
land <- st_crop(world, box)

p <- ggplot()+
  geom_raster(aes(x=x,y=y,fill=layer),data=rasdf)+
  geom_sf(fill=grey(0.8),data=land)+
  # scale_fill_viridis_c('Transits month-1 km-2', trans="log10", direction = 1,
  #                      breaks = trans_breaks("log10", function(x) 10^x),
  #                      labels = trans_format("log10", math_format(10^.x)))+
  scale_fill_distiller(palette = "Spectral",
                       name=expression(Traffic~density~(transits~month^-1~km^-2)),
                       trans="log10", direction = -1,
                       breaks = trans_breaks("log10", function(x) 10^x, n=4),
                       labels = trans_format("log10", math_format(10^.x))) +
  coord_sf(expand=c(0,0))+
  labs(x='',y='')+
  theme_article(base_size=16) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_text(size = 15),
        legend.key.height=grid::unit(0.6,"cm"),
        legend.key.width=grid::unit(2,"cm"),
        legend.text = element_text(size = 13),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  # legend title
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, ticks = TRUE))


# Export plot
p_png = paste0(out_dir, "/", sprintf("%s_dens.png", jvar))
ggsave(p_png, p, width=22, height=18, units="cm", dpi=300)
