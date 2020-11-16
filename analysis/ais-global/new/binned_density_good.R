# binned plots


library(raster)
library(dplyr)
library(egg)
library(data.table)
source("https://raw.githubusercontent.com/dmarch/abigoos/master/R/utils.R")

# set input directory
input_dir <- "data/out/ais-global/density/"

# select variables to compare
vars <- c("COUNT", "FISHING", "PASSENGER", "CARGO", "TANKER", "OTHER")


dates_post <- c(
  seq.Date(as.Date("2020-01-01"), as.Date("2020-06-01"), by = "month")
)

dates_pre <- c(
  seq.Date(as.Date("2019-01-01"), as.Date("2019-06-01"), by = "month")
)

dates <- c(dates_pre, dates_post)


# import ocean mask
ocean <- raster("data/out/ais-global/oceanmask_mol.nc")
total_area <- binsurf(ocean)


# get density values for each month and ship cateogory
cnt <- 1
data_list <- list()

for(j in 1:length(vars)){
  
  jvar <- vars[j]
  
  for(i in 1:length(dates)){
    
    idate <- dates[i]
    
    tif_file <- list.files(input_dir, full.names = TRUE, recursive=TRUE, pattern=sprintf("%s_%s_dens_moll.tif$", format(idate, "%Y%m%d"), jvar))
    r <- raster(tif_file)
    
    data_list[[cnt]] <- data.frame(var = jvar, date = idate, dens = na.omit(values(r)))
    cnt <- cnt+1
  }
}

# combine data
data <- rbindlist(data_list)

data$year <- year(data$date)
data$month <- month(data$date)
data$log <- log10(data$dens)


# calculate media per month
med <- data %>%
  group_by(month) %>%
  summarize(med = median(dens))


p <- ggplot(data, aes(month, dens))
p <- p + geom_bin2d(bins = 60,  drop=FALSE)+#(bins=60,  drop=FALSE) +color ="white",
  scale_y_log10(expand = c(0, 0),
                #limits = c(0.001299857, 100),
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_fill_viridis_c(name = "count", trans = "log",
                       breaks = breaks_log(n = 5)) +
  scale_x_continuous(expand = c(0, 0), breaks=c(1:7)) +
  # add line
  geom_line(data = med, aes(x=month, y=med)) +
  theme_article()

expand_scale(add = c(0.2,0))
c(0, 0)
## set output folder
## here we store plots
output_plot <- "results/ais-global/bins"
if (!dir.exists(output_plot)) dir.create(output_plot, recursive = TRUE)



# export figure
out_file <- paste0(output_plot, "/test.png")
ggsave(out_file, p, width=10, height=8, units = "cm")




ggplot(data, aes(x=month, y=log)) +
  stat_bin2d(aes(fill = after_stat(count)), binwidth = c(0.5,0.5))


##### Histograms

# https://cran.r-project.org/web/packages/ggridges/vignettes/introduction.html

library(ggridges)


hist(data$log)
abline(v=log10(med$med))



data$month <- as.factor(data$month)
data$year <- as.factor(data$year)

ggplot(data, aes(x = log, y = month, fill=year)) +
  geom_density_ridges(aes(fill = year), scale = 0.9,
                      quantile_lines = T, quantiles = 2, alpha=.3) +
  #scale_fill_brewer(palette = 4) +
  #facet_wrap(~year) +
  scale_fill_cyclical(values = c("blue", "green")) +
  theme_article() + theme(legend.position = "none")



#### Beanplot

library(beanplot)

#data$month <- format(data$date, "%b")
data$month <- month(data$date)
data$my <- paste(data$month, data$year)
data$my <- as.factor(data$my)

## set output folder
## here we store plots
output_plot <- "results/ais-global/bins"
if (!dir.exists(output_plot)) dir.create(output_plot, recursive = TRUE)

# export figure
out_file <- paste0(output_plot, "/", jvar, "_beans.png")
png(out_file, width=3000, height=1500, res=300)

beanplot(log ~ my, data = data,
         main = "", ylab = "Density (log10)", side = "both",
         border = NA, col = list('#9ecae1', "#3182bd"),
         bw="nrd0",
         overallline = "median",
         what=c(0,1,1,0),
         beanlines="median",
         ll = 0.1,
         maxstripline=0.01,
         beanlinewd=1,
         frame.plot = T,
         boxwex = 0.9#,
         )
#legend("topright", fill = c('#9ecae1', "#3182bd"),legend = c("2019", "2020"))


dev.off()


### Boxplot

data$month <- format(data$date, "%b")
data$month <- factor(data$month, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun"))
data$year <- factor(data$year, levels=c("2019", "2020"))


pd = position_dodge(width = 0.8)

p <- ggplot(filter(data, var != "COUNT"), aes(x= month, y = dens, fill=year)) +
  stat_boxplot(geom="errorbar", position=pd, width=0.2, lwd=0.3) +
  geom_boxplot(width=0.7, position=pd, lwd=0.2, outlier.shape = 21, outlier.alpha = 0.2, outlier.size = 1.5) +
  scale_y_log10(expand = c(0.05, 0.05),
                #limits = c(1e-3, 1e0),
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_fill_manual(values=c('#9ecae1', "#3182bd"))+
  labs(y = "Traffic density", x = "") +
  facet_wrap(var ~ .,  ncol = 1) +
  theme_bw(base_size = 12, base_family = "") +
  theme_article()+
  theme(legend.position =  "none")

p2 <- tag_facet(p, open = "", close = "")


## set output folder
## here we store plots
output_plot <- "results/ais-global/bins"
if (!dir.exists(output_plot)) dir.create(output_plot, recursive = TRUE)

# export multi-panel plot
out_file <- paste0(output_plot, "/boxplot_densities.png")
ggsave(out_file, p2, width=20, height=20, units = "cm")

