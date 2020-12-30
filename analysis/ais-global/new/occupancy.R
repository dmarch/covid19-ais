#-------------------------------------------------------------
# occupancy analysis
#-------------------------------------------------------------
# For a given area (global, EEZ, FAO, etc) calculate the occupancy for each ship category at each time stamp.
# We use Mollweide projection to calculate area
# Estimate for both median and hotpost
# Then, calculate difference between equivalent periods in 2019 and 2020
# We could use same approach for T-AIS in June under different spatial resolutions
# to provide a better understanding on how the change of support affect the differences



library(dplyr)
library(raster)
library(lubridate)
library(ggplot2)
library(egg)
library(reshape2)
source("https://raw.githubusercontent.com/dmarch/abigoos/master/R/utils.R")

# set input directory
input_dir <- "data/out/ais-global/density/"

# set output directory
output_dir <- "data/out/ais-global/occupancy/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# select variables to compare
vars <- c("COUNT", "FISHING", "PASSENGER", "CARGO", "TANKER", "OTHER")

# select months to process
dates <- c(
  seq.Date(as.Date("2019-01-01"), as.Date("2019-07-01"), by = "month"),
  seq.Date(as.Date("2020-01-01"), as.Date("2020-07-01"), by = "month")
) #%>% format("%Y%m%d")


## Calculate occupancy

data_list <- list()
cnt <- 1

for (j in 1:length(vars)){
  
  # get ship type
  jvar = vars[j]
  
  for (i in 1:length(dates)){
    
    print(paste("Processing month", i, "from", length(dates)))
    
    # get date
    idate <- dates[i]
    
    # import density raster
    tif_file <- list.files(input_dir, full.names = TRUE, recursive = TRUE, pattern = sprintf("%s_%s_dens_moll.tif$", format(idate, "%Y%m%d"), jvar))
    r1 <- raster(tif_file)
    
    # calculate occupancy (any ship density)
    r1q <- quantile(r1, prob = 0)
    r1ext <- reclassify(r1, c(-Inf, r1q, NA, r1q, Inf, 1))  # create a ocean mask
    occ_km2 <- binsurf(r1ext)  # km2
    
    # calculate occupancy (high density)
    r1q <- quantile(r1, prob = 0.8)
    r1ext <- reclassify(r1, c(-Inf, r1q, NA, r1q, Inf, 1))  # create a ocean mask
    occ_high_km2 <- binsurf(r1ext)  # km2
    
    # create data.frame
    df <- data.frame(var = jvar, date = idate, occ_km2, occ_high_km2)
    data_list[[cnt]] <- df
    cnt <- cnt+1
  }
}

# combine data
data <- data.table::rbindlist(data_list)

# export data
outfile <- paste0(output_dir, "occupancy_global_ais.csv")
write.csv(data, outfile, row.names=FALSE)
data <- read.csv(outfile)
data$date <- ymd(data$date)

# prepare data
data$year <- as.character(year(data$date))
data$month <- month(data$date)

# reorder vessel categories for plots
data$var <- factor(data$var, levels=c("COUNT", "CARGO", "TANKER", "PASSENGER", "FISHING", "OTHER"))

# filter data per study period
data <- filter(data, month <= 6)


## Plot data
p1 <- ggplot(data,  aes(x = month)) + #filter(data, var == "COUNT")
  geom_line(aes(y = occ_high_km2, color = year), size = 1) +
  scale_color_manual(values=c('#9ecae1', "#3182bd"))+
  #scale_x_date(date_breaks = "2 month", date_labels = "%b") +
  scale_x_continuous(breaks = 1:7) +
  ylab("") + xlab("") +
  facet_wrap(var ~ ., ncol = 6)+#, scales = "free") + # , 
  theme_article() +
  theme(legend.position = "none") +
  guides(fill = FALSE)



#-----------------------------------------------------------------
# 2. Calculate percentage change in relation to 2019
#-----------------------------------------------------------------

# transform from long to wide format
dataw <- dcast(data, var + month ~ year, value.var = "occ_km2")

# calculate change metrics
change <- dataw %>%
  mutate(
    delta = `2020` - `2019`,
    per = (delta/`2019`),
    perlog = log(`2020`/`2019`),
    change_positive = per > 0
  ) %>%
  mutate(monthAbb = month.abb[month])

# levels for plot
change$monthAbb <- factor(change$monthAbb, levels=month.abb[1:7])
#change$var <- factor(change$var, levels=c("All vessels", "Cargo", "Tanker", "Fishing", "Other vessels"))
levels(change$var) <- c("All vessels", "Cargo", "Tanker", "Passsenger","Fishing", "Other vessels")

p2 <- ggplot(change, mapping=aes(x = monthAbb, y = per, fill = change_positive)) +
  geom_col(alpha=1, width=0.8) +
  #ylab("Relative change (L%) in occupancy") + xlab("Month") +
  #scale_x_continuous( breaks = c(1,3,5), labels = month.abb[c(1, 3,5)]) +  # , labels = month.abb[1:6]
  #scale_x_continuous(breaks = 1:6) +
  scale_fill_manual(values=c("#e34a33", "#9ecae1")) +
  scale_y_continuous(labels = scales::percent, position = "left")+
  facet_wrap(var ~ ., ncol = 3) +
  xlab("") + ylab("Relative change in occupancy") +
  theme_article() +
  guides(fill = FALSE)

# add tags
#p2 <- tag_facet(p2, open = "", close = "")

# export multi-panel plot
out_file <- paste0(output_dir, "occupancy_relative.png")
ggsave(out_file, p2, width=16, height=10, units = "cm")

#-----------------------------------------------------------------
# 2. Compare count with expected
#-----------------------------------------------------------------

dataw$expected <- dataw$`2019` * 1.0253204  # see "tournadre_occupancy"

# calculate change metrics
change <- dataw %>%
  mutate(
    delta = `2020` - expected,
    per = 100*(delta/expected),
    perlog = 100*log(`2020`/expected),
    change_positive = per > 0
  )


p3 <- ggplot(change, mapping=aes(x = month, y = per, fill = change_positive)) +
  geom_col(alpha=1, width=0.8) +
  #ylab("Relative change (L%) in occupancy") + xlab("Month") +
  #scale_x_continuous( breaks = c(1,3,5), labels = month.abb[c(1, 3,5)]) +  # , labels = month.abb[1:6]
  scale_x_continuous(breaks = 1:6) +
  scale_fill_manual(values=c("#e34a33", "#9ecae1")) +
  facet_wrap(var ~ ., ncol = 6) +
  theme_article() +
  guides(fill = FALSE)




# define layout
lay <- rbind(c(1),
             c(2),
             c(3))

# plot
p <- grid.arrange(p1, p2, p3, layout_matrix = lay)


# export multi-panel plot
out_file <- paste0(output_dir, "occupancy_global_ais.png")
ggsave(out_file, p, width=25, height=15, units = "cm")
