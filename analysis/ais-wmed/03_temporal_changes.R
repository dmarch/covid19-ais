#---------------------------------------------------------------------------
# 03_temporal_changes                                                        
#---------------------------------------------------------------------------
# This script calculates temporal changes during COVID-19 in relation to
# 2019
#---------------------------------------------------------------------------

## load libraries
library(dplyr)
library(lubridate)
library(stringr)
library(ggplot2)
library(reshape2)
library(forecastML)
library(zoo)
library(janitor)
library(imputeTS)
library(egg)
library(gtable)
library(grid)
library(gridExtra)
library(rworldmap)
library(cowplot)
library(raster)
library(sf)
library(flextable)
source("scr/miniglobe.R")


## Import data
cnt <- read.csv("data/out/ais-wmed/vessels_day.csv")
cnt$date <- as.Date(cnt$date)


## Set output plots
output_data <- "results/ais-wmed"
if (!dir.exists(output_data)) dir.create(output_data, recursive = TRUE)



#-----------------------------------------------------------------
# Calculate averages
#-----------------------------------------------------------------

# Calculate 7-day moving average (and sd)
# Reduces missalignment in day data between years due to influence of weekends
mave <- arrange(cnt, date) %>%
  group_by(type) %>%
  mutate(vessels_ma7 = rollapply(vessels, width = 7, FUN = "mean", align="center", fill=NA),
         vessels_ma30 = rollapply(vessels, width = 30, FUN = "mean", align="center", fill=NA),
         vessels_sd30 = rollapply(vessels, width = 30, FUN = "sd", align="center", fill=NA))#,

# Filter out extremes with no moving averange
mave <- mave[c(first(which(!is.na(mave$vessels_ma7))):last(which(!is.na(mave$vessels_ma7)))),]

# recalculate temporal variables
mave$week <- week(mave$date)
mave$doy <- strftime(mave$date, format = "%j")
mave$year <- as.character(year(mave$date))
mave$year_week <-  paste(mave$year, str_pad(mave$week, 2, pad = "0"), sep="_")
mave$date2020 <- as.Date(as.numeric(mave$doy)-1, origin = "2020-01-01")

# reorder vessel categories for plots
mave$type <- factor(mave$type, levels=c("All", "Cargo", "Tanker", "Passenger", "Fishing", "Recreational", "Other"))



#-----------------------------------------------------------------
# 1. Plot number of vessels per day during 2020
#-----------------------------------------------------------------

last_doy <- 365  # 15th May 2020
mave_sel <- filter(mave, doy <= last_doy, year >= 2019)

# create ribbon
rib <- mave_sel %>%
  group_by(date2020, type) %>%
  filter(date2020 <= max(date)) %>%
  summarize(ymin = min(vessels_ma7),
            ymax = max(vessels_ma7))


# plot all vessels
p1 <- ggplot(filter(mave_sel, type=="All"), aes(x = date2020)) + # removed group=year
  geom_ribbon(data = filter(rib, type=="All"), aes(ymin = ymin, ymax = ymax), fill="#3182bd", alpha=.2, linetype=0) +
  geom_line(aes(y = vessels_ma7, color = year), size = 1) +
  geom_vline(xintercept = as.Date("2020-03-11"), linetype="dotted") +
  scale_color_manual(values=c('#9ecae1', "#3182bd"))+
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  ylab("") + xlab("") +
  annotate(geom = "text", x = first(mave_sel$date2020), y = Inf, label = "a", hjust = 1.2, vjust = 2, fontface = 2, family = "") +
  theme_article() +
  theme(legend.position = c(0.9, 0.9), legend.title = element_blank()) +
  guides(fill = FALSE)

# plot by type
p2 <- ggplot(filter(mave_sel, type != "All"), aes(x = date2020)) +
  geom_ribbon(data = filter(rib, type !="All"), aes(ymin = ymin, ymax = ymax), fill="#3182bd", alpha=.2, linetype=0) +
  geom_line(aes(y = vessels_ma7, color = year), size = 1) +
  geom_vline(xintercept = as.Date("2020-03-11"), linetype="dotted") +
  scale_color_manual(values=c('#9ecae1', "#3182bd"))+
  scale_x_date(date_breaks = "2 month", date_labels = "%b") +
  ylab("") + xlab("") +
  facet_wrap(type ~ ., ncol = 2, scales = "free") +
  theme_article() +
  theme(legend.position = "none") +
  guides(fill = FALSE)

# add letters
p2 <- tag_facet(p2, x = first(mave_sel$date2020), y = -Inf, 
                vjust = -10.5, hjust=0.5,
                open = "", close = "",
                tag_pool = letters[-1])

# define layout
lay <- rbind(c(1,1),
             c(2,2),
             c(2,2))

# plot
p <- grid.arrange(p1, p2, layout_matrix = lay, left = textGrob("Number of vessels", rot = 90, vjust = 1))

# export multi-panel plot
out_file <- paste0(output_data, "/daily_average.png")
ggsave(out_file, p, width=20, height=20, units = "cm")





#-----------------------------------------------------------------
# 2. Calculate percentage change in relation to 2019
#-----------------------------------------------------------------

# calculate percent change
change <- mave_sel %>%
  filter(date2020 <= last(date)) %>%
  group_by(date2020, type) %>%
  summarize(change_per = ((last(vessels_ma7) - first(vessels_ma7))  / first(vessels_ma7))*100,
            change_perlog = 100*log((last(vessels_ma7))  / first(vessels_ma7)),
            change_positive = change_per > 0)

# calculate maximum decline per sector
max_decline <- change %>%
  filter(date2020 > as.Date("2020-03-11")) %>%
  group_by(type) %>%
  summarize(max_drop = min(change_per, na.rm=TRUE),
            median_drop = median(change_per, na.rm=TRUE),
            last_drop = last(change_per)) %>%
  merge(., change, by = 'type') %>%
  filter(max_drop == change_per) %>%
  dplyr::select(type, max_drop, date2020, median_drop, last_drop) %>%
  rename(date_max = date2020) %>%
  mutate(date_last = last(change$date2020))


#-----------------------------------------------------------------
# 2.1. Generate table
#-----------------------------------------------------------------

# set names for categories
vars_names <- c("Cargo", "Tanker", "Passenger", "Fishing", "Recreational", "Other", "All")
vars_renames <- c("Cargo", "Tanker", "Passenger", "Fishing", "Recreational", "Other", "All vessels")

# reorder according to types
max_decline$type <- factor(max_decline$type, levels = vars_names)
max_decline <- arrange(max_decline, type)

# rename categories
max_decline$type <- as.character(max_decline$type)
max_decline$type <- vars_renames

# reorder vessel categories for plots
max_decline$type <- factor(max_decline$type, levels=c("All vessels", "Cargo", "Tanker", "Passenger", "Fishing", "Recreational", "Other"))


# create flextable
ft <- flextable(max_decline)  # data is a data.frame

# set labels for headers (data.frame columns)
# now only those labels without complex text. To break lines, use "\n"
ft <- set_header_labels(ft, type = "Vessel category",
                        max_drop = "(%)",
                        date_max = "Date",
                        last_drop = "(%)",
                        median_drop = "(%)",
                        date_last = "Date")

# add rows on top to include groups
# also include vertical lines to separate groups
ft <- add_header_row(ft, values = c("", "Maximal", "Median", "Most recent"), colwidths = c(1,2,1,2))
ft <- theme_booktabs(ft)

# set decimal formats
ft <- colformat_num(ft, j = c(2,4), digits = 1)

# general configuration of table
ft <- bold(ft, part = "header") # bold header
ft <- fontsize(ft, part = "all", size = 11)  # font size
ft <- set_table_properties(ft, width = 0.4, layout = "autofit")  # autofit to width

# align columns
ft <- align(ft, i = NULL, j = c(2:5), align = "center", part = "all")  # center columns except vessel category

# export table for Word
save_as_docx(ft, path = paste0(output_data, "/maxdrop.docx"))

# export as image
img_file <- paste0(output_data, "/maxdrop.png")
save_as_image(ft, path = img_file)


#-----------------------------------------------------------------
# 2.2. Plot changes
#-----------------------------------------------------------------

# plot (all vessels)
p1 <- ggplot(filter(change, type=="All")) +
  geom_rect(aes(xmin=as.Date("2020-03-11"), xmax=as.Date("2020-06-22"), ymin=-Inf, ymax=Inf), fill="grey90", alpha=0.5) +
  geom_col(aes(x = date2020, y = change_per, fill = change_positive), alpha=1, width=1) +
  geom_point(data = filter(max_decline, type=="All vessels"), aes(x=date_max, y=max_drop), shape=1, size=3, stroke=1)+
  ylab("") + xlab("") +
  #geom_vline(xintercept = as.Date("2020-03-11"), linetype="dotted") +
  scale_fill_manual(values=c("#e34a33", "#9ecae1")) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  geom_text(data = filter(max_decline, type=="All vessels"), aes(x=as.Date("2020-10-11"), y=-60,label=paste(format(round(max_drop, digits=1), nsmall = 1) , "%")))+
  geom_text(data = filter(max_decline, type=="All vessels"), aes(x=as.Date("2020-05-01"), y=15, label=paste(format(round(median_drop, digits=1), nsmall = 1), "%")))+
  theme_article() +
  guides(fill = FALSE)


# plot (categories)
p2 <- ggplot(filter(change, type!="All")) +
  geom_rect(aes(xmin=as.Date("2020-03-11"), xmax=as.Date("2020-06-22"), ymin=-Inf, ymax=Inf), fill="grey90", alpha=0.5) +
  geom_col(aes(x = date2020, y = change_per, fill = change_positive), alpha=1, width=1) +
  geom_point(data = filter(max_decline, type!="All vessels"), aes(x=date_max, y=max_drop), shape=1, size=3, stroke=1)+
  ylab("") + xlab("") +
  #geom_vline(xintercept = as.Date("2020-03-11"), linetype="dotted") +
  scale_x_date(date_breaks = "2 month", date_labels = "%b") +
  scale_fill_manual(values=c("#e34a33", "#9ecae1")) +
  geom_text(data = filter(max_decline, type!="All vessels"), aes(x=as.Date("2020-09-11"), y=-80, label=paste(format(round(max_drop, digits=1), nsmall = 1) , "%")))+
  geom_text(data = filter(max_decline, type!="All vessels"), aes(x=as.Date("2020-05-01"), y=100, label=paste(format(round(median_drop, digits=1), nsmall = 1), "%")))+
  facet_wrap(type ~ ., ncol = 2) +
  theme_article() +
  guides(fill = FALSE)


p2 <- tag_facet(p2, x = first(mave_sel$date2020), y = -Inf, 
                vjust = -12, hjust=0,
                open = "", close = "",
                tag_pool = letters[-1])


lay <- rbind(c(1,1),
             c(2,2),
             c(2,2))

p <- grid.arrange(p1, p2, layout_matrix = lay, left = textGrob("Relative change (%) in number of vessels", rot = 90, vjust = 1))


# export multi-panel plot
out_file <- paste0(output_data, "/daily_change.png")
ggsave(out_file, p, width=20, height=20, units = "cm")



#-----------------------------------------------------------------
# 3. plot number of vessels per day since 2016
#-----------------------------------------------------------------

mave_sel <- mave

# plot all vessels
p1 <- ggplot(filter(mave_sel, type=="All"), aes(x = date, group = year)) +
  geom_ribbon(aes(ymin = vessels_ma30 - vessels_sd30, ymax = vessels_ma30 + vessels_sd30), fill="#3182bd", alpha=.2, linetype=0) +
  geom_line(aes(y = vessels_ma30), size = 1, colour="#3182bd") +
  geom_vline(xintercept = as.Date("2020-03-11"), linetype="dotted") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  ylab("") + xlab("") +
  annotate(geom = "text", x = first(mave_sel$date), y = Inf, label = "a", hjust = 1.2, vjust = 2, fontface = 2, family = "") +
  theme_article() +
  theme(legend.position = c(0.9, 0.9), legend.title = element_blank()) +
  guides(fill = FALSE)


# plot by type
p2 <- ggplot(filter(mave_sel, type != "All"), aes(x = date, group = type)) +
  geom_ribbon(aes(ymin = vessels_ma30 - vessels_sd30, ymax = vessels_ma30 + vessels_sd30), fill="#3182bd", alpha=.2, linetype=0) +
  geom_line(aes(y = vessels_ma30), colour="#3182bd", size = 0.5) +
  geom_vline(xintercept = as.Date("2020-03-11"), linetype="dotted") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  ylab("") + xlab("") +
  facet_wrap(type ~ ., scales = "free", ncol = 2) +
  theme_article() +
  theme(legend.position = "none") +
  guides(fill = FALSE)

p2 <- tag_facet(p2, x = first(mave_sel$date), y = -Inf, 
                vjust = -10, hjust=0,
                open = "", close = "",
                tag_pool = letters[-1])

lay <- rbind(c(1,1),
             c(2,2),
             c(2,2))

p <- grid.arrange(p1, p2, layout_matrix = lay, left = textGrob("Number of vessels", rot = 90, vjust = 1))

# export multi-panel plot
out_file <- paste0(output_data, "/daily_average_allperiod.png")
ggsave(out_file, p, width=20, height=20, units = "cm")



#-----------------------------------------------------------------
# 4. Globe map with study area
#-----------------------------------------------------------------

bb <- st_as_sfc(st_bbox(c(xmin = -1.60, xmax = 9.93, ymax = 35.56, ymin = 43.45), crs = st_crs(4326)))
p <- miniglobe(lon=0, lat=25, pol=bb)
out_file <- paste0(output_data, "/globe.png")
ggsave(out_file, p, width=10, height=10, units = "cm")





