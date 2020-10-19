#---------------------------------------------------------------------------
# stringency            Process and Analyse the stringency index
#---------------------------------------------------------------------------

## load libraries
library(lubridate)
library(ggplot2)
library(dplyr)
library(egg)
library(sf)
library(forcats)
library(tmap)

## set output folder
## here we store plots
output_plot <- "results/stringency"
if (!dir.exists(output_plot)) dir.create(output_plot, recursive = TRUE)

## here we store data
output_data <- "data/out/stringency"
if (!dir.exists(output_data)) dir.create(output_data, recursive = TRUE)

## set last date to process
last_date <- "2020-07-31"


#--------------------------------------------------------------------------------
# 1. Import data
#--------------------------------------------------------------------------------
# There are two formats available

# Option 1: Latest from Github
latest_csv <- "https://raw.githubusercontent.com/OxCGRT/covid-policy-tracker/master/data/OxCGRT_latest.csv"
data <- read.csv(latest_csv)
data$Date <- ymd(data$Date)  # parse date
data$Month <- month(data$Date)

# Option 2: Latest from website
# latest_csv <- "data/input/OxCGRT/covid-stringency-index.csv"
# data <- read.csv(latest_csv)
# data$Date <- mdy(data$Date)
# data$Month <- month(data$Date)
# data$CountryName <- data$Entity
# data$StringencyIndex <- data$Government.Response.Stringency.Index...0.to.100..100...strictest..


# Filter data untill 31st July
data <- filter(data, Date <= as.Date(last_date))

# Remove sub-regions
data <- filter(data, RegionName == "")

# Store a version for further reproducibility
write.csv(data, paste0(output_data, "/stringency.csv"), row.names=FALSE)


#--------------------------------------------------------------------------------
# 2. Plot timeline of Individual countries
#--------------------------------------------------------------------------------

## select countries of interest
WMedCountries <- c("Spain", "France", "Italy", "China", "India", "United States", "South Korea")

# Plot time series all countries in the same plot
p1 <- ggplot(filter(data, CountryName %in% WMedCountries), aes(x = Date, group = CountryName)) +
  geom_line(aes(y = StringencyIndex, color = CountryName), size = 1) +
  geom_vline(xintercept = as.Date("2020-03-11"), linetype="dotted") +
  scale_colour_brewer(palette="Set1") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  ylim(c(0,100))+
  ylab("Stringency Index") + xlab("") +
  theme_article() +
  theme(legend.position = c(0.85, 0.2), legend.title = element_blank())

# Faceted by country
p2 <- ggplot(filter(data, CountryName %in% WMedCountries), aes(x = Date, group = CountryName)) +
  geom_line(aes(y = StringencyIndex, color = CountryName), size = 1) +
  geom_vline(xintercept = as.Date("2020-03-11"), linetype="dotted") +
  scale_colour_brewer(palette="Set2") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  ylab("Stringency Index") + xlab("") +
  facet_wrap(CountryName ~ ., ncol = 3) +
  theme_article() +
  theme(legend.position = "none")

# Barplot
data_sub <-filter(data, CountryName %in% WMedCountries)
data_sub$si_group <- cut(data_sub$StringencyIndex, 5)
data_sub$CountryName <- factor(data_sub$CountryName, levels = c("China", "South Korea", "Italy", "France", "India", "United States", "Spain"))
data_sub$CountryName <- fct_rev(data_sub$CountryName)
p3 <- ggplot(data_sub, aes(Date, CountryName, color = si_group, group=rev(CountryName))) +
  geom_line(size = 10) +
  scale_color_brewer(palette = "RdYlBu", direction=-1, na.value = "grey70") +
  labs(x=NULL, y=NULL) +
  #geom_vline(xintercept = seq.Date(as.Date("2020-01-01"), as.Date("2020-06-30"), "day"), color = "white", lwd=0.01) + 
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  theme_article() +
  theme(legend.position = "none")

# export plot
out_file <- paste0(output_plot, "/country_timeline.png")
ggsave(out_file, p3, width=13, height=8, units = "cm")



#--------------------------------------------------------------------------------
# 2. Plot timeline of average
#--------------------------------------------------------------------------------

# filter coastal countries
eez <- st_read("data/input/marine_regions/World_EEZ_v11_20191118_gpkg/eez_v11.gpkg")
eez_iso <- unique(eez$ISO_SOV1)
data <- data %>% filter(CountryCode %in% eez_iso)
length(unique(data$CountryCode)) # 133 countries with EEZ

# calculate average and SD
dailyAvgGlobal <- data %>%
  group_by(Date) %>%
  filter(Date < as.Date(last_date)) %>%
  summarize(StringencyIndexAvg = mean(StringencyIndex, na.rm=TRUE),
            StringencyIndexSd = sd(StringencyIndex, na.rm=TRUE))

# plot
p <- ggplot(dailyAvgGlobal, aes(x = Date)) +
  geom_rect(aes(xmin=as.Date("2020-04-01"), xmax=as.Date("2020-04-30"), ymin=-Inf, ymax=Inf), fill="grey90", alpha=0.5) +
  geom_ribbon(aes(ymin = StringencyIndexAvg - StringencyIndexSd, ymax = StringencyIndexAvg + StringencyIndexSd), alpha=.2, linetype=0, fill = "steelblue") +
  geom_line(aes(y = StringencyIndexAvg), size = 1, color="steelblue") +  # 'steelblue', "#00BFC4"
  geom_vline(xintercept = as.Date("2020-03-11"), linetype="dotted") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  ylab("Stringency Index") + xlab("") +
  theme_article() +
  guides(fill = FALSE)

# export plot
out_file <- paste0(output_plot, "/daily_global_avg.png")
ggsave(out_file, p, width=10, height=8, units = "cm")

# Calculate average and sd
april <- data %>%
  filter(Date >= as.Date("2020-04-01") & Date < as.Date("2020-04-30"))
mean(april$StringencyIndex, na.rm=TRUE)  # 79.7002
sd(april$StringencyIndex, na.rm=TRUE)  # 15.34479


#--------------------------------------------------------------------------------
# 2. Plot map
#--------------------------------------------------------------------------------

# Calculate median
avgGlobal <- data %>%
  group_by(CountryName, CountryCode) %>%
  filter(Month == 4) %>%
  summarize(StringencyIndexAvg = median(StringencyIndex),
            StringencyIndexSd = sd(StringencyIndex),
            StringencyIndexMin = min(StringencyIndex),
            StringencyIndexMax = max(StringencyIndex))

# import country map
data("World")

# combine with countries with EEZ by ISO code
map_data <- left_join(World, avgGlobal, by = c("iso_a3" = "CountryCode"))
map_data <- map_data %>% filter(name != "Antarctica")

# plot data (average)
p1 <- tm_shape(map_data) +
  tm_polygons("StringencyIndexAvg",
              title=("Stringency Index \nMonthly median (April 2020)"),
              palette = "-RdYlBu", 
              border.col = "grey60", 
              border.alpha = 0.3) +
  tm_layout(legend.title.size=1)

# export plot
out_file <- paste0(output_plot, "/april_global_map_med.png")
tmap_save(tm = p1, filename = out_file, width=22, height=10, units = "cm")
