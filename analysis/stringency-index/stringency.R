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
p1 <- ggplot(data, aes(x = Date, group = CountryName)) +
  geom_line(aes(y = StringencyIndex), size = 1, color="grey60", alpha=0.5) +
  geom_vline(xintercept = as.Date("2020-03-11"), linetype="dotted") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  ylim(c(0,100))+
  ylab("Stringency Index") + xlab("") +
  theme_article() +
  theme(legend.position = c(0.85, 0.2), legend.title = element_blank())


# Plot time series all SELECTED countries in the same plot
p2 <- ggplot(filter(data, CountryName %in% WMedCountries), aes(x = Date, group = CountryName)) +
  geom_line(aes(y = StringencyIndex, color = CountryName), size = 1) +
  geom_vline(xintercept = as.Date("2020-03-11"), linetype="dotted") +
  scale_colour_brewer(palette="Set1") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  ylim(c(0,100))+
  ylab("Stringency Index") + xlab("") +
  theme_article() +
  theme(legend.position = c(0.85, 0.2), legend.title = element_blank())

# Faceted by country
p3 <- ggplot(filter(data, CountryName %in% WMedCountries), aes(x = Date, group = CountryName)) +
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
data_sub$si_group <- cut(data_sub$StringencyIndex, 5) # 5
data_sub$CountryName <- factor(data_sub$CountryName, levels = c("China", "South Korea", "Italy", "France", "India", "United States", "Spain"))
data_sub$CountryName <- fct_rev(data_sub$CountryName)
p4 <- ggplot(data_sub, aes(Date, CountryName, color = si_group, group=rev(CountryName))) +
  geom_line(size = 10.5, color="grey20") +
  geom_line(size = 10) +
  scale_color_brewer(palette = "YlGnBu", direction=1, na.value = "grey70") + #brewer.ylgnbu(10). RdYlBu
  labs(x=NULL, y=NULL) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b", expand = c(0,0)) +
  theme_article() +
  theme(legend.position = "none")

# export plot
out_file <- paste0(output_plot, "/country_timeline.png")
ggsave(out_file, p4, width=11, height=8, units = "cm")


# Heat map
library(RColorBrewer)
library(viridis)
library(pals)
textcol <- "grey40"


# for each country calculate the first date when they start increasing stringency index
first_increase <- data %>%
                    filter(StringencyIndex > 50) %>%
                    group_by(CountryName) %>%
                    summarize(first_date = first(Date)) %>% 
                    arrange(first_date)
                    

data_sub <-filter(data, CountryName %in% first_increase$CountryName)
data_sub$CountryName <- factor(data_sub$CountryName, levels = first_increase$CountryName)
data_sub$CountryName <- fct_rev(data_sub$CountryName)

# create discrete categories
data_sub$si_group <- cut(data_sub$StringencyIndex, 10)
data_sub <-filter(data_sub, CountryName %in% WMedCountries)


data_sub <- data_sub %>%
        mutate(si_group=cut(StringencyIndex, breaks=c(-1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100),
                            labels=seq(10,100,10)))


m4 <- m3 %>%
  # convert state to factor and reverse order of levels
  mutate(state=factor(state,levels=rev(sort(unique(state))))) %>%
  # create a new variable from count
  mutate(countfactor=cut(count,breaks=c(-1,0,1,10,100,500,1000,max(count,na.rm=T)),
                         labels=c("0","0-1","1-10","10-100","100-500","500-1000",">1000"))) %>%
  # change level order
  mutate(countfactor=factor(as.character(countfactor),levels=rev(levels(countfactor))))



#basic ggplot
p5 <- ggplot(data_sub, aes(x = Date, y = CountryName, fill = si_group))+
  #add border white colour of line thickness 0.25
  geom_tile(colour="grey90",size=0.1)+
  #remove extra space
  scale_y_discrete(expand=c(0,0))+
  #define new breaks on x-axis
  scale_x_date(expand=c(0,0), date_breaks = "1 month", date_labels = "%b") +
  #colors
  scale_fill_manual(values=brewer.ylgnbu(5), na.value="grey90")+
  #set a base size for all fonts
  theme_grey(base_size=8)+
  # labels and title
  labs(x="",y="",title="Stringency index per country")+
  # title legend
  #guides(fill=guide_legend(label.position = "bottom", title="Stringency Index"))+
  guides(fill = guide_legend(title="Stringency Index",
                             reverse = FALSE,
                             title.position = "top",
                             label.position = "bottom",
                             keywidth = 3,
                             nrow = 1)) +
  coord_fixed(ratio = 20) +
  #coord_equal() +
  theme(legend.position="bottom",legend.direction="horizontal",
        legend.spacing.x = unit(0, 'cm'),
      legend.title=element_text(colour=textcol),
      legend.margin=margin(grid::unit(0,"cm")),
      legend.text=element_text(colour=textcol,size=7,face="bold"),
      legend.key.height=grid::unit(0.8,"cm"),
      legend.key.width=grid::unit(0.2,"cm"),
      axis.text.x=element_text(size=10,colour=textcol),
      axis.text.y=element_text(vjust=0.2,colour=textcol),
      axis.ticks=element_line(size=0.4),
      plot.background=element_blank(),
      panel.border=element_blank(),
      plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),
      plot.title=element_text(colour=textcol,hjust=0,size=14,face="bold"))

# export figure
#out_file <- paste0(output_plot, "/heatmap.png")
#ggsave(out_file, p5, width=13, height=8, units = "cm")



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
  #geom_rect(aes(xmin=as.Date("2020-04-01"), xmax=as.Date("2020-04-30"), ymin=-Inf, ymax=Inf), fill="grey90", alpha=0.5) +
  geom_ribbon(aes(ymin = StringencyIndexAvg - StringencyIndexSd, ymax = StringencyIndexAvg + StringencyIndexSd), alpha=.2, linetype=0, fill = "steelblue") +
  geom_line(aes(y = StringencyIndexAvg), size = 1, color="steelblue") +  # 'steelblue', "#00BFC4"
  geom_vline(xintercept = as.Date("2020-03-11"), linetype="dotted") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b", expand = c(0,0)) +
  ylab("Stringency Index") + xlab("") +
  theme_article() +
  guides(fill = FALSE)

# export plot
out_file <- paste0(output_plot, "/daily_global_avg.png")
ggsave(out_file, p, width=10, height=8, units = "cm")

# Calculate average and sd
april <- data %>%
  filter(Date >= as.Date("2020-04-01") & Date < as.Date("2020-04-30"))
mean(april$StringencyIndex, na.rm=TRUE)  # 80.00187
sd(april$StringencyIndex, na.rm=TRUE)  # 14.78759


#--------------------------------------------------------------------------------
# 2. Plot map for one month
#--------------------------------------------------------------------------------

# import country map
# data("World")
# 
# # Calculate median
# avgGlobal <- data %>%
#   group_by(CountryName, CountryCode) %>%
#   filter(Month == 4) %>%
#   summarize(StringencyIndexAvg = median(StringencyIndex),
#             StringencyIndexSd = sd(StringencyIndex),
#             StringencyIndexMin = min(StringencyIndex),
#             StringencyIndexMax = max(StringencyIndex))
# 
# 
# # combine with countries with EEZ by ISO code
# map_data <- left_join(World, avgGlobal, by = c("iso_a3" = "CountryCode"))
# map_data <- map_data %>% filter(name != "Antarctica")
# 
# # plot data (average)
# p1 <- tm_shape(map_data) +
#   tm_polygons("StringencyIndexAvg",
#               title=("Stringency Index \nMonthly median (April 2020)"),
#               palette = "YlGnBu", #brewer.ylgnbu(10) # YlGnBu, -RdYlBu
#               border.col = "grey60", 
#               border.alpha = 0.3) +
#   tm_layout(legend.title.size=1)

# export plot
#out_file <- paste0(output_plot, "/april_global_map_med.png")
#tmap_save(tm = p1, filename = out_file, width=22, height=10, units = "cm")



#--------------------------------------------------------------------------------
# 2. Plot map for all months
#--------------------------------------------------------------------------------

source("scr/fun_common.R")  # bb()
box <- bb(xmin = -180, xmax = 180, ymin = -90, ymax = 90, crs="+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

# import country map
data("World")
#data(countriesHigh, package = "rworldxtra", envir = environment())
#World <- st_as_sf(countriesHigh)

for(i in 1:7){
  
  # Calculate median
  avgGlobal <- data %>%
    group_by(CountryName, CountryCode) %>%
    filter(Month == i) %>%
    summarize(StringencyIndexAvg = mean(StringencyIndex),
              StringencyIndexSd = sd(StringencyIndex),
              StringencyIndexMin = min(StringencyIndex),
              StringencyIndexMax = max(StringencyIndex))
  
  
  # combine with countries with EEZ by ISO code
  map_data <- left_join(World, avgGlobal, by = c("iso_a3" = "CountryCode"))
  map_data <- map_data %>% filter(name != "Antarctica")
  
  # plot data (average)
  p1 <- tm_shape(box) +
    tm_borders()+
    tm_shape(map_data) +
    tm_polygons("StringencyIndexAvg",
                title=(paste0("Stringency Index \nMonthly median (",month.abb[i], " 2020)")),
                palette = "YlGnBu", # brewer.ylgnbu(10)
                breaks = seq(0,100,20),
                border.col = "grey10", 
                border.alpha = 0.3,
                legend.show = TRUE,
                legend.is.portrait = F) +
    tm_layout(title = paste(month.abb[i], 2020), title.size = 2, title.position = c("center","bottom"),
              frame = F,
              legend.show=FALSE,
              legend.outside = TRUE, legend.outside.position = "bottom", legend.title.size=1)

  
  # export plot
  out_file <- paste0(output_plot, paste0("/", month.abb[i], "_global_avg.png"))
  tmap_save(tm = p1, filename = out_file, width=22, height=10, units = "cm")
}


# plot with lengend
p1 <- tm_shape(box) +
  tm_borders()+
  tm_shape(map_data) +
  tm_polygons("StringencyIndexAvg",
              title=(paste0("Stringency Index (Monthly mean)")),
              palette = "YlGnBu", # brewer.ylgnbu(10)
              breaks = seq(0,100,20),
              border.col = "grey10", 
              border.alpha = 0.3,
              legend.show = TRUE,
              legend.is.portrait = T) +
  tm_layout(#title = paste(month.abb[i], 2020), title.size = 2, title.position = c("center","bottom"),
            frame = F,
            legend.show=TRUE,
            legend.outside = TRUE, legend.outside.position = "bottom",
            legend.title.size=1, legend.text.size = 1)

# export plot
out_file <- paste0(output_plot, paste0("/", month.abb[i], "_global_avg_legend.png"))
tmap_save(tm = p1, filename = out_file, width=22, height=10, units = "cm")
