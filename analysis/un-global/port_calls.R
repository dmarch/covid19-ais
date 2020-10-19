

# Is there a correlation with traffic density maps from S-AIS globally? (during S-AIS months)
# Is there a correlation with number of active vessels from T-AIS in the Med? (all period available). Select ports using study area.
# Is there a correlation with the Stringency Index? Weekly variability. Plot time series for specific ports.

# Make plot of all ports.


library(data.table)
library(lubridate)
library(dplyr)
library(ggplot2)
library(egg)
library(reshape)
library(sf)

# create output directory
out_dir <- "data/out/un/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


# import data
# original file was encoded in UTF-16. Encoded to UTF-8
# Solution: https://rud.is/b/2019/02/20/i-just-wanted-the-data-turning-tableau-tidyverse-tears-into-smiles-with-base-r-an-encoding-detective-story/
port_call_file <- "data/input/un/Port_Map_Datos_completos_data.csv"
data <- read.csv(port_call_file, sep="\t", dec=",", fileEncoding = "UTF-16LE", stringsAsFactors=FALSE)

# parse date-time
data$Date.Entry <- dmy(data$Date.Entry)
data$month <- floor_date(data$Date.Entry, "month")
data$ym <- format(data$Date.Entry, "%Y%m")

# filter data by date range
data <- filter(data, Date.Entry >= as.Date("2019-01-01"), Date.Entry <= as.Date("2020-07-31"))

# reclass vessel types
unique(data$vessel_type)  # [1] "Container"     "Dry bulk"      "General cargo" "Tanker"        "Vehicles"   
data$type <- "CARGO"
data$type[data$vessel_type == "Tanker"] <- "TANKER"

# select ports (n = 1,189)
# TODO: Summarize and calculate avg of port call per port
ports <- data %>%
          dplyr::select(Port.Name, Lon, Lat, CountryISO) %>%
          distinct()

# export port database
write.csv(ports, paste0(out_dir, "un_ports.csv"), row.names = FALSE)




#-------------------------------------------------------------------------------
# Summarize port calls per port and month


port_monthly <- data %>%
  group_by(Port.Name, Lon, Lat, CountryISO, ym, Year, month) %>%
  summarise(port_calls = sum(Port.Calls))

port_monthly$month2020 <- port_monthly$month
year(port_monthly$month2020) <- 2020

# filter data from January till July
port_monthly <- filter(port_monthly, month2020 <= as.Date("2020-07-01"))

# transform from long to wide format
port_monthly <- dcast(port_monthly, Port.Name + CountryISO + Lon + Lat + month2020 ~Year, value.var="port_calls")

# calculate difference metrics
port_monthly$delta <- port_monthly$'2020' - port_monthly$'2019'
port_monthly$per <- 100 * (port_monthly$delta / port_monthly$'2019')
port_monthly$perlog <- 100 * log(port_monthly$'2020' / port_monthly$'2019')


# transfrom again from long to wide
port_monthly2 <- dcast(port_monthly, Port.Name + CountryISO + Lon + Lat ~ month2020, value.var="perlog")

# export port database
names(port_monthly2)[5:11] <- c("jan", "feb", "mar", "apr", "may", "jun", "jul")
port_sf <- st_as_sf(port_monthly2, coords = c("Lon", "Lat"), crs = 4326)
st_write(port_sf, paste0(out_dir, "un_ports_perlog.gpkg"))
write.table(port_monthly2, paste0(out_dir, "un_ports_perlog.csv"), row.names = FALSE, sep=";", dec=",")


#-------------------------------------------------------------------------------
# Temporal variation of overall port calls

data_global <- data %>%
  group_by(Date.Entry) %>%
  summarise(port_calls = sum(Port.Calls))

# Plot time series all countries in the same plot
p1 <- ggplot(filter(data_global), aes(x = Date.Entry)) +
  geom_line(aes(y = port_calls), size = 0.2) +
  geom_vline(xintercept = as.Date("2020-03-11"), linetype="dotted") +
  scale_colour_brewer(palette="Set1") +
  #scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  #ylim(c(0,100))+
  ylab("Port Calls") + xlab("") +
  theme_article() #+
#theme(legend.position = c(0.85, 0.2), legend.title = element_blank())



#-------------------------------------------------------------------------------
# Monthly variation of overall port calls

data_monthly <- data %>%
  group_by(ym, Year, month) %>%
  summarise(port_calls = sum(Port.Calls))

data_monthly$month2020 <- data_monthly$month
year(data_monthly$month2020) <- 2020

data_monthly$Year <- as.character(data_monthly$Year)

p <- ggplot(filter(data_monthly, month2020 <= as.Date("2020-07-01")), aes(x=month2020, y=port_calls, fill=Year)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual(values=c('#9ecae1', "#3182bd"))+
  ylab("Port Calls") + xlab("Month") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  theme_article() #+



newdata <- dcast(data_monthly, month2020~Year, value.var="port_calls")
newdata$delta <- newdata$'2020' - newdata$'2019'
newdata$per <- 100 * (newdata$delta / newdata$'2019')

#-------------------------------------------------------------------------------


data_week <- data %>%
  group_by(Port.Name, Lon, Lat, type, CountryISO, Date.Entry) %>%
  summarise(port_calls = sum(Port.Calls))

## select ports of interest
sel_ports <- c("Barcelona", "Shanghai", "Rotterdam")

# Plot time series all countries in the same plot
p1 <- ggplot(filter(data_week, Port.Name %in% sel_ports, type == "CARGO"), aes(x = Date.Entry, group = Port.Name)) +
  geom_line(aes(y = port_calls, color = Port.Name), size = 0.2) +
  geom_vline(xintercept = as.Date("2020-03-11"), linetype="dotted") +
  scale_colour_brewer(palette="Set1") +
  #scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  #ylim(c(0,100))+
  ylab("Port Calls") + xlab("") +
  facet_wrap(Port.Name ~ ., ncol = 3, scales = "free") +
  theme_article() #+
  #theme(legend.position = c(0.85, 0.2), legend.title = element_blank())


### Summarize per country 
data_country <- data %>%
  group_by(CountryISO, Country, type, Date.Entry) %>%
  summarise(port_calls = sum(Port.Calls))

sel_countries <- c("Spain", "France", "Italy", "China", "India", "U.S.A.", "South Korea", "United Kingdom", "Germany")

# Plot time series all countries in the same plot
p2 <- ggplot(filter(data_country, Country %in% sel_countries, type == "CARGO"), aes(x = Date.Entry, group = Country)) +
  geom_line(aes(y = port_calls, color = Country), size = 0.2) +
  geom_vline(xintercept = as.Date("2020-03-11"), linetype="dotted") +
  scale_colour_brewer(palette="Set1") +
  #scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  #ylim(c(0,100))+
  ylab("Port Calls") + xlab("") +
  facet_wrap(Country ~ ., ncol = 3, scales = "free") +
  theme_article() #+
#theme(legend.position = c(0.85, 0.2), legend.title = element_blank())





#-------------------------------------------------------------------------------
## Analysis 1. Assess correlation between port calls and traffic density
# Traffic density is monthly and port calls weekly.
# For each port call date, get first day of month. Filter Jan-July 2019-2020. format("%Y%m%d")
# Calculate total number of port calls per month.
# extract density data for each port and each month

# set input directory
input_dir <- "data/out/ais-global/density/"

# get first day of month
# we use it to match S-AIS monthly density maps
data$month <- floor_date(data$Date.Entry, "month") %>% format("%Y%m%d")

# select months to process
dates <- c(
  seq.Date(as.Date("2019-01-01"), as.Date("2019-07-01"), by = "month"),
  seq.Date(as.Date("2020-01-01"), as.Date("2020-06-01"), by = "month")  # ADD Jul 2020 when ready
) %>% format("%Y%m%d")

# select variables to compare
vars <- c("CARGO", "TANKER")


# filter data and calculate total number of port calls per month
data_month <- data %>%
            filter(month %in% dates) %>%
            group_by(Port.Name, Lon, Lat, type, CountryISO, month) %>%
            summarise(port_calls = sum(Port.Calls))
            
# extract density data
# use lon/lat projected data
extract_data <- list()
cnt <- 1

for (i in 1:length(dates)){
  for (j in 1:length(vars)){
  
    print(paste("Processing month", i, "var", j))
  
    # select date and variable
    idate <- dates[i]
    jvar <- vars[j]

    # import raster
    tif <- list.files(input_dir, full.names = TRUE, recursive = TRUE, pattern = sprintf("%s_%s_dens.tif$", idate, jvar))
    r <- raster(tif)
    
    # select port call data for date and variable
    data_sel <- filter(data_month, month == idate, type == jvar)
    
    # extract vessel density
    data_sel$density <- extract(r, cbind(data_sel$Lon, data_sel$Lat))

    # append to list
    extract_data[[cnt]] <- data_sel
    cnt <- cnt + 1
  }
}

# combine all data
data_comb <- rbindlist(extract_data)

plot(data_comb$port_calls, data_comb$density)


# TODO
# OVERLAP IN QGIS THE DENSITY MAP WITH A CSV OF ALL PORTS
# COMBINE PORTS THAT OVERLAP WITHIN A GRID CELL (SUM PORT CALLS)
# EXTRACT PER PORT IS NOT RIGHT APPROACH BECAUSE SPATIAL SCALE OF DENSITY MAPS (E.G. >1 PORT PER CELL)
# INSTEAD, RASTERIZE AND CALCULATE SUM OF VESSELS PER GRID CELL







