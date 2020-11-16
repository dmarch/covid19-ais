#---------------------------------------------------------------------
# MMSI trend plot
#---------------------------------------------------------------------
# Plot the annual trend of MMSI data
#
# Other plots
# Number of merchant fleet (>= 100 gross tons. 98,140)

# Load libraries
library(ggplot2)
library(egg)
library(dplyr)

# import data
# data provided by Simon Chesworth (Exact Earth)
# 'mmsi' total number of AIS transponders with > 20 messages
mmsi <- read.csv("data/input/exacEarth/mmsi_per_year.csv")
mmsi <- mmsi %>%
  rename(num_ships = mmsi) %>%
  mutate(type = "Total AIS")

# import number of merchant ships
ships <- read.csv("data/input/unctad/num_ships.csv")
ships <- ships %>%
      mutate(type = "Merchant")

# import number of fishing ships
# data from 38 countries
fishing <- read.csv("data/input/oecd/FISH_FLEET_15112020011208145.csv")
fishing <- fishing %>%
              group_by(Year) %>%
              summarize(num_ships = sum(Value, na.rm=TRUE)) %>%
              rename(year = Year) %>%
              mutate(type = "All fishing")

# import AIS fishing vessels 2012-2016
gfw <- read.csv("data/input/gfw/fishing-vessels-v1.csv")
gfw <- data.frame(
  year = c(2012:2016),
  num_ships = c(
    sum(gfw$active_2012 == "true"),
    sum(gfw$active_2013 == "true"),
    sum(gfw$active_2014 == "true"),
    sum(gfw$active_2015 == "true"),
    sum(gfw$active_2016 == "true")
  ),
  type = "GFW"
)

# function for y-axis
# Adapted from here: https://stackoverflow.com/questions/10762287/how-can-i-format-axis-labels-with-exponents-with-ggplot2-and-scales
scientific_10 <-  function(x) {
  ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scales::scientific_format()(x)))))
}

# plot
ggplot(data = mmsi, aes(x = year, y = num_ships)) +
  # set column type
  geom_line(color="#3182bd", size=2) + # '#9ecae1', "#3182bd"
  geom_point(size=2, fill="white", shape=21, colour="#3182bd", stroke=2) +
  # add years in x-axis
  scale_x_continuous(breaks=c(2010:2020)) +
  # change format in y-axis
  scale_y_continuous(label=scientific_10) +
  # labels
  xlab("") +
  ylab("Number of AIS transponders") +
  # set theme
  theme_article()



# plot
ggplot(data = ships, aes(x = year, y = num_ships)) +
  # set column type
  geom_line(color="#3182bd", size=2) + # '#9ecae1', "#3182bd"
  geom_point(size=2, fill="white", shape=21, colour="#3182bd", stroke=2) +
  # add years in x-axis
  scale_x_continuous(breaks=c(2010:2020)) +
  # change format in y-axis
  scale_y_continuous(label=scientific_10) +
  # labels
  xlab("") +
  ylab("Number of ships (>= 100 gross tons)") +
  # set theme
  theme_article()


# plot
ggplot(data = filter(fishing, year >= 2010), aes(x = year, y = num_ships)) +
  # set column type
  geom_line(color="#3182bd", size=2) + # '#9ecae1', "#3182bd"
  geom_point(size=2, fill="white", shape=21, colour="#3182bd", stroke=2) +
  # add years in x-axis
  scale_x_continuous(breaks=c(2010:2020)) +
  # change format in y-axis
  scale_y_continuous(label=scientific_10) +
  # labels
  xlab("") +
  ylab("Number of fishing ships") +
  # set theme
  theme_article()


# plot
ggplot(data = gfw, aes(x = year, y = num_ships)) +
  # set column type
  geom_line(color="#3182bd", size=2) + # '#9ecae1', "#3182bd"
  geom_point(size=2, fill="white", shape=21, colour="#3182bd", stroke=2) +
  # add years in x-axis
  scale_x_continuous(breaks=c(2010:2020)) +
  # change format in y-axis
  scale_y_continuous(label=scientific_10) +
  # labels
  xlab("") +
  ylab("Number of fishing ships") +
  # set theme
  theme_article()


### Combine data

data <- rbind(mmsi, ships, gfw)

ggplot(data = data, aes(x = year, y = num_ships, fill = type)) +
  # set column type
  geom_line(size=2) + # '#9ecae1', "#3182bd"
  geom_point(size=2, fill="white", shape=21, stroke=2) +
  # add years in x-axis
  scale_x_continuous(breaks=c(2010:2020)) +
  # change format in y-axis
  scale_y_continuous(label=scientific_10) +
  # labels
  xlab("") +
  ylab("Number of fishing ships") +
  # set theme
  theme_article()
