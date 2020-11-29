#-----------------------------------------------------------------------------
# model_count_EEZ       Model count data per EEZ
#-----------------------------------------------------------------------------
# For each month amd country, add Stringency Index
# Issues found in the data: too large dataset. Switch to density median with log-transform data.

library(dplyr)
library(lubridate)

#---------------------------------------------
# Import and prepare count data 
#---------------------------------------------
# import data
data <- read.csv("data/out/ais-global/eez/eez_ais_count.csv")

# filter by ship type
jvar <- "PASSENGER"
jdata <- filter(data, var == jvar)
rm(data)

# parse data
jdata$date <- parse_date_time(jdata$date, "Ymd")

# select main territories (TERRITORY1 name = SOVEREIGN1 name)
jdata <- jdata %>%
  mutate(TERRITORY1 = as.character(TERRITORY1),
         SOVEREIGN1 = as.character(SOVEREIGN1)) %>%
  filter(TERRITORY1 == SOVEREIGN1)



#---------------------------------------------
# Import and prepare Stringency Idenx
#---------------------------------------------
# import Stringency index
# see stringency.R
si <- read.csv("data/out/stringency/stringency.csv")

# Calculate median
si_sum <- si %>%
  group_by(CountryCode, Month) %>%
  summarize(SiMed = median(StringencyIndex)) %>%
  rename(date = Month) %>%
  mutate(date = as.Date(paste(2020, date, 1, sep="-")))


#---------------------------------------------
# Combine Count data and Stringency Idenx
#---------------------------------------------

# combine with countries with EEZ by ISO code
# first, filter countries without SI
# then combine data by date and country code
jdata2 <- jdata %>%
  filter(ISO_SOV1 %in% unique(si_sum$CountryCode)) %>%
  left_join(si_sum, by = c("ISO_SOV1" = "CountryCode", "date" = "date"))

# Create a binary variable to define lockdown
jdata2 <- jdata2 %>%
  mutate(SiMed = replace(SiMed, is.na(SiMed), 0)) %>%
  mutate(covid = ifelse(SiMed > 0, "after", "before"))

# create time
jdata2 <- jdata2 %>%
  mutate(time = as.numeric(difftime(jdata2$date, min(jdata2$date), units="days")))


#---------------------------------------------
# Model
#---------------------------------------------
# ZIP Poisson, using count data. Negative binomial if overdispersed
# Offset: eez_km2, to control differences in cell size
# Random effects: ISO_SOV1, Month
# Fixed effect: Year, SiMed, covid, Year * covid
# DiD

library(mgcv)

jdata2$covid <- as.factor(jdata2$covid)
jdata2$covid <- relevel(jdata2$covid, "before")
jdata2$ISO_SOV1 <- as.factor(jdata2$ISO_SOV1)


m7 <- lme(data=jdata2, vessels ~ covid, corr=corARMA(form=~1|ISO_SOV1/time, q=2), random= ~1|ISO_SOV1)
m1 <- gam(data=jdata2, vessels ~ covid, family = ziP)


