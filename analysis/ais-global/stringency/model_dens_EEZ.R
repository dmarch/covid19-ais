#-----------------------------------------------------------------------------
# model_dens_EEZ       Model density data per EEZ
#-----------------------------------------------------------------------------
# For each month amd country, add Stringency Index


library(dplyr)
library(lubridate)

#---------------------------------------------
# Import and prepare density data 
#---------------------------------------------
# import data
data <- read.csv("data/out/ais-global/eez/eez_ais_density.csv")

# parse date and derive month-year
data$date <- parse_date_time(data$date, "Ymd")
data$month <- month(data$date)
data$year <- year(data$date)

# select main territories (TERRITORY1 name = SOVEREIGN1 name)
data <- data %>%
  mutate(TERRITORY1 = as.character(TERRITORY1),
         SOVEREIGN1 = as.character(SOVEREIGN1)) %>%
  filter(TERRITORY1 == SOVEREIGN1)

# filter by ship type
jvar <- "PASSENGER"
jdata <- filter(data, var == jvar)
#rm(data)






#---------------------------------------------
# Import and prepare Stringency Idenx
#---------------------------------------------
# import Stringency index
# see stringency.R
si <- read.csv("data/out/stringency/stringency.csv")

# Calculate median
si_sum <- si %>%
  group_by(CountryCode, Month) %>%
  summarize(SiAvg = mean(StringencyIndex)) %>%
  rename(month = Month) %>%
  mutate(date = as.Date(paste(2020, month, 1, sep="-")))


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
  mutate(SiAvg = replace(SiAvg, is.na(SiAvg), 0)) %>%  ## CARE: FIRST REPLACE THOSE NA IN 2020. NOT ALL COUNTRIES HAVE DATA.
  mutate(covid = ifelse(SiAvg > 20, "after", "before"))

# create time
jdata2 <- jdata2 %>%
  mutate(time = as.numeric(difftime(jdata2$date, min(jdata2$date), units="days")))



#---------------------------------------------
# Prepare response variable
#---------------------------------------------
# Median density per EEZ and date

jdata3 <- jdata2 %>%
  group_by(ISO_SOV1, TERRITORY1, date, SiAvg, covid, time) %>%
  filter(vessels_km2 > 0) %>%
  summarize(medDens = median(vessels_km2),
            avgDens = mean(vessels_km2),
            surface = sum(eez_km2),
            ncells = n())

jdata3$month <- month(jdata3$date)
jdata3$year <- year(jdata3$date)
jdata3$avgDensLog <- log10(jdata3$avgDens)

#---------------------------------------------
# Model
#---------------------------------------------
# ZIP Poisson, using count data. Negative binomial if overdispersed
# Offset: eez_km2, to control differences in cell size
# Random effects: ISO_SOV1, Month
# Fixed effect: Year, SiMed, covid, Year * covid
# DiD

library(mgcv)
library(lme4)

jdata3$covid <- as.factor(jdata3$covid)
jdata3$covid <- relevel(jdata3$covid, "before")
jdata3$ISO_SOV1 <- as.factor(jdata3$ISO_SOV1)



m1 <- glm(vessels ~ covid, family = gaussian, data = jdata3)
m2 <- lme(data=jdata3, vessels ~ covid, corr=corARMA(form=~1|ISO_SOV1/time, q=2), random= ~1|ISO_SOV1)
m3 <- lme(data=jdata3, vessels ~ SiMed, corr=corARMA(form=~1|ISO_SOV1/time, q=2), random= ~1|ISO_SOV1)
m4 <- lme(data=jdata3, vessels ~ SiMed * month, corr=corARMA(form=~1|ISO_SOV1/time, q=2), random= ~1|ISO_SOV1)

m5 <- lme(data=jdata3, vessels ~ SiMed, random = ~ 1|ISO_SOV1|time)

m6 <- lme(data=jdata3, vessels ~ year * month, random= ~1|ISO_SOV1)
m6 <- lme(data=jdata3, vessels ~ year:month)
m7 <- lme(data=jdata3, vessels ~ SiMed * month, random= ~1|ISO_SOV1)
m8 <- lme(data=jdata3, vessels ~ covid * month, random= ~1|ISO_SOV1)

library(sjPlot)
library(sjmisc)
library(ggplot2)
require(lattice)  
library(ggeffects)
plot_model(m6, type = "pred", terms = c("month", "year"))
p <- ggpredict(m6, terms = c("month", "year"))
plot(p)
qqmath(m6, id=0.05)


ggCaterpillar(ranef(n6, condVar=TRUE))  ## using ggplot2
qqmath(ranef(m6, condVar=TRUE))  


m5 <- lme(data=jdata3, vessels ~ SiMed * month, random = ~ 1|ISO_SOV1)
m6 <- lme(data=jdata3, vessels ~ covid * month, random = ~ 1|ISO_SOV1)


m1 <- glm(data=jdata3, vessels ~ covid*month)
m2 <- glm(data=jdata3, vessels ~ SiMed + ISO_SOV1)
plot_model(m1, type = "pred", terms = c("month", "covid"))
plot_model(m2, type = "pred", terms = c("SiMed"))
plot_model(m5, type = "pred", terms = c("SiMed", "month"))
plot_model(m6, type = "pred", terms = c("covid", "month"), show.ci = FALSE)


p <- ggpredict(m5, terms = c("SiMed", "month"))
plot(p)

p <- ggpredict(m5, terms = c("month", "SiMed"))
plot(p)

p <- ggpredict(m6, terms = c("month", "covid"))
plot(p)


p <- ggpredict(m2, terms = c("SiMed"))
plot(p)


## Caterpillar plot
lattice::dotplot(ranef(m5, postVar = TRUE))



library(pdp)

# Compute partial dependence data for lstat and rm
pd <- partial(m5, pred.var = c("SiMed", "month"))

# Default PDP
pdp1 <- plotPartial(pd)

# Add contour lines and use a different color palette
rwb <- colorRampPalette(c("red", "white", "blue"))
pdp2 <- plotPartial(pd, contour = TRUE, col.regions = rwb)

# 3-D surface
pdp3 <- plotPartial(pd, levelplot = FALSE, zlab = "vessels", colorkey = TRUE, 
                    screen = list(z = -20, x = -60))

# Figure 5
grid.arrange(pdp1, pdp2, pdp3, ncol = 3)




#---------------------------------------------
# Model al grid level
#---------------------------------------------

jdata2$month <- month(jdata2$date)
m1 <- glmer(data=jdata2, vessels ~ covid * month + (1|ISO_SOV1/idcell), offset=log(eez_km2), family = poisson(link = "log"))

# rescale variables?


# refs
# https://www.dallasfed.org/~/media/documents/institute/wpapers/2020/0384.pdf
# https://www.nature.com/articles/s41598-020-75848-2   3d plots




#---------------------------------------------
# Model differences
#---------------------------------------------
# change from baseline as a dependent variable in mixed models, specifically mixed models for repeat measurements (MMRM).

# transform from long to wide format
library(reshape2)
wide <- jdata3 %>%
  dcast(ISO_SOV1 + TERRITORY1 + month ~ year, value.var=c("avgDens")) %>%
  left_join(si_sum, by = c("ISO_SOV1" = "CountryCode", "month" = "month")) %>%
  filter(!is.na(SiAvg)) %>%
  rename(dens2019 = `2019`, dens2020 = `2020`) %>%
  mutate(delta = dens2020 - dens2019,
         per = 100*(delta/dens2019))


# model
m1 <- glm(delta ~ SiAvg, family = gaussian, data = wide)
p <- ggpredict(m1, terms = c("SiAvg"))
plot(p)

m2 <- glm(dens2020 ~ SiAvg + dens2019, family = gaussian, data = wide)
p <- ggpredict(m2, terms = c("SiAvg"))
plot(p)

m3 <- glm(delta ~ SiAvg + ISO_SOV1, family = gaussian, data = wide)
p <- ggpredict(m3, terms = c("SiAvg"))
plot(p)


m4 <- glm(per ~ SiAvg + dens2019, family = gaussian, data = wide)
p <- ggpredict(m4, terms = c("SiAvg"))
plot(p)


library(lmerTest)
m3 <- lmer(data=wide, delta ~ poly(SiAvg, 2) + (1|ISO_SOV1))
p <- ggpredict(m3, terms = c("SiAvg [all]"))
plot(p)

m4 <- lmer(data=wide, delta ~ SiAvg + (1|ISO_SOV1))
p <- ggpredict(m4, terms = c("SiAvg"))
plot(p)

m5 <- lmer(dens2020 ~ SiAvg + dens2019 + month + (1|ISO_SOV1), data = wide)
p <- ggpredict(m5, terms = c("dens2019"))
plot(p)


m5 <- lmer(per ~ SiAvg + (1|ISO_SOV1), data = wide)
p <- ggpredict(m5, terms = c("SiAvg"))
plot(p)
