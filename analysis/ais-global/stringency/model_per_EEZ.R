#-----------------------------------------------------------------------------------
# model_delta_EEZ.R
#-----------------------------------------------------------------------------------
# 1. Select EEZs for analysis
# 2. Prepare data combination for extraction
# 3. Extract delta values for each EEZ and month
# 4. Filter EEZ without good data
# 5. Generate plot of mean values with SD
# 6. Assess effect of Stringency Index



library(sf)
library(dplyr)
library(raster)
library(doParallel)
library(foreach)
library(parallel)
library(exactextractr)
library(lubridate)
library(ggplot2)
library(egg)
library(dplyr)
library(reshape2)
library(lmerTest)
library(lme4)
library(sjPlot)
library(sjmisc)
library(ggplot2)
require(lattice)  
library(ggeffects)
library(pdp)
library(xlsx)
library(merTools)
library(HLMdiag)
library(DHARMa)



## Set output plots
output_data <- "results/ais-global/eez"
if (!dir.exists(output_data)) dir.create(output_data, recursive = TRUE)


#--------------------------------------------
# 1. Select EEZs for analysis
#--------------------------------------------

# import EEZ
# source: marineregions
eez <- st_read("data/input/marine_regions/World_EEZ_v11_20191118_gpkg/eez_v11_covid.gpkg")
eez <- eez %>%
  mutate(TERRITORY1 = as.character(TERRITORY1),
         SOVEREIGN1 = as.character(SOVEREIGN1)) %>%
  filter(!POL_TYPE %in% c("Joint regime", "Overlapping claim"), # remove joint regimes
         TERRITORY1!="Antarctica",  # remove Antarctica
         TERRITORY1 == SOVEREIGN1,  # select main territories (excludes overseas)
         AREA_KM2>(769*3)) # remove EEZ smaller than 3 grid sizes at equator

# transform EEZ to mollweide
poly_moll <- st_transform(eez, crs = st_crs('+proj=moll'))


#--------------------------------------------------
# 2. Prepare data combination for extraction
#--------------------------------------------------


# select variables to summarize
vars <- c("COUNT", "FISHING", "PASSENGER", "CARGO", "TANKER", "OTHER")

# set dates
dates_post <- c(
  seq.Date(as.Date("2020-01-01"), as.Date("2020-06-01"), by = "month")
)

dates_pre <- c(
  seq.Date(as.Date("2019-01-01"), as.Date("2019-06-01"), by = "month")
)

# Prepare combinations
combinations <- data.frame(
  var = rep(vars, each=length(dates_pre)),
  date_post = rep(dates_post, length(vars)),
  date_pre =rep(dates_pre, length(vars)))

# set input directory
input_dir <- "data/out/ais-global/delta/"




#--------------------------------------------------
# 3. Extract delta values for each EEZ and month
#--------------------------------------------------

# Prepare cluster for parallel computing
cl <- makeCluster(12)
registerDoParallel(cl)

# Extract data for each combination  
data <- bind_rows(foreach(i = 1:nrow(combinations), .packages=c("dplyr", "raster", "exactextractr")) %dopar% { 
  
  # select file
  ivar <- combinations$var[i]
  idate_post <- combinations$date_post[i]
  idate_pre <- combinations$date_pre[i]
  tif_file <- list.files(input_dir, full.names = TRUE, recursive=TRUE, pattern=sprintf("%s_%s_%s_per_mol.tif", ivar, format(idate_post, "%Y%m%d"), format(idate_pre, "%Y%m%d")))
  
  # import raster
  rast <- raster(tif_file)
  
  # summarize data per polygon
  # weighted by the fraction of each cell that is covered by the polygon
  extr_mean <- exact_extract(rast, poly_moll, "mean")
  extr_med <- exact_extract(rast, poly_moll, "median")
  extr_sd <- exact_extract(rast, poly_moll, "stdev")
  extr_count <- exact_extract(rast, poly_moll, "count")
  extr_q25 <- exact_extract(rast, poly_moll, "quantile", quantiles = c(0.25))
  extr_q75 <- exact_extract(rast, poly_moll, "quantile", quantiles = c(0.75))
  
  # generate data.frame
  df <- data.frame(id = poly_moll$MRGID, var = ivar, date = idate_post, measure = "delta",
                   mean = extr_mean, median = extr_med, stdev = extr_sd, count = extr_count,
                   q25 = extr_q25, q75 = extr_q75)
  df
})

stopCluster(cl)  # Stop cluster


#--------------------------------------------------
# 4. Filter EEZ without good data
#--------------------------------------------------

## combine with EEZ
eez <- eez %>% st_drop_geometry() # remove geometry to convert to data.frame
data <- data %>%
  left_join(dplyr::select(eez, MRGID, TERRITORY1, SOVEREIGN1, ISO_SOV1, AREA_KM2), by = c("id" = "MRGID"))

# select combinations of EEZ and variable with at least 3 cells and data for all months (n=128)
select_countries <- data %>%
  group_by(id, var) %>%
  filter(count >= 3) %>%
  summarize(n_months = n()) %>%
  filter(n_months == 6)

length(unique(select_countries$id)) # 143

# filter data
data <- data %>%
  right_join(select_countries, by = c("id" = "id", "var" = "var"))


#--------------------------------------------------
# 5. Generate plot of mean values with SD
#--------------------------------------------------

# individual country profiles
sel_countries <- c("Spain", "France", "Italy", "China", "India", "United States", "South Korea", "Indonesia", "United Kingdom", "Belgium")


# calculate monthly summaries
month_sum <- data %>%
  group_by(date, var) %>%
  summarize(n = n(),
            avg = mean(mean, na.rm=TRUE),
            avg = mean(median, na.rm=TRUE))
            
# Plot time series all countries in the same plot
p1 <- ggplot(filter(data), aes(x = date)) +
  geom_line(aes(y = mean, group = TERRITORY1), size = 0.7, color="grey50", alpha=0.2) +
  #geom_line(data=filter(month_sum), aes(x=date, y = delta_avg), size = 1, color="black", alpha=0.8) +
  #geom_ribbon(data=month_sum, aes(x=date, ymin = delta_avg-delta_sd, ymax = delta_avg+delta_sd), fill="#3182bd", alpha=.2, linetype=0) +
  geom_hline(yintercept = 0, linetype="dotted") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b", expand = c(0, 0)) +
  ylab(expression(Absolute~change~(Delta~vessels~km^-2~month^-1))) +
  xlab("") +
  ylim(c(-100,100))+
  facet_wrap(var ~ ., ncol =2) +
  theme_article() +
  theme(legend.position = c(0.85, 0.2), legend.title = element_blank())


data$change_positive <- data$mean > 0
data$month <- month(data$date)

p <- ggplot(filter(data, var == "COUNT", TERRITORY1 %in% sel_countries), mapping=aes(x = month, y = mean, fill = change_positive)) +
  geom_errorbar(aes(ymin = mean-(mean < 0)*stdev, ymax = mean+(mean > 0)*stdev), width = 0.4, size=0.3) +
  #geom_col(alpha=1, width=0.8, size=2) +
  geom_bar(aes(colour=change_positive), stat = "identity", size=0.5, width=0.8) +
  #ylab(expression(Absolute~change~(Delta~vessel~transits~km^-2~month^-1))) +
  ylab(expression(atop(Mean~absolute~change, paste((Delta~vessel~transits~km^-2~month^-1))))) +
  xlab("Month") +
  scale_x_continuous(breaks = 1:6) +
  scale_fill_manual(values=c("#e34a33", "#9ecae1")) +
  scale_color_manual(values=c("#e34a33", "#9ecae1")) +
  facet_wrap(TERRITORY1 ~ ., ncol = 3) +
  theme_article() +
  theme(legend.position =  "none") +
  guides(fill = FALSE)


# mean and stdev
p3 <- ggplot(filter(data, var == "COUNT", TERRITORY1 %in% sel_countries), aes(x = date)) +
  geom_ribbon(aes(ymin = mean-stdev, ymax = mean+stdev, group = TERRITORY1, fill = TERRITORY1), alpha=.2, linetype=0) +
  geom_line(aes(y = mean, group = TERRITORY1, color = TERRITORY1), size = 1) +
  geom_hline(yintercept = 0, linetype="dotted") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b", expand = c(0, 0)) +
  #scale_colour_brewer(palette="Set2") +
  ylab(expression(Absolute~change~(Delta~vessels~km^-2~month^-1))) +
  xlab("") +
  facet_wrap(TERRITORY1 ~ ., ncol = 3) +
  theme_article() +
  theme(legend.position = "none")

# median and IQR
p3 <- ggplot(filter(data, var == "COUNT", TERRITORY1 %in% sel_countries), aes(x = date)) +
  geom_ribbon(aes(ymin = q25, ymax = q75, group = TERRITORY1, fill = TERRITORY1), alpha=.2, linetype=0) +
  geom_line(aes(y = median, group = TERRITORY1, color = TERRITORY1), size = 1) +
  geom_hline(yintercept = 0, linetype="dotted") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b", expand = c(0, 0)) +
  #scale_colour_brewer(palette="Set2") +
  ylab(expression(Absolute~change~(Delta~vessels~km^-2~month^-1))) +
  xlab("") +
  facet_wrap(TERRITORY1 ~ ., ncol = 3) +
  theme_article() +
  theme(legend.position = "none")



#-------------------------------------------------------
# 6. Assess effect of Stringency Index
#-------------------------------------------------------



# Import and prepare Stringency Idenx
#---------------------------------------------
# import Stringency index
# see stringency.R
si <- read.csv("data/out/stringency/stringency.csv")

# Calculate median
si_sum <- si %>%
  group_by(CountryCode, Month) %>%
  summarize(SiAvg = mean(StringencyIndex)) %>%
  mutate(date = as.Date(paste(2020, Month, 1, sep="-"))) %>%
  dplyr::select(-Month)

# Calculate class
si_sum$SiClass <- cut(si_sum$SiAvg, 4)
levels(si_sum$SiClass) <- c("low", "midlow", "midhigh", "high")




# Import World Bank indicators
#---------------------------------------------
# income groups and regions
# source: https://datahelpdesk.worldbank.org/knowledgebase/articles/906519-world-bank-country-and-lending-groups
classfile <- "data/input/worldbank/CLASS.xls"
class <- read.xlsx(classfile, sheetName = "List of economies", header=TRUE, startRow = 5, endRow = 224)
class <- dplyr::select(class, Code, Income.group, Region) %>% rename(ISO = Code, income = Income.group, region = Region)
class <- class[-1,]  # remove first row

# rename factors
class$income <- as.character(class$income)
class$income[class$income == "High income"] <- "High"
class$income[class$income == "Upper middle income"] <- "Upper middle"
class$income[class$income == "Lower middle income"] <- "Lower middle"
class$income[class$income == "Low income"] <- "Low"
class$income <- factor(class$income, levels=c("High", "Upper middle", "Lower middle", "Low"))

# Combinedata and Stringency Idenx
#---------------------------------------------
sdata <- data %>%
  filter(ISO_SOV1 %in% unique(si_sum$CountryCode)) %>%
  left_join(si_sum, by = c("ISO_SOV1" = "CountryCode", "date" = "date")) %>%
  left_join(class, by = c("ISO_SOV1" = "ISO"))# %>%
  #mutate(SiAvg = replace(SiAvg, is.na(SiAvg), 0)) %>%  ## CARE: FIRST REPLACE THOSE NA IN 2020. NOT ALL COUNTRIES HAVE DATA.
  #mutate(covid = ifelse(SiAvg > 25, "after", "before"))


# select combinations of EEZ with data for all months (n=128)
select_countries <- sdata %>%
  group_by(id, var) %>%
  filter(!is.na(SiAvg)) %>%
  summarize(n_months = n()) %>%
  filter(n_months == 6)

length(unique(select_countries$id)) # 123

# filter data
sdata <- sdata %>%
  right_join(select_countries, by = c("id" = "id", "var" = "var"))




#---------------------------------------------
# Model delta
#---------------------------------------------
# change from baseline as a dependent variable in mixed models, specifically mixed models for repeat measurements (MMRM).


model_list <- list()

for(j in 1:length(vars)){
  
  # filter data by category
  jvar <- vars[j]
  jdata <- filter(sdata, var == jvar)
  
  # model
  m1 <- lmer(median ~ SiAvg*income + (1|region/ISO_SOV1), data = jdata)

  # append model to list
  model_list[[j]] <- m1
}




