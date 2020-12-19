#-----------------------------------------------------------------------------
# model_dens_EEZ       Model density data per EEZ
#-----------------------------------------------------------------------------
# For each month amd country, add Stringency Index

# https://www.nature.com/articles/s41598-020-75848-2
# then calculated the country level monthly averaged NO2, AOD, and XCO2 values: https://www.nature.com/articles/s41467-020-18922-7#Sec11
# On average, pollution concentrations decreased by 23.0% for NO2, 15.4% for PM2.5, and 12.5% for CO during January–March 2020 relative to the same period in 2019.

# random effect for percentange change in mobility: https://www.nature.com/articles/s41598-020-76763-2

# baseline for air quality. https://advances.sciencemag.org/content/6/28/eabc2992

library(dplyr)
library(lubridate)
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
library(sf)
#library(mixedup) # remotes::install_github('m-clark/mixedup')


## Set output plots
output_data <- "results/ais-global/eez"
if (!dir.exists(output_data)) dir.create(output_data, recursive = TRUE)



#---------------------------------------------
# Import and prepare density data 
#---------------------------------------------
# import data
# density data at grid cell level
data <- read.csv("data/out/ais-global/eez/eez_ais_density.csv")


# import EEZ
# source: marineregions
eez <- st_read("data/input/marine_regions/World_EEZ_v11_20191118_gpkg/eez_v11_covid.gpkg")
eez <- eez %>%
  mutate(TERRITORY1 = as.character(TERRITORY1),
         SOVEREIGN1 = as.character(SOVEREIGN1)) %>%
  filter(!POL_TYPE %in% c("Joint regime"), # remove joint regimes
         TERRITORY1!="Antarctica",  # remove Antarctica
         TERRITORY1 == SOVEREIGN1,  # select main territories (excludes overseas)
         AREA_KM2>(769*3)) # remove EEZ smaller than 3 grid sizes at equator

# select EEZ:
# main territories (TERRITORY1 name = SOVEREIGN1 name)
data <- data %>%
  mutate(date = parse_date_time(date, "Ymd")) %>%
  filter(MRGID %in% eez$MRGID)


# filter by ship type
vars <- c("COUNT",  "CARGO", "TANKER","PASSENGER", "FISHING","OTHER")

# create list to append all models
model_list <- list()
change_list <- list()

# import country map
library(tmap)
data("World")


## LOOP FOR EACH SHIP CATEGORY
for(j in 1:length(vars)){
  
# select data for category (j)
jvar <- vars[j]
jdata <- filter(data, var == jvar)

# Summary statistics
sdata <- jdata %>%
  group_by(ISO_SOV1, TERRITORY1, date) %>%
  filter(vessels_km2 > 0) %>%
  summarize(medDens = median(vessels_km2),
            avgDens = mean(vessels_km2),
            surface = sum(eez_km2),
            ncells = n()) %>%
  mutate(month = month(date),
         year = year(date),
         avgDensLog = log10(avgDens))

# select countries with at least 3 cells and data for all months (n=128)
select_countries <- sdata %>%
  group_by(ISO_SOV1, TERRITORY1) %>%
  filter(ncells > 3) %>%
  summarize(n_months = n()) %>%
  filter(n_months == 14)

# filter dataset per selected countries
jdata <- jdata %>% filter(ISO_SOV1 %in% select_countries$ISO_SOV1)
sdata <- sdata %>% filter(ISO_SOV1 %in% select_countries$ISO_SOV1)


# calculate change
change <- sdata %>%
  dcast(ISO_SOV1 + TERRITORY1 + month ~ year, value.var=c("avgDens")) %>%
  rename(dens2019 = `2019`, dens2020 = `2020`) %>%
  mutate(delta = dens2020 - dens2019,
         per = 100*(delta/dens2019),
         date = as.Date(paste(2020, month, 1, sep="-")),
         var = jvar)


#---------------------------------------------
# Plots
#---------------------------------------------

library(ggplot2)
library(egg)


# calculate monthly summaries
month_sum <- change %>%
  group_by(date) %>%
  summarize(n = n(),
            delta_avg = mean(delta, na.rm=TRUE),
            delta_sd = sd(delta, na.rm=TRUE),
            delta_sem = delta_sd/sqrt(n),
            delta_lci = delta_avg + qt((1-0.95)/2, df=n-1) * delta_sem,
            delta_uci = delta_avg - qt((1-0.95)/2, df=n-1) * delta_sem,
            per_avg = mean(per, na.rm=TRUE),
            per_sd = sd(per, na.rm=TRUE))

# country profiles (relative change per month)
# Plot time series all countries in the same plot
p1 <- ggplot(change, aes(x = date)) +
  geom_line(aes(y = delta, group = TERRITORY1), size = 0.7, color="grey50", alpha=0.2) +
  geom_line(data=month_sum, aes(x=date, y = delta_avg), size = 1, color="black", alpha=0.8) +
  geom_ribbon(data=month_sum, aes(x=date, ymin = delta_avg-delta_sd, ymax = delta_avg+delta_sd), fill="#3182bd", alpha=.2, linetype=0) +
  geom_hline(yintercept = 0, linetype="dotted") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b", expand = c(0, 0)) +
  ylab(expression(Absolute~change~(Delta~vessels~km^-2~month^-1))) +
  xlab("") +
  theme_article() +
  theme(legend.position = c(0.85, 0.2), legend.title = element_blank())

# export multi-panel plot
out_file <- paste0(output_data, "/", jvar, "_eez_profiles_delta.png")
ggsave(out_file, p1, width=12, height=10, units = "cm")

p2 <- ggplot(change, aes(x = date)) +
  geom_line(aes(y = per, group = TERRITORY1), size = 0.7, color="grey50", alpha=0.2) +
  geom_line(data=month_sum, aes(x=date, y = per_avg), size = 1, color="black", alpha=0.8) +
  geom_ribbon(data=month_sum, aes(x=date, ymin = per_avg-per_sd, ymax = per_avg+per_sd), fill="#3182bd", alpha=.2, linetype=0) +
  geom_hline(yintercept = 0, linetype="dotted") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b", expand = c(0, 0)) +
  ylab("Relative change (%)") +
  xlab("") +
  theme_article() +
  theme(legend.position = c(0.85, 0.2), legend.title = element_blank())

# export multi-panel plot
out_file <- paste0(output_data, "/", jvar, "_eez_profiles_per.png")
ggsave(out_file, p2, width=12, height=10, units = "cm")


# individual country profiles
sel_countries <- c("Spain", "France", "Italy", "China", "India", "United States", "South Korea", "Indonesia", "United Kingdom")


p3 <- ggplot(filter(change, TERRITORY1 %in% sel_countries), aes(x = date)) +
  geom_line(aes(y = delta, group = TERRITORY1, color = TERRITORY1), size = 1) +
  geom_hline(yintercept = 0, linetype="dotted") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b", expand = c(0, 0)) +
  #scale_colour_brewer(palette="Set2") +
  ylab(expression(Absolute~change~(Delta~vessels~km^-2~month^-1))) +
  xlab("") +
  facet_wrap(TERRITORY1 ~ ., ncol = 3) +
  theme_article() +
  theme(legend.position = "none")

# export multi-panel plot
out_file <- paste0(output_data, "/", jvar, "_eez_individual_profiles_delta.png")
ggsave(out_file, p3, width=18, height=14, units = "cm")


p4 <- ggplot(filter(change, TERRITORY1 %in% sel_countries), aes(x = date)) +
  geom_line(aes(y = per, group = TERRITORY1, color = TERRITORY1), size = 1) +
  geom_hline(yintercept = 0, linetype="dotted") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b", expand = c(0, 0)) +
  #scale_colour_brewer(palette="Set2") +
  ylab("Relative change (%)") +
  xlab("") +
  facet_wrap(TERRITORY1 ~ ., ncol = 3) +
  theme_article() +
  theme(legend.position = "none")

out_file <- paste0(output_data, "/", jvar, "_eez_individual_profiles_per.png")
ggsave(out_file, p4, width=18, height=14, units = "cm")




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
  mutate(date = as.Date(paste(2020, Month, 1, sep="-"))) %>%
  dplyr::select(-Month)

# Calculate class
si_sum$SiClass <- cut(si_sum$SiAvg, 4)
levels(si_sum$SiClass) <- c("low", "midlow", "midhigh", "high")



#---------------------------------------------
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


#---------------------------------------------
# Combine Count data and Stringency Idenx
#---------------------------------------------

# combine with countries with EEZ by ISO code
# first, filter countries without SI
# then combine data by date and country code
sdata <- sdata %>%
  filter(ISO_SOV1 %in% unique(si_sum$CountryCode)) %>%
  left_join(si_sum, by = c("ISO_SOV1" = "CountryCode", "date" = "date")) %>%
  left_join(class, by = c("ISO_SOV1" = "ISO")) %>%
  mutate(SiAvg = replace(SiAvg, is.na(SiAvg), 0)) %>%  ## CARE: FIRST REPLACE THOSE NA IN 2020. NOT ALL COUNTRIES HAVE DATA.
  mutate(covid = ifelse(SiAvg > 25, "after", "before"))

# add data into change
change <- change %>%
  filter(ISO_SOV1 %in% unique(si_sum$CountryCode)) %>%
  left_join(si_sum, by = c("ISO_SOV1" = "CountryCode", "date" = "date")) %>%
  left_join(class, by = c("ISO_SOV1" = "ISO")) %>%
  filter(!is.na(SiAvg))  ## CARE: FIRST REPLACE THOSE NA IN 2020. NOT ALL COUNTRIES HAVE DATA.

# add data into cell data
jdata <- jdata %>%
  filter(ISO_SOV1 %in% unique(si_sum$CountryCode)) %>%
  left_join(si_sum, by = c("ISO_SOV1" = "CountryCode", "date" = "date")) %>%
  left_join(class, by = c("ISO_SOV1" = "ISO")) %>%
  mutate(SiAvg = replace(SiAvg, is.na(SiAvg), 0)) %>%  ## CARE: FIRST REPLACE THOSE NA IN 2020. NOT ALL COUNTRIES HAVE DATA.
  mutate(covid = ifelse(SiAvg > 25, "after", "before"))


# append change df to list
change_list[[j]] <- change


#---------------------------------------------
# Model delta
#---------------------------------------------
# change from baseline as a dependent variable in mixed models, specifically mixed models for repeat measurements (MMRM).

# model
m1 <- lmer(per ~ SiAvg*income + (1|region/ISO_SOV1), data = change)

# append model to list
model_list[[j]] <- m1


#---------------------------------------------
# Residuals (plot map and moran)
#---------------------------------------------

# # get residuals
# change$resid <- residuals(m1)
# 
# # summarize residuals per country
# resid_eez <- change %>%
#         group_by(ISO_SOV1) %>%
#         summarize(resid_avg = mean(resid),
#                   resid_sd = sd(resid))
#         
# # plot average
# 
# 
# # combine with countries with EEZ by ISO code
# map_data <- left_join(World, resid_eez, by = c("iso_a3" = "ISO_SOV1"))
# map_data <- left_join(World, change, by = c("iso_a3" = "ISO_SOV1"))
# map_data <- map_data %>% filter(name != "Antarctica")
# 
# # plot data (average)
# library(pals)
# p1 <- tm_shape(map_data) +
#   tm_polygons("delta",
#               #style="quantile",
#               title=("Residuals (average)"),
#               #palette = "YlGnBu", #brewer.ylgnbu(10) # YlGnBu, -RdYlBu
#               border.col = "grey60", 
#               border.alpha = 0.3) +
#   tm_layout(legend.title.size=1) +
#   tm_facets(by="month")
# 
# p2 <- tm_shape(map_data) +
#   tm_polygons("resid_sd",
#               title=("Residuals (sd)"),
#               palette = brewer.ylgnbu(10), #brewer.ylgnbu(10) # YlGnBu, -RdYlBu
#               border.col = "grey60", 
#               border.alpha = 0.3) +
#   tm_layout(legend.title.size=1)

}


#----------------------------------------------------------------------
# Summary table for LMM
#----------------------------------------------------------------------
out_file <- paste0(output_data, "/", jvar, "_model.doc")
out_file <- paste0(output_data, "/model.doc")

tab_model(m1, p.val = "kr", show.df = FALSE, show.icc = FALSE, show.re.var = FALSE, file=out_file)
tab_model(model_list[[1]], p.val = "kr", show.df = FALSE, show.icc = FALSE, show.re.var = FALSE, file=out_file)


tab_model(model_list,
          p.val = "kr",
          collapse.ci = TRUE,
          show.df = FALSE,
          show.icc = FALSE,
          show.re.var = TRUE,
          p.style = "stars",
          string.ci = "Conf. Int (95%)",
          string.est = "Estimates (Conf. Int 95%)",
          dv.labels = c("All vessels", "Cargo", "Tanker", "Passenger", "Fishing", "Other"),
          file=out_file)



#----------------------------------------------------------------------
# Plots for LMM
#----------------------------------------------------------------------

for(j in 1:length(vars)){
  
  # extract model
  m <- model_list[[j]]

  # create interaction plot
  p <- ggpredict(m, terms = c("SiAvg", "income"))
  
  # output file
  jvar <- vars[j]
  p <- plot(p, facet = TRUE) +
    labs(x="Stringency Index",
         y="Relative change (%)",
         title="") +
    theme_article()
  
  #
  out_file <- paste0(output_data, "/", jvar, "_int.png")
  ggsave(out_file, p, width=10, height=10, units = "cm")
}




#----------------------------------------------------------------------
# K-means (per country)
#----------------------------------------------------------------------

library(factoextra)
library(NbClust)

# combine all change dataframe
change_all <- rbindlist(change_list)

change_cast <- change_all %>%
  dcast(ISO_SOV1 + var ~ month, value.var=c("delta"))

# filter rows with NA
change_cast <- change_cast %>% filter(!is.na(`7`))

# convert to matrix
change_mat <- as.matrix(change_cast[, 3:9])  


#### Select number of clusters -----------

# Elbow method
fviz_nbclust(change_mat, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method")

# Silhouette method
fviz_nbclust(change_mat, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")

# Gap statistic
# nboot = 50 to keep the function speedy. 
# recommended value: nboot= 500 for your analysis.
# Use verbose = FALSE to hide computing progression.
set.seed(123)
fviz_nbclust(change_mat, kmeans, nstart = 25,  method = "gap_stat", nboot = 50)+
  labs(subtitle = "Gap statistic method")

NbClust(data = change_mat, diss = NULL, distance = "euclidean",
        min.nc = 2, max.nc = 15, method = "kmeans")

# Initialise some variables
kRange <- seq(from=2, to=10, by=1)
intra <- rep(NA, length(kRange)) # intracluster sum-of-squares
inter <- rep(NA, length(kRange)) # intercluster sum-of-squares

# Loop across desired range of ks
for (k in kRange) {
  mdl <- kmeans(x=change_mat, centers=k)
  intra[k-1] <- mdl$tot.withinss # it’s (k-1) because k starts from 2
  inter[k-1] <- mdl$betweenss 
}

# Plot inter/intercluster sum-of-squares as a function of $k$
yMin <- min(c(intra , inter)) 
yMax <- max(c(intra , inter))
plot(kRange, inter, type="o", pch=1, col="blue", lwd=3, lty=1, 
     ylim=c(yMin, yMax), xlab="k", ylab='Sum-of-squares')
points(kRange, intra, type="o", pch=1, col="red", lwd=3)
legend("right", c("inter-cluster", "intra-cluster"), 
       bty="n", col=c("blue", "red"), pch=1, lwd=3, lty=1)

#### Select number of clusters -----------




# k-means
res_km <- kmeans(change_mat, 5, nstart = 10)

# assign class
change_cast$class <- as.character(res_km$cluster)

# combine into original format
change_all <- change_all %>%
  left_join(dplyr::select(change_cast, ISO_SOV1, var, class), by=c("ISO_SOV1", "var")) %>%
  filter(!is.na(class))

# prepare centroids
centers <- melt(data.table(ID = 1:nrow(res_km$centers),
                           class = as.character(1:nrow(res_km$centers)),
                           res_km$centers),
                id.vars = c("ID", "class"),
                variable.name = "Time",
                variable.factor = FALSE
)
centers[, Time := as.integer(gsub("V", "", Time))]



# plot the results
#http://ofdataandscience.blogspot.com/2013/03/capital-bikeshare-time-series-clustering.html
change_all$grp <- paste(change_all$ISO_SOV1, change_all$var, sep="_")

ggplot(change_all) +
  geom_line(aes(month, delta, group = grp, colour=class),  size=0.7, alpha = 0.5) +
  geom_line(data = centers, aes(Time, value), alpha = 0.80, colour="black", size = 1) +
  scale_color_brewer(palette="Set2") +
  geom_hline(yintercept = 0, linetype="dotted") +
  labs(x = "Month", y = "Relative change (%)") +
  #ylim(c(-50,50))+
  #facet_wrap(var ~ class, ncol = 6, scales="free")+
  facet_wrap(. ~ class, ncol = 2, scales="free")+
  theme_article()


#----------------------------------------------------------------------
# dtwclust (per cell)
#----------------------------------------------------------------------

library(dtwclust)

hc_sbd <- tsclust(change_mat, type = "h", k = 20L,
                  preproc = zscore, seed = 899,
                  distance = "sbd", centroid = shape_extraction,
                  control = hierarchical_control(method = "average"))
plot(hc_sbd)


plot(hc_sbd, type = "sc")


#----------------------------------------------------------------------
# K-means (per cell)
#----------------------------------------------------------------------

data2 <- data

data2$month <- month(data2$date)
data2$year <- year(data2$date)

# calculate change
change <- data2 %>%
  dcast(idcell + var + month ~ year, value.var=c("vessels_km2")) %>%
  rename(dens2019 = `2019`, dens2020 = `2020`) %>%
  mutate(delta = dens2020 - dens2019,
         per = 100*(delta/dens2019),
         date = as.Date(paste(2020, month, 1, sep="-")))


change_cast <- change %>%
  dcast(idcell + var ~ month, value.var=c("per"))

# filter rows with NA
change_cast <- change_cast %>% filter(!is.na(`7`))


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
# Model delta
#---------------------------------------------
# change from baseline as a dependent variable in mixed models, specifically mixed models for repeat measurements (MMRM).

#change$perT <- asin(sqrt(change$per))

m1 <- lmer(per ~ SiAvg*income + (1|region/ISO_SOV1), data = change)
p <- ggpredict(m1, terms = c("SiAvg"))
p <- ggpredict(m1, terms = c("SiAvg", "month"))
p <- ggpredict(m1, terms = c("SiAvg", "income"))

p <- ggpredict(m1, terms = c("income"))
p <- ggpredict(m1, terms = c("ISO_SOV1"))
plot(p)

# table
out_file <- paste0(output_data, "/", jvar, "_model.doc")
tab_model(m1, p.val = "kr", show.df = FALSE, show.icc = FALSE, show.re.var = FALSE, file=out_file)
tab_model(model_list)


# Check residuals
plot(m1)
Plot.Model.F.Linearity<-plot(resid(m1), change$per) 


## Caterpillar plot
lattice::dotplot(ranef(m1, condVar = TRUE))


# confidence intervals
confint(m1)

rr1 <- ranef(m1, condVar = TRUE)

dotplot(rr1$`ISO_SOV1:region`)

dd <- as.data.frame(rr1$ISO_SOV1)


dd <- as.data.frame(rr1)
ggplot(dd, aes(y=grp,x=condval)) +
    geom_point() + facet_wrap(~term,scales="free_x") +
    geom_errorbarh(aes(xmin=condval -2*condsd,
                       xmax=condval +2*condsd), height=0)


# interval estimates
predictInterval(m1)   # for various model predictions, possibly with new data
REsim(m1)             # mean, median and sd of the random effect estimates
plotREsim(REsim(m1))  # plot the interval estimates

# predict
predict_with_re <- predict(m1)


re = ranef(m1)
qplot(x = re$region, geom = 'density')
























#---------------------------------------------
# Model change from baseline
#---------------------------------------------
# https://core.ac.uk/download/pdf/3165006.pdf

change$dens2020log <- log(change$dens2020)
change$dens2019log <- log(change$dens2019)

m1 <- lmer(dens2020log ~ dens2019log + SiAvg + income + (SiAvg|ISO_SOV1), data = change)
m1 <- lmer(dens2020log ~ dens2019log + SiAvg*income + (1|region/ISO_SOV1), data = change)


p <- ggpredict(m1, terms = c("SiAvg"))
p <- ggpredict(m1, terms = c("ISO_SOV1"))
p <- ggpredict(m1, terms = c("SiAvg", "income"))
p <- ggpredict(m1, terms = c("SiAvg", "income"))

plot(p)


p <- ggpredict(m1, terms = c("SiAvg", "income"), type="re")
plot(p)

# Check residuals
plot(m1)

# check variance
Model.F.Res<- residuals(m1) #extracts the residuals and places them in a new column in our original data table
Abs.Model.F.Res <-abs(Model.F.Res) #creates a new column with the absolute value of the residuals
Model.F.Res2 <- Abs.Model.F.Res^2 #squares the absolute values of the residuals to provide the more robust estimate
change$Model.F.Res2 <- Model.F.Res2
Levene.Model.F <- lm(Model.F.Res2 ~ ISO_SOV1, data=change) #ANOVA of the squared residuals
anova(Levene.Model.F) #displays the results
# no Homogeneity of Variance


#---------------------------------------------
# Model density
#---------------------------------------------
#sdata$month <- as.factor(sdata$month)
m1 <- lmer(avgDensLog ~ covid*month + (1|region/ISO_SOV1), data = sdata)
p <- ggpredict(m1, terms = c("SiAvg"))
p <- ggpredict(m1, terms = c("covid"))
p <- ggpredict(m1, terms = c("month"))
p <- ggpredict(m1, terms = c("SiAvg", "month"))
p <- ggpredict(m1, terms = c("covid", "month"))
p <- ggpredict(m1, terms = c("income"))
plot(p)


m1 <- lmer(avgDensLog ~ covid*month + income + (1|region/ISO_SOV1), data = sdata)

plot(m1)
plotREsim(REsim(m1))  # plot the interval estimates

# Check residuals
# https://ademos.people.uic.edu/Chapter18.html
Plot.Model.F.Linearity<-plot(resid(m1), jdata$avgDensLog) 

# check variance
jdata$Model.F.Res<- residuals(m1) #extracts the residuals and places them in a new column in our original data table
jdata$Abs.Model.F.Res <-abs(jdata$Model.F.Res) #creates a new column with the absolute value of the residuals
jdata$Model.F.Res2 <- jdata$Abs.Model.F.Res^2 #squares the absolute values of the residuals to provide the more robust estimate
Levene.Model.F <- lm(Model.F.Res2 ~ ISO_SOV1, data=jdata) #ANOVA of the squared residuals
anova(Levene.Model.F) #displays the results
# no Homogeneity of Variance


Plot.m1 <- plot(m1) #creates a fitted vs residual plot
Plot.m1


