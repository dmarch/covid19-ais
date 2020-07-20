

library(ggplot2)
library(egg)
library(cowplot)
library(ggrepel)

# create output directory
out_dir <- "results/ais-global/eez/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


#----------------------------------------------------
# Part 1. Prepare data
#----------------------------------------------------

# import changes per EEZ
# see 07a_summary_eez.R
inputfile <- "data/out/ais-global/eez_change.csv"
df <- read.csv(inputfile)

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

# combine with changes
change <- merge(df, class, all.x=TRUE, by.x="ISO_SOV1", by.y="ISO")



#----------------------------------------------------
# Part 2. Boxplot per income class
#----------------------------------------------------

# prepare data
eez_data <- change
eez_data$var <- factor(eez_data$var, levels =c("COUNT", "CARGO", "TANKER", "PASSENGER", "FISHING", "OTHER"))

# boxplot
p <- ggplot(filter(eez_data, !is.na(income), var != "COUNT"), aes(factor(income), perlog)) +
  geom_hline(yintercept = 0, linetype="dotted") +
  stat_boxplot(geom="errorbar", width=0.2, lwd=0.3) +
  geom_boxplot(aes(fill = income), lwd=0.2) +
  scale_fill_brewer(palette="Blues", direction = -1) +
  #ylim(-100,100)+
  labs(y = "Relative change (%)", x = "Income class") +
  facet_wrap(var ~ .,  ncol = 1) +
  theme_bw(base_size = 12, base_family = "") +
  theme_article()+
  theme(legend.position =  "none")

p <- tag_facet(p, open = "", close = "")

# export multi-panel plot
out_file <- paste0(out_dir, "change_income_boxplot.png")
ggsave(out_file, p, width=10, height=20, units = "cm")




#----------------------------------------------------
# Part 3. Boxplot per region
#----------------------------------------------------

# boxplot
p <- ggplot(filter(eez_data, !is.na(income), var != "COUNT"), aes(factor(region), perlog)) +
  geom_hline(yintercept = 0, linetype="dotted") +
  stat_boxplot(geom="errorbar", width=0.2, lwd=0.3) +
  geom_boxplot(aes(fill = region), lwd=0.2) +
  scale_fill_brewer(palette="Set3") +
  #ylim(-100,100)+
  labs(y = "Relative change (%)", x = "Income class") +
  facet_wrap(var ~ .,  ncol = 1) +
  theme_bw(base_size = 12, base_family = "") +
  theme_article()+
  guides(fill=guide_legend(ncol=2)) +
  theme(legend.position="bottom",
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      legend.title=element_blank())

p <- tag_facet(p, open = "", close = "", tag_pool = letters[-c(1:5)])

# export multi-panel plot
out_file <- paste0(out_dir, "change_region_boxplot.png")
ggsave(out_file, p, width=10, height=20, units = "cm")




#----------------------------------------------------
# Part 4. Bubble plots
#----------------------------------------------------

p <- ggplot(filter(eez_data, var== "COUNT", !is.na(income)), aes(x=perlog, y=dens2020)) +
  geom_point(aes(colour = income, size=AREA_KM2), alpha = 0.6, stroke = 0, shape=16) +
  scale_color_brewer(palette="PuOr")+
  scale_size(range = c(2, 10)) +
  geom_vline(xintercept = 0, linetype="dotted", color="black") +
  scale_y_continuous(trans='log10') +
  ylab(expression(Traffic~density~(vessels/km^2))) + xlab("Change in traffic density (%)") +
  # geom_text_repel(data=filter(change, var== "COUNT", perlog > 25 | perlog < -50 | dens2020 > 1 | dens2020 < 0.005),
  #           aes(x=perlog, y=dens2020,label=TERRITORY1), size = 2, hjust=0,vjust=0) +
  labs(colour="Income", size=expression(Area~(km^2))) +
  facet_wrap(var ~ .,  ncol = 4) +
  theme_article() +
  background_grid(size.major = 0.3, color.major = "grey90") +
  guides(colour = guide_legend(override.aes = list(size=5)))

# export multi-panel plot
out_file <- paste0(out_dir, "change_all_bubble.png")
ggsave(out_file, p, width=20, height=10, units = "cm")



# same but per sector
p <- ggplot(filter(eez_data, var!= "COUNT", !is.na(income)), aes(x=perlog, y=dens2020)) +
  geom_point(aes(colour = income, size=AREA_KM2), alpha = 0.6, stroke = 0, shape=16) +
  scale_color_brewer(palette="PuOr")+
  scale_size(range = c(2, 10)) +
  geom_vline(xintercept = 0, linetype="dotted", color="black") +
  scale_y_continuous(trans='log10') +
  ylab(expression(Traffic~density~(vessels/km^2))) + xlab("Change in traffic density (%)") +
  labs(colour="Income", size=expression(Area~(km^2))) +
  facet_wrap(var ~ .,  ncol = 2) +
  theme_article() +
  background_grid(size.major = 0.3, color.major = "grey90") +
  guides(colour = guide_legend(override.aes = list(size=5)))

p <- tag_facet(p, open = "", close = "", tag_pool = letters[-1])

out_file <- paste0(out_dir,"change_sector_bubble.png")
ggsave(out_file, p, width=20, height=15, units = "cm")



