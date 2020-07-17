#---------------------------------------------------------------------------
# 06_summary_global     Summarize global changes
#---------------------------------------------------------------------------
# Estimate the proportion (%) of the ocean that decreased and increased


## libraries
library(doParallel)
library(foreach)
library(parallel)
library(exactextractr)
library(sf)
library(raster)
library(tmap)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(xlsx)
library(egg)
library(cowplot)
library(wbstats)
library(flextable)
library(officer)
library(webshot)
source("scr/fun_common.R")
source("https://raw.githubusercontent.com/dmarch/abigoos/master/R/utils.R")


# set input dir
input_delta<- "data/out/ais-global/delta/"
input_dens<- "data/out/ais-global/density/"


# create output directory
out_dir <- "results/ais-global/summary"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# select variables to summarize
vars <- c("CARGO", "TANKER", "PASSENGER", "FISHING", "OTHER", "COUNT")
vars_names <- c("Cargo", "Tanker", "Passenger", "Fishing", "Other", "All vessels")

# import ocean mask
ocean <- raster("data/out/ais-global/oceanmask_mol.nc")
total_area <- binsurf(ocean)

# define study area code
code <- "global"


#----------------------------------------------------
# Extract data from shipping rasters
#----------------------------------------------------
# We sumarize data using average.

# Prepare cluster for parallel computing
cl <- makeCluster(6)
registerDoParallel(cl)

# Extract data for each combination  
data <- bind_rows(foreach(i = 1:length(vars), .packages=c("dplyr", "raster", "exactextractr")) %dopar% { 
  
  # select file
  # mask by study area (only applies for MPA)
  ivar <- vars[i]
  r1 <- raster(paste0(input_dens, sprintf("20190401_%s_dens_moll.tif", ivar))) %>% mask(ocean)
  r2 <- raster(paste0(input_dens, sprintf("20200401_%s_dens_moll.tif", ivar))) %>% mask(ocean)
  delta <- raster(paste0(input_delta, sprintf("%s_delta_mol.tif", ivar))) %>% mask(ocean)
  
  # occupancy (any vessel)
  r1q <- quantile(r1, prob = 0)
  r1ext <- reclassify(r1, c(-Inf, r1q, NA, r1q, Inf, 1))  # create a ocean mask
  r2q <- quantile(r2, prob = 0)
  r2ext <- reclassify(r2, c(-Inf, r2q, NA, r2q, Inf, 1))  # create a ocean mask
  area2 <- binsurf(r2ext)
  occ_all_2019 <- 100 * (binsurf(r1ext)/total_area)
  occ_all_2020 <- 100 * (binsurf(r2ext)/total_area)
  occ_all_change <- occ_all_2020 - occ_all_2019
  
  
  # area increase/decrease
  diff_area <- sum((r1ext * 2), r2ext, na.rm=TRUE)  
  scell <- (res(diff_area)[1] * res(diff_area)[2]) / 1000000    # Calculate surface of one cell (in km2)
  df_area <- freq(diff_area) %>% data.frame()  # 0, no traffic; 1 (increase in 2020), 2(decrease in 2020), 3(no change)
  area_increase <- (df_area$count[df_area$value == 1] * scell)
  area_decrease <- (df_area$count[df_area$value == 2] * scell)
  prop_area_increase <- 100 * (area_increase / binsurf(r1ext))
  prop_area_decrease <- 100 * (area_decrease / binsurf(r1ext))
  
  
  # occupancy (high traffic)
  r1q <- quantile(r1, prob = 0.8)
  r1ext <- reclassify(r1, c(-Inf, r1q, NA, r1q, Inf, 1))  # create a ocean mask
  r2q <- quantile(r2, prob = 0.8)
  r2ext <- reclassify(r2, c(-Inf, r2q, NA, r2q, Inf, 1))  # create a ocean mask
  occ_high_2019 <- 100 * (binsurf(r1ext)/total_area)
  occ_high_2020 <- 100 * (binsurf(r2ext)/total_area)
  occ_high_change <- occ_high_2020 - occ_high_2019
  high_dens_2019 <- r1q
  high_dens_2020 <- r2q
  
  # calculate surface of cells with negative values
  rneg <- delta
  rneg[rneg >= 0] <- NA
  rneg[rneg < 0] <- 1
  neg_area <- binsurf(rneg)
  prop_decrease <-  100*(neg_area/total_area)
  
  # calculate surface of cells with positive values
  rpos <- delta
  rpos[rpos <= 0] <- NA
  rpos[rpos > 0] <- 1
  pos_area <- binsurf(rpos)
  prop_increase <-  100*(pos_area/total_area)
  
  # calculate % of the ocean with no change
  prop_nochange <- 100 - sum(prop_decrease, prop_increase)
  
  # global average change
  avg_change <- values(delta) %>% mean(na.rm=TRUE)
  
  
  # generate data.frame
  df <- data.frame(var = ivar,
                   prop_decrease = round(prop_decrease, 2),
                   prop_increase = round(prop_increase, 2),
                   prop_nochange = round(prop_nochange, 2),
                   avg_change = round(avg_change, 5),
                   area2 = round(area2/1e6, 1),
                   area_increase = round(area_increase/1e6, 1),
                   prop_area_increase = round(prop_area_increase, 2),
                   area_decrease = round(area_decrease/1e6, 1),
                   prop_area_decrease = round(prop_area_decrease, 2),
                   occ_all_2019, occ_all_2020, occ_all_change,
                   occ_high_2019, occ_high_2020, occ_high_change,
                   high_dens_2019, high_dens_2020)
  
  df
})

stopCluster(cl)  # Stop cluster

# save
write.csv(data, paste0(out_dir, sprintf("/%s_summary.csv", code)), row.names=FALSE)


#----------------------------------------------------
# Create table
#----------------------------------------------------

data <- read.csv(paste0(out_dir, sprintf("/%s_summary.csv", code)))

# Prepare table
data$var <- vars_names  # change names
data <- data[,c(1:10)]  # select columns
data$avg_change <- data$avg_change * 1000

# create flextable
ft <- flextable(data)  # data is a data.frame

# set decimal formats
ft <- colformat_num(ft, j = c(2,3,4,8,10), digits = 1)
ft <- colformat_num(ft, j = 5, digits = 2)

# set labels for headers (data.frame columns)
# now only those labels without complex text. To break lines, use "\n"
ft <- set_header_labels(ft, var = "Vessel \ncategory",
                        prop_decrease = "Decrease \n(%)",
                        prop_increase = "Increase \n(%)",
                        prop_nochange = "No change \n(%)",
                        prop_area_increase = "New area \n(%)",
                        prop_area_decrease = "Reduced area \n(%)")

# set labels for complex text
# include superscripts
ft <- compose(ft, part = "header", j = "area2",
              value = as_paragraph("Area \n(x10", as_sup("6"), " km",as_sup("2"),")"))
ft <- compose(ft, part = "header", j = "area_increase",
              value = as_paragraph("New area \n(x10", as_sup("6"), " km",as_sup("2"),")"))
ft <- compose(ft, part = "header", j = "area_decrease",
              value = as_paragraph("Reduced area \n(x10", as_sup("6"), " km",as_sup("2"),")"))
ft <- compose(ft, part = "header", j = "avg_change",
              value = as_paragraph("Average change \n(x10", as_sup("-3"), " vessels/", "km",as_sup("2"),")"))

# add rows on top to include groups
# also include vertical lines to separate groups
ft <- add_header_row(ft, values = c("", "Traffic density", "Occupancy"), colwidths = c(1,4,5))
ft <- theme_booktabs(ft)
ft <- vline(ft, i = NULL, j = c(1,5), border = fp_border(color = "black", style = "solid", width = 1), part = "all")
ft <- align(ft, i = NULL, j = c(2:10), align = "center", part = "all")  # center columns except vessel category

# general configuration of table
ft <- bold(ft, part = "header") # bold header
ft <- fontsize(ft, part = "all", size = 11)  # font size
ft <- set_table_properties(ft, width = 0.6, layout = "autofit")  # autofit to width

# export table for Word
save_as_docx(ft, path = paste0(out_dir, sprintf("/%s_summary.docx", code)))

# export as image
img_file <- paste0(out_dir, sprintf("/%s_summary.png", code))
save_as_image(ft, path = img_file)

