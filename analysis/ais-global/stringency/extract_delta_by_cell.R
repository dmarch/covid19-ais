#------------------------------------------------------------------------------
# extract_dens_by_EEZ
#------------------------------------------------------------------------------
# This script extract vessel density information per EEZ for each variable and timing
# The output is a data.frame where each row represent a grid cell assigned to a EEZ


## libraries
library(raster)
library(dplyr)
library(sf)
library(fasterize)
library(data.table)
library(foreach)
library(parallel)
library(doParallel)
source("scr/fun_common.R")
source("scr/fun_table.R")
source("scr/fun_ais.R")


## Prepare clusters
cores <- 1 # detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)


#----------------------------------------------------
# 1. Set input parameters
#----------------------------------------------------

# set input directory
input_dir <- "data/out/ais-global/delta/"

# output directory
out_dir <- "data/out/ais-global/eez/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# select variables to summarize
vars <- c("COUNT", "FISHING", "PASSENGER", "CARGO", "TANKER", "OTHER")

# set dates
dates_post <- c(
  seq.Date(as.Date("2020-01-01"), as.Date("2020-07-01"), by = "month")
) %>% format("%Y%m%d")

dates_pre <- c(
  seq.Date(as.Date("2019-01-01"), as.Date("2019-07-01"), by = "month")
) %>% format("%Y%m%d")




#----------------------------------------------------
# 4. Process data
#----------------------------------------------------

data_list <- list()
counter <- 1

for (i in 1:length(dates_post)){  ## process monthly data
  for (j in 1:length(vars)){  ## process variable 
    
    
    # set data to import
    idate_post <- dates_post[i]
    idate_pre <- dates_pre[i]
    jvar <- vars[j]
    print(paste("processing", idate_post, jvar))

    # import raster
    tif_file <- list.files(input_dir, full.names = TRUE, recursive=TRUE, pattern=sprintf("%s_%s_%s_delta.tif", jvar, idate_post, idate_pre))
    cnt <- raster(tif_file)
    
    # extract data
    df <- as.data.frame(cnt, xy=TRUE, na.rm=TRUE)
    names(df)[3] <-  dates_post[i]
    df <- mutate(df, var = jvar, month = month(dates_post[i]))
    
    ## Add cell ID
    df$idcell <- cellFromXY(cnt, cbind(df$x, df$y))
    
    # reorder
    df <- relocate(df, idcell, x, y, var, everything())
    data_list[[counter]] <- df
    counter <- counter + 1
  }
}



# combine lists
rbindlist()