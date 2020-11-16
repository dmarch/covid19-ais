# summarize deltas
# calculate average and coefficient of variation of deltas during 2020


library(stringr)
library(raster)
library(lubridate)
source("scr/fun_common.R")  # bb()
source("scr/fun_ais.R")




# set input folder
input_dir <- "data/out/ais-global/delta/"

# create output directory
out_dir <- "data/out/ais-global/delta/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)



# locate all delta files

input_dir <- "data/out/ais-global/delta/"
tif_file <- list.files(input_dir, full.names = TRUE, recursive=TRUE, pattern="delta_mol.tif$")

count_files <- which(str_detect(tif_file, "COUNT")==TRUE)
tif_file <- tif_file[count_files]


jvar <- "COUNT"

dates_post <- c(
  seq.Date(as.Date("2020-01-01"), as.Date("2020-06-01"), by = "month")
) %>% format("%Y%m%d")

dates_pre <- c(
  seq.Date(as.Date("2019-01-01"), as.Date("2019-06-01"), by = "month")
) %>% format("%Y%m%d")



s <- stack()

for(i in 1:length(dates_post)){
  
  idate_post <- dates_post[i]
  idate_pre <- dates_pre[i]
  
  tif_file <- list.files(input_dir, full.names = TRUE, recursive=TRUE, pattern=sprintf("%s_%s_%s_delta_mol.tif", jvar, idate_post, idate_pre))
  r <- raster(tif_file)
  s <- stack(s, r)
  
}



u <- mean(s, na.rm=TRUE)
s <- calc(s, fun=sd, na.rm=TRUE)
cv <- s/u

i=1
median(values(s[[i]]), na.rm=TRUE)
hist(values(s[[i]]), breaks=500, xlim = c(-0.05, 0.05))
abline(v = median(values(s[[i]]), na.rm=TRUE), col="red", lwd=1, lty=1)

