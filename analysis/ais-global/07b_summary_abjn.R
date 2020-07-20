#-----------------------------------------------------------------------------------------
# summary_abjn          Summarize data by Areas Beyond National Jurisdiction (ABJN)
#-----------------------------------------------------------------------------------------


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
library(reshape2)
library(flextable)
library(officer)
source("scr/fun_common.R")
source("scr/fun_table.R")
source("https://raw.githubusercontent.com/dmarch/abigoos/master/R/utils.R")  # bb()


# set input dir
input_dir <- "data/out/ais-global/density/"

# create output directory
out_dir <- "results/ais-global/abjn/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# select variables to summarize
vars <- c("COUNT", "FISHING", "PASSENGER", "CARGO", "TANKER", "OTHER")
dates <- c("20190401", "20200401")






#----------------------------------------------------
# Part 1. Prepare polygon and combination table
#----------------------------------------------------

# import polygon
poly <- st_read("data/input/marine_regions/Intersect_EEZ_IHO_v4_2020/Intersect_EEZ_IHO_v4_2020.shp")


# reproject (use clip_to_globe function from OHI-Science)
# https://github.com/OHI-Science/ohiprep_v2020/blob/18880e6ab7eb1607f754660b67d4b3b46f2a6043/globalprep/spp/v2020/_setup/common_fxns.R
poly_moll <- poly %>%
  filter(MRGID_EEZ == 0) %>%
  clip_to_globe() %>%  # ensures no coordinates outside +-180, +-90
  st_transform("+proj=moll")
plot(poly_moll["IHO_SEA"])

# Prepare combinations
combinations <- data.frame(
  var = rep(vars, each=length(dates)),
  date = rep(dates, length(vars)))


#----------------------------------------------------
# Part 2. Calculate occupancy and global average on ABJN
#----------------------------------------------------

# import density maps
dens20 <- raster("data/out/ais-global/density/20200401_COUNT_dens_moll.tif")
dens19 <- raster("data/out/ais-global/density/20190401_COUNT_dens_moll.tif")

# rasterize EEZ
rpoly <- fasterize(poly_moll, dens20)

# calculate total area
total_area <- binsurf(rpoly)

# calculate area of overlap
rdens <- dens20/dens20
ov <- rpoly * rdens
ov_area <- binsurf(ov)
100 * (ov_area / total_area) # 73.89404

# calculate global averages
avg20 <- values(dens20 * rpoly) %>% mean(na.rm=TRUE) 
avg19 <- values(dens19 * rpoly) %>% mean(na.rm=TRUE) 
(avg20-avg19) * 1000  # 0.4400944
100*(avg20-avg19)/avg19  # 2.441447



#----------------------------------------------------
# Part 3. Extract data from shipping rasters
#----------------------------------------------------

# Prepare cluster for parallel computing
cl <- makeCluster(10)
registerDoParallel(cl)

# Extract data for each combination  
data <- bind_rows(foreach(i = 1:nrow(combinations), .packages=c("dplyr", "raster", "exactextractr")) %dopar% { 
  
  # select file
  ivar <- combinations$var[i]
  idate <- combinations$date[i]
  iyear <- substr(idate, 1, 4)
  tif_file <- list.files(input_dir, full.names = TRUE, pattern=sprintf("%s_%s_dens_moll.tif$", idate, ivar))
  
  # import raster
  rast <- raster(tif_file)
  
  # summarize data per polygon
  # weighted by the fraction of each cell that is covered by the polygon
  extr <- exact_extract(rast, poly_moll, "mean")
  
  # generate data.frame
  df <- data.frame(MRGID = poly_moll$MRGID, var = ivar, year = iyear, measure = "dens", value = extr)
  df
})

stopCluster(cl)  # Stop cluster



#----------------------------------------------------
# Part 4. Calculate differences
#----------------------------------------------------

change <- data %>%
  group_by(MRGID, var) %>%
  arrange(year) %>%
  summarize(delta = last(value)-first(value),
            per = 100*(delta/first(value)),
            perlog = 100*log(last(value)/first(value)))

change$MRGID <- as.character(change$MRGID)


#----------------------------------------------------
# Part 5. Plot map
#----------------------------------------------------

# combine with polygon map
poly_moll$MRGID <- as.character(poly_moll$MRGID)
poly_ais <- left_join(poly_moll, change, by = c("MRGID" = "MRGID"))

# transform to data.frame
df <- data.frame(poly_ais)
df <- dplyr::select(df, -geometry)

# import landmask
data(countriesHigh, package = "rworldxtra", envir = environment())
countriesHigh <- spTransform(countriesHigh, CRS("+proj=moll +ellps=WGS84"))
land <- st_as_sf(countriesHigh)

# create box
box <- bb(xmin = -180, xmax = 180, ymin = -90, ymax = 90, crs="+proj=moll +ellps=WGS84")

# common scale to compare with other regions
breaks <- c(-0.8, -0.4, -0.1, -0.05, -0.02, -0.01, -0.005, 0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.4)
labels <- c("-0.8 to -0.4", "-0.4 to -0.1", "-0.1 to -0.05", "-0.05 to -0.02", "-0.02 to -0.01", "-0.01 to -0.005", "-0.005 to 0.00",
            "0.00 to 0.005", "0.005 to 0.01", "0.01 to 0.02", "0.02 to 0.05", "0.05 to 0.1", "0.1 to 0.4")

for(i in 1:length(vars)){
  
  # select variable
  ivar <- vars[i]
  
  # plot variable
  p1 <- tm_shape(filter(poly_ais, var == ivar), projection="+proj=moll +ellps=WGS84") +
    tm_polygons("delta",
                style = "fixed",
                breaks = breaks,
                labels = labels,
                midpoint = 0,
                palette = "RdBu",
                border.col = "grey60",
                border.alpha = 0.3,
                textNA = "No data",
                colorNA = "white") +
    tm_shape(land, projection="+proj=moll +ellps=WGS84") +
    tm_polygons(col="grey60", border.alpha = 0) +
    tm_shape(box) +
    tm_polygons(alpha=0, border.col = "grey60") +
    tm_layout(frame = FALSE, legend.title.size=1, legend.outside = TRUE)
  
  
  # export plot
  out_file <- paste0(out_dir, sprintf("%s_abjn_change.png", ivar))
  tmap_save(tm = p1, filename = out_file, width=25, height=12, units = "cm")
}



#----------------------------------------------------
# Part 6. Table
#----------------------------------------------------

# transform to wide format
data_delta <- dcast(df, IHO_SEA ~ var, value.var = "delta")
data_per <- dcast(df, IHO_SEA ~ var, value.var = "per")

# combine delta with %
data <- merge(data_delta, data_per, by=c("IHO_SEA"))

# filter countries with no data (Caspian Sea)
data <- filter(data, !is.na(COUNT.x))

# Order table
data <- arrange(data, COUNT.x)

# Order columns
data <- dplyr::select(data, "IHO_SEA",
                      "COUNT.x", "COUNT.y",
                      "CARGO.x", "CARGO.y",
                      "TANKER.x", "TANKER.y",
                      "PASSENGER.x", "PASSENGER.y",
                      "FISHING.x", "FISHING.y",
                      "OTHER.x", "OTHER.y")

# set indices oer columns
idx_abs <- c(3, 5, 7, 9, 11, 13)-1
idx_rel <- idx_abs + 1

# transform values to 10^-3
data[,idx_abs] <- data[,idx_abs] * 1000  # select columns

#### Create flextable

# typology
typology <- data.frame(
  col_keys = names(data)[1:13],
  what = c("High Seas", rep(c("All vessels", "Cargo", "Tanker", "Passenger", "Fishing", "Other"),each=2)),
  measure = c("High Seas", rep(c("Absolute change", "Relative change"), 6)),
  units = c("High Seas", rep(c("delta", "(%)"), 6)),
  stringsAsFactors = FALSE)

# create table
ft <- flextable(data)
ft <- set_header_df(ft, mapping = typology, key = "col_keys" )

# merge headers
ft <- merge_h(ft, part = "header")
ft <- merge_v(ft, part = "header")

# set NA values and digits
ft <- colformat_num(ft, j = idx_abs, digits = 2, na_str = "NA")
ft <- colformat_num(ft, j = idx_rel, digits = 1, na_str = "NA")

# label units
ft <- compose(ft, part = "header", i=3, j = idx_abs,
              value = as_paragraph("(x10", as_sup("-3"), " vessels/", "km",as_sup("2"),")"))

# theme of the table
ft <- theme_zebra(ft, odd_header = "transparent", even_header = "transparent")
ft <- color(ft, i = 1:2, color = "#007FA6", part = "header")
ft <- color(ft, i = 3, color = "gray40", part = "header")
ft <- hline(ft, border = fp_border(width = .75, color = "#007FA6"), part = "body" )
ft <- hline(ft, border = fp_border(width = 2, color = "#007FA6"), part = "header" )
ft <- bold(ft, part = "header") # bold header
ft <- bold(ft, i=3, bold=FALSE, part = "header") # bold header

# set properties
ft <- set_table_properties(ft, width = 1, layout = "autofit")  # autofit to width

# add header
ft <- add_header_lines(ft, values = "Suppl. Data 2.")
ft <- fontsize(ft, size = 7, part = "all")


#### Export to docx

# create a new document object
doc=new.word.doc()

# add the report title and an empty line
add.title(doc, "Supplementary Data 2")
add.text(doc, "Average change for each vessel category for each High Seas subregion. Estimated changes correspond to the difference between April 2020, when more stringent lockdown measures took place globally, and April 2019, as a reference period. Data ordered by increasing absolute change for `All vessels` category.")
add.empty.line(doc)

# add table
add.table(doc, ft)

#set the orientation back to portrait
end.landscape(doc)
#start.landscape(doc)
# generate the Word document using the print function
print(doc, target=paste0(out_dir, "suppData2.docx"))
