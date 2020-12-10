# choke_points


# import chokepoints
# generate buffer
# overlap with cells
# extract density and delta values at each time step and vessel type.
# reuse summary


library(sf)
library(raster)
library(doParallel)
library(foreach)
library(parallel)
library(grid)
library(gridExtra)
library(dplyr)
library(ggplot2)
library(lubridate)
library(egg)
library(reshape2)

# set input dir
input_dir <- "data/out/ais-global/density/"

## Set output plots
output_data <- "results/chokepoints"
if (!dir.exists(output_data)) dir.create(output_data, recursive = TRUE)


# import chockepoints
choke_file <- "data/input/chokepoints/chokepoints.csv"
choke <- read.csv(choke_file)

# add nudge for plot
choke$x_nudge <- c(-40, 0, 0, 40, 40, 0, 50, 0, -40, 0, 0, 0)
choke$y_nudge <- c(0, -20, -20, 30, 20, -20, 0, -40, 20, 40, 20, 20)

# reorder levels by longitude
#levels(choke$name) <- choke$name[order(choke$lon)]
#reorder(choke$name, order(choke$lon))
choke$name <- factor(choke$name, levels = choke$name[order(choke$lon)])


# remove chokepoints with smaller densities
choke <- filter(choke, !name %in% c("Tsugaru Strait", "Yucatan Channel"))

# create spatial buffer
radius <- 20000  # in meters
choke_sf <- st_as_sf(choke, coords = c("lon", "lat"), crs = 4326)
#choke_proj <- st_transform(choke_sf, crs = st_crs('+proj=moll'))
#choke_buff <- st_buffer(choke_proj, radius)
choke_buff <- st_buffer(choke_sf, 0.5)

# export geopackage
st_write(choke_buff, "data/input/chokepoints/chokepoints_buff_geo.gpkg", append = FALSE)



#------------------------------------------------------------------------------
# Part 1. Plot location of chokepoints
#------------------------------------------------------------------------------

# prepare world map
world <- map_data("world")
world.sf <- sf::st_as_sf(world, coords = c("long", "lat"), crs = 4326) %>% 
  group_by(group) %>% 
  summarize(do_union = FALSE) %>%
  st_cast("POLYGON") %>% 
  ungroup()

# plot world map
world <- ggplot() +
  geom_sf(data = world.sf, colour = "gray90", fill = "gray90") + 
  theme(panel.background = element_rect(fill = 'white'))

# plot chokepoints
p1 <- world +
  geom_sf(data = choke_sf, colour = "#3182bd", alpha = .8, size = 4) +
  ggrepel::geom_text_repel(
    data = choke_sf,
    aes(label = name, geometry = geometry),
    stat = "sf_coordinates",
    min.segment.length = 0,
    colour = "grey30",
    segment.color = "grey30",
    size=4,
    nudge_x = choke_sf$x_nudge,
    nudge_y = choke_sf$y_nudge) +
  coord_sf(ylim = c(-60, 90), datum = NA, expand = FALSE) +
  xlab("") + ylab("")

out_file <- paste0(output_data, "/map.png")
ggsave(out_file, p1, width=20, height=10, units = "cm")



#----------------------------------------------------
# Part 2. Prepare chokepoints and combination table
#----------------------------------------------------


# select variables to summarize
vars <- c("COUNT", "FISHING", "PASSENGER", "CARGO", "TANKER", "OTHER")

# select dates
dates <- c(
  seq.Date(as.Date("2019-01-01"), as.Date("2019-07-01"), by = "month"),
  seq.Date(as.Date("2020-01-01"), as.Date("2020-07-01"), by = "month")
)

# transform chokepoint into Mollweide
poly_moll <- choke_buff %>% st_transform("+proj=moll")

# Prepare combinations
combinations <- data.frame(
  var = rep(vars, each=length(dates)),
  date = rep(dates, length(vars)))



#----------------------------------------------------
# Part 3. Extract data from shipping rasters
#----------------------------------------------------

# Prepare cluster for parallel computing
cl <- makeCluster(12)
registerDoParallel(cl)

# Extract data for each combination  
data <- bind_rows(foreach(i = 1:nrow(combinations), .packages=c("dplyr", "raster", "exactextractr")) %dopar% { 
  
  # select file
  ivar <- combinations$var[i]
  idate <- combinations$date[i]
  #iyear <- substr(idate, 1, 4)
  tif_file <- list.files(input_dir, full.names = TRUE, recursive = TRUE, pattern=sprintf("%s_%s_dens_moll.tif$", format(idate, "%Y%m%d"), ivar))
  
  # import raster
  rast <- raster(tif_file)
  
  # summarize data per polygon
  # weighted by the fraction of each cell that is covered by the polygon
  extr <- exact_extract(rast, poly_moll, "mean")

  # generate data.frame
  df <- data.frame(id = poly_moll$id, var = ivar, date = idate, measure = "dens", value = extr)
  df
})

stopCluster(cl)  # Stop cluster

# combine with dataset
data <- merge(choke, data, by = "id")




#----------------------------------------------------
# Part 4. Plot densities
#----------------------------------------------------

data$month <- month(data$date)
data$year <- as.character(year(data$date))

# plot by type
p <- ggplot(filter(data, var == "COUNT", month <= 7), aes(x = month)) +
  geom_line(aes(y = value, color = year), size = 1) +
  scale_color_manual(values=c('#9ecae1', "#3182bd"))+
  scale_x_continuous( breaks = c(1,3,5), labels = month.abb[c(1, 3,5)]) +  # , labels = month.abb[1:6]
  ylab("") + xlab("") +
  facet_wrap(name ~ ., ncol = 5, scales = "free") + # 
  theme_article() +
  theme(legend.position = "none") +
  guides(fill = FALSE)




#----------------------------------------------------
# Part 4. Calculate differences
#----------------------------------------------------


# transform from long to wide format
data <- dcast(data, id + name + var + month ~ year, value.var = "value")

# calculate change metrics
change <- data %>%
          mutate(
            delta = `2020` - `2019`,
            per = 100*(delta/`2019`),
            perlog = 100*log(`2020`/`2019`),
            change_positive = per > 0
          )

change$var <- factor(change$var, levels = c("COUNT", "CARGO", "TANKER", "PASSENGER", "FISHING", "OTHER"))


# plot all vessels

p2 <- ggplot(filter(change, var == "COUNT", month <= 7), mapping=aes(x = month, y = delta, fill = change_positive)) +
  geom_col(alpha=1, width=0.8) +
  #ylab("Relative change (%) in traffic density")+
  ylab(expression(Absolute~change~(Delta~vessels~km^-2~month^-1))) +
  xlab("Month") +
  #scale_x_continuous( breaks = c(1,3,5), labels = month.abb[c(1, 3,5)]) +  # , labels = month.abb[1:6]
  scale_x_continuous(breaks = 1:7) +
  scale_fill_manual(values=c("#e34a33", "#9ecae1")) +
  facet_wrap(name ~ ., ncol = 5) +
  theme_article() +
  guides(fill = FALSE)

out_file <- paste0(output_data, "/monthly_delta.png")
ggsave(out_file, p2, width=19, height=8, units = "cm")




library(ggforce)
p3 <- ggplot(filter(change,  var != "COUNT", month <= 7), mapping=aes(x = month, y = delta, fill = change_positive)) +
  geom_col(alpha=1, width=0.8) +
  #ylab("Relative change (%) in traffic density") + xlab("Month") +
  ylab(expression(Absolute~change~(Delta~vessels~km^-2~month^-1))) +
  #scale_x_continuous( breaks = c(1,3,5), labels = month.abb[c(1, 3,5)]) +  # , labels = month.abb[1:6]
  scale_x_continuous(breaks = 1:7) +
  scale_fill_manual(values=c("#e34a33", "#9ecae1")) +
  facet_wrap(name ~ var, ncol = 5) +
  theme_article() +
  theme(strip.text.x = element_blank(), strip.background = element_blank()) +
  guides(fill = FALSE)

out_file <- paste0(output_data, "/monthly_change.png")
ggsave(out_file, p3, width=16, height=18, units = "cm")

# 
# p <- ggplot(filter(change, var != "COUNT", month <= 6), mapping=aes(x = month, y = perlog, fill = change_positive)) +
#   geom_col(alpha=1, width=0.9) +
#   ylab("") + xlab("") +
#   #geom_vline(xintercept = as.Date("2020-03-11"), linetype="dotted") +
#   #scale_x_date(date_breaks = "2 month", date_labels = "%b") +
#   scale_x_continuous( breaks = c(1,3,5), labels = month.abb[c(1, 3,5)]) +  # , labels = month.abb[1:6]
#   scale_fill_manual(values=c("#e34a33", "#9ecae1")) +
#   facet_wrap(var ~ name, ncol = 9) +
#   theme_article() +
#   guides(fill = FALSE)
# 
# 
# p2 <- ggplot(filter(change, month <= 6), mapping=aes(x = month)) +
#   geom_line(aes(y = delta, color = var), size = 1) +
#   ylab("") + xlab("") +
#   #geom_vline(xintercept = as.Date("2020-03-11"), linetype="dotted") +
#   #scale_x_date(date_breaks = "2 month", date_labels = "%b") +
#   scale_x_continuous( breaks = c(1,3,5), labels = month.abb[c(1, 3,5)]) +  # , labels = month.abb[1:6]
#   scale_fill_manual(values=c("#e34a33", "#9ecae1")) +
#   facet_wrap(name ~ ., ncol = 4) +
#   theme_article() +
#   guides(fill = FALSE)




