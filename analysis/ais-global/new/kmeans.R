library(TSrepr)
library(ggplot2)
library(data.table)






change2 <- change %>%
  dcast(ISO_SOV1 ~ month, value.var=c("per"))

change3 <- as.matrix(change2[, 2:8])  

res_km <- kmeans(change3, 5, nstart = 10)

change2$class <- as.character(res_km$cluster)


change4 <- change %>%
  left_join(dplyr::select(change2, ISO_SOV1, class), by="ISO_SOV1")



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
ggplot(change4) +
  facet_wrap(~class, ncol = 2) +
  geom_line(aes(month, per, group = ISO_SOV1, colour=class),  size=0.7, alpha = 0.3) +
  geom_line(data = centers, aes(Time, value, colour=class), alpha = 0.90, size = 1) +
  scale_color_brewer(palette="Set2") +
  geom_hline(yintercept = 0, linetype="dotted") +
  labs(x = "Time", y = "Load (normalised)") +
  theme_article()

# plot map

# import country map
library(tmap)
data("World")

# combine with countries with EEZ by ISO code
map_data <- left_join(World, change4, by = c("iso_a3" = "ISO_SOV1"))
map_data <- map_data %>% filter(name != "Antarctica")

# plot data (average)
p1 <- tm_shape(map_data) +
  tm_polygons("class",
              title=("Cluster"),
              palette = "Set2", #brewer.ylgnbu(10) # YlGnBu, -RdYlBu
              border.col = "grey60", 
              border.alpha = 0.3) +
  tm_layout(legend.title.size=1)

# export plot
#out_file <- paste0(output_plot, "/april_global_map_med.png")
#tmap_save(tm = p1, filename = out_file, width=22, height=10, units = "cm")

https://cran.r-project.org/web/packages/TSrepr/vignettes/TSrepr_representations_use_case.html

