#### Cosmpolitan lineage in HCMC - Map plot #### 

library(ggplot2)
library(patchwork)

## Import HCMC+ shapefile
map <- sf::st_read("VNM_ADM_1.geojson")

# Add centroids
map$mid <- sf::st_centroid(map$geometry)

# Subset only south Vietnam
map_cosmo <- sf::st_crop(map, xmin = 105.5, xmax = 108, ymin = 10, ymax = 12.5)

# Add marker column for provinces with Cosmopolitan genomes
for(i in 1:nrow(map_cosmo)){
  map_cosmo$genome[i] <- ifelse(
    map_cosmo$NAME_1[i] %in% cosmo_df$location,
    "Sampled",
    "Unsampled"
  )
}

## Import Cosmopolitan clade info
cosmo <- readr::read_csv("Sample list_1.csv") |>
  dplyr::filter(is.na(Location) == FALSE)

# Clean names
cosmo$Location[cosmo$Location == "HCMC"] <- "Hồ Chí Minh"
cosmo$Location[cosmo$Location == "TÂY NINH"] <- "Tây Ninh"

# Summarise clades per province
cosmo_df <- as.data.frame(unclass(table(cosmo[,c(1,3)]))) |>
  tibble::rownames_to_column() |>
  dplyr::rename(location = rowname)

# Add centroid for each location in Cosmo data frame
for(i in 1:nrow(cosmo_df)){
  cosmo_df$long[i] <-
    sf::st_coordinates(map$mid[map$NAME_1 == cosmo_df$location[i]])[,1]
}

for(i in 1:nrow(cosmo_df)){
  cosmo_df$lat[i] <-
    sf::st_coordinates(map$mid[map$NAME_1 == cosmo_df$location[i]])[,2]
}

## Plots
cols <- c("#C97064", "#BCA371", "#A6B07E", "#f2f2f2", "#cccccc")

main <- ggplot(map_cosmo) +
  annotate("rect", xmin = 106.5, xmax = 108, ymin = 10, ymax = 11,
           alpha = 0.5, fill = "#a3cef1") +
  annotate("rect", xmin = 105.5, xmax = 106.5, ymin = 10, ymax = 12.5,
           alpha = 0.8, fill = "#7f7f7f") +
  annotate("rect", xmin = 106.5, xmax = 108, ymin = 11, ymax = 12.5,
           alpha = 0.8, fill = "#7f7f7f") +
  geom_sf(aes(fill = genome), colour = "#7f7f7f", linewidth = 0.4) +
  geom_label(data = cosmo_df,
            aes(x = long, y = lat, label = location, fontface = 2),
            nudge_y = -0.19, size = 3) +
  annotate("rect", xmin = 105.5, xmax = 108, ymin = 10, ymax = 12.5,
           color = "black", fill = NA, linewidth = 0.6) +
  scatterpie::geom_scatterpie(aes(x = long, y = lat, group = location, r = 0.12),
                              data = cosmo_df,
                              cols = LETTERS[1:3]) +
  scale_fill_manual(values = cols,
                    name = "Cosmopolitan\nclade",
                    breaks = c("A", "B", "C", "Sampled", "Unsampled"),
                    labels = c("Clade A", "Clade B", "Clade C", "", ""),
                    guide =
                      guide_legend(override.aes =
                                     list(linetype = c(1, 1, 1, 0, 0),
                                          fill = c(cols[1:3], "white", "white")
                                          ))) +

  theme_void()


inset <-ggplot(map) +
  annotate("rect", xmin = 102, xmax = 110, ymin = 8, ymax = 24,
           color = "black", fill = "#f2f2f2", linewidth = 0.4) +
  geom_sf(colour = "#7f7f7f", fill = "white", linewidth = 0.2) +
  annotate("rect", xmin = 105.5, xmax = 108, ymin = 10, ymax = 12.5,
           color = "darkred", fill = NA, linewidth = 0.6) +
  theme_void() +
  theme()

plot <- main + inset_element(inset, -0.5, 0.4, 0.5, 1)

## Export plot
ggsave(plot = plot, "Cosmo_HCMC_map.png", dpi = 300,
       height = 7.15, width = 10, bg = "white")
ggsave(plot = plot, "Cosmo_HCMC_map.pdf", dpi = 300,
       height = 7.15, width = 10, bg = "white")
