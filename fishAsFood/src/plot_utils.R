generate_continent_map <- function(country_data, out_file) {
  
  countries_sf <- rnaturalearth::ne_countries(returnclass = 'sf')
  
  country_data_sf <- countries_sf %>%
    dplyr::right_join(country_data)
  
  fish_continents <- unique(country_data_sf[["continent"]])
  fish_continent_sf <- countries_sf |>
    filter(continent %in% fish_continents) |>
    sf::st_transform(crs = "ESRI:54030")
  
  colors <- c("#648E8E","#845c93","#b24d4b","#a27846","#1D3867","#899bb7")
  names(colors) <- c("Africa", "Oceania", "North America", "Europe", "Asia", 
                     "South America")
  
  ggplot() +
    geom_sf(data = fish_continent_sf,
            aes(fill = continent),
            color = NA)  +
    scale_fill_manual(values = colors) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    theme_void() +
    theme(legend.position = "None")
  
  ggsave(out_file, height = 1, width = 2, dpi = 300)
}
