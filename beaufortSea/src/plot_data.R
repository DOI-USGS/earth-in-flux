plot_assemblages <- function(data_in, year){
  
  # set up data for plotting using plotcircles package syntax
  year_df <- data_in |> filter(decade == year)
  year_pack <- circleProgressiveLayout(year_df$mean_abundance, sizetype = "area")
  year_join <- cbind(year_df, year_pack)
  year_plot <- circleLayoutVertices(year_pack, npoints = 50)
  year_plot$species <- rep(year_df$species, each = 51)
  year_plot$color <- rep(year_df$hexcode, each = 51)
  year_plot$name <- rep(year_df$name, each = 51)
  year_plot$genus <- rep(year_df$genus, each = 51)
  
  ggplot() + 
    # Make the bubbles
    geom_polygon(data = year_plot, 
                 aes(x, y, group = id, fill = color), 
                 colour = "grey80", linewidth = 0.3) +
    # Add text in the center of each bubble + control its size
    geom_text(data = year_join, aes(x, y, size = mean_abundance, label = name)) +
    scale_size_continuous(range = c(1,3)) +
    scale_fill_identity() +
    # General theme:
    theme_void() + 
    theme(legend.position="none") +
    coord_equal()
}


plot_timeline <- function(data_in){
  ggplot(data = data_in, aes(y = decade, x = pct_abundance, fill = hexcode)) +
    geom_bar(position = "stack", stat = "identity", orientation = "y",
             width = 15, color = NA) +
    scale_fill_identity() +
    scale_y_reverse() +
    theme_minimal() +
    theme(axis.title = element_blank())
}


plot_species_trend <- function(data_in, species_name){
  # Filter data to the species
  data_species <- data_in |> filter(epithet == {{ species_name }}) 
  
  image_path <- sprintf("images/%s", unique(data_species$image_name))
  
  hexcode <- as.character(unique(data_species$hexcode))
  
  ggplot(data_species, aes(x = decade, y = pct_abundance)) +
    geom_smooth(method = "lm", se = FALSE, color = hexcode) +
    geom_image(image = image_path,
               aes(size = I(pct_abundance/100))) +
    ylim(c(0, max(data_species$pct_abundance, na.rm = TRUE)))+
    ylab("Relative Abundance (%)") +
    xlab("Year (A.D.)") +
    theme_bw() +
    theme(legend.position = "none")
}

save_plot <- function(plot_grob, save_name, width, height){
  ggsave(plot = plot_grob,
         file = save_name,
         width = width, height = height, unit = "px")
  return(save_name)
}