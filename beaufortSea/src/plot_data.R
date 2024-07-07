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
    scale_size_continuous(range = c(1,6)) +
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
    ylab("Year") +
    xlab("Relative Abundance")
}
