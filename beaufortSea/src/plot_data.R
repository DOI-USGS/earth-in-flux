
plot_species_trend <- function(data_in, species_name){
  # Filter data to the species
  data_species <- data_in |> filter(epithet == {{ species_name }}) 
  
  image_path <- sprintf("images/in/%s", unique(data_species$image_name))
  
  hexcode <- as.character(unique(data_species$hexcode))
  
  image_size <- unique(data_species$image_size)
  
  image_y <- unique(data_species$image_y)
  
  ylim <- unique(data_species$ylim)
  
  ggplot(data_species, aes(x = year, y = value)) +
    geom_point(color = hexcode) +
    geom_image(image = image_path, size = image_size,
               y = image_y, x = 150) +
    #ylim(c(0, ylim)) +
    ggtitle("Relative Abundance (%)") +
    xlab("Year (A.D.)") +
    scale_y_continuous(limit = c(0, ylim), 
                       breaks = c(0, 25, 50, 75, 100)) +
    #ggtitle(name) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.title.y = element_blank(),
          plot.title = element_text(family = annotation_font, size = 20, angle = 0, hjust = -0.080),
          axis.title.x = element_text(family = annotation_font, size = 20),
          panel.grid.major = element_line(linewidth = 0.4),
          panel.grid.minor = element_line(linewidth = 0.4))
}

save_plot <- function(plot_grob, save_name, width, height){
  ggsave(plot = plot_grob,
         file = save_name,
         width = width, height = height, unit = "px")
  return(save_name)
}