threat_map <- function(in_dat, threat_category, threat_pal, hybas_habitat_types, proj){
# hybas_habitat_types = p2_hybas_habitat_types_sf
filtered_df <- in_dat |> 
  left_join(hybas_habitat_types) |> 
  st_as_sf() |> 
  dplyr::filter(ThreatCategory == threat_category) |> 
  # remove visual bug with robinson projection
  st_wrap_dateline()

proj_df <- st_transform(filtered_df, crs = st_crs(proj))

if(threat_category == "Habitat"){
  pal <- threat_pal$Habitat_pal
} else if(threat_category == "Pollution"){
  pal <- threat_pal$Pollution_pal
} else if(threat_category == "Exploitation"){
  pal <- threat_pal$Exploitation_pal
} else if(threat_category == "Invasive species"){
  pal <- threat_pal$Invasive_pal
} else if(threat_category == "Climate and weather"){
  pal <- threat_pal$Climate_pal
}

threat_map <- ggplot()+
  geom_sf(data = proj_df, aes(geometry = Shape, fill = MeanWeightedThreatMetric, color = MeanWeightedThreatMetric))+
  scale_fill_gradientn(
    colors = colorRampPalette(c(rev(unlist(pal))))(100),
    limits = c(0, max(proj_df$MeanWeightedThreatMetric, na.rm = T)),
    na.value = "gray80",
    breaks = c(0 + max(proj_df$MeanWeightedThreatMetric, na.rm = T)/10, 
               #max(habitat_data$MeanWeightedThreatMetric)/2, 
               max(proj_df$MeanWeightedThreatMetric, na.rm = T) - max(proj_df$MeanWeightedThreatMetric, na.rm = T)/10),
    labels = c("Lower", "Higher")
  )+
  scale_color_gradientn(
    colors = colorRampPalette(c(rev(unlist(pal))))(100), 
    na.value= "gray80"
  )+
  guides(color = "none")+
  guides(fill = guide_colorbar(title = "Mean Threat",
                               title.position = "top",
                               direction = "horizontal",
                               barwidth = 7,
                               barheight = 1))+
  theme_void()+
  theme(
    #legend.position = c(0.1, 0.21),
    legend.ticks = element_blank(),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 11)
  )

return(threat_map)

}

subThreat_map <- function(in_dat, threat_category, threat_pal, subcat_habitat, subcat_pollution, subcat_climate, proj, hybas_habitat_types){
  
  filtered_df <- in_dat |> 
    left_join(hybas_habitat_types) |> 
    st_as_sf() |> 
    dplyr::filter(ThreatCategory == threat_category) |> 
    # remove visual bug with robinson projection
    st_wrap_dateline()
  
  proj_df <- st_transform(filtered_df, crs = st_crs(proj))
  
  if(threat_category %in% subcat_habitat){
    pal <- threat_pal$Habitat_pal
  } else if(threat_category %in% subcat_pollution){
    pal <- threat_pal$Pollution_pal
  } else if(threat_category == "Overfishing"){
    pal <- threat_pal$Exploitation_pal
  } else if(threat_category == "Invasive non-native species"){
    pal <- threat_pal$Invasive_pal
  } else if(threat_category %in% subcat_climate){
    pal <- threat_pal$Climate_pal
  }
  
  threat_map <- ggplot()+
    geom_sf(data = proj_df, aes(geometry = Shape, fill = MeanWeightedThreatMetric, color = MeanWeightedThreatMetric))+
    scale_fill_gradientn(
      colors = colorRampPalette(c(rev(unlist(pal))))(100),
      limits = c(0, max(proj_df$MeanWeightedThreatMetric, na.rm = T)),
      na.value = "gray80",
      breaks = c(0 + max(proj_df$MeanWeightedThreatMetric, na.rm = T)/10, 
                 #max(habitat_data$MeanWeightedThreatMetric)/2, 
                 max(proj_df$MeanWeightedThreatMetric, na.rm = T) - max(proj_df$MeanWeightedThreatMetric, na.rm = T)/10),
      labels = c("Lower", "Higher")
    )+
    scale_color_gradientn(
      colors = colorRampPalette(c(rev(unlist(pal))))(100), 
      na.value= "gray80"
    )+
    guides(color = "none")+
    guides(fill = guide_colorbar(title = "Mean Threat",
                                 title.position = "top",
                                 direction = "horizontal",
                                 barwidth = 7,
                                 barheight = 1))+
    theme_void()+
    theme(
      #legend.position = c(0.1, 0.21),
      legend.ticks = element_blank(),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 11)
    )
  
  return(threat_map)
  
}
# in_dat = p2_mean_weighted_threats, legend_png = p3_legend_png, threat_category = threat_cat
## note about out_file - need to find a way to save the cowplot version of the legend without crowding the images folder
## maybe by saving the initial legend png in the findex out folder and then saving 
## the cowplot version to the earth-in-flux parent directory "src/assets/images/" folder
cowplot_legend <- function(in_dat, legend_png, threat_category, out_file){
  
  threat_df <- in_dat |> 
    filter(ThreatCategory == threat_category)
  
  min_val <- round(min(threat_df$MeanWeightedThreatMetric, na.rm = T), digits = 2)
  max_val <- round(max(threat_df$MeanWeightedThreatMetric, na.rm = T), digits = 2)
  
  # Define colors
  background_color = NA
  font_color = "#ffffff"
  
  # The background canvas for your viz (DO NOT EDIT)
  canvas <- grid::rectGrob(
    x = 0, y = 0, 
    width = 9, height = 9,
    gp = grid::gpar(fill = background_color, alpha = 1, col = background_color)
  )
  
  # margin for plotting (DO NOT EDIT)
  margin = 0.04
  
  # Load in USGS logo (also a black logo available)
  legend <- magick::image_read(legend_png) 
  
  final_legend <- ggdraw(ylim = c(0,1), # 0-1 scale makes it easy to place viz items on canvas
                         xlim = c(0,1)) +
    # a background  (DO NOT EDIT)
    draw_grob(canvas,
              x = 0, y = 1,
              height = 9, width = 16,
              hjust = 0, vjust = 1) +
    draw_image(legend, 
               x = 0.08,
               y = 0.08,
               width = 0.77, 
               hjust = 0, vjust = 0, 
               halign = 0, valign = 0)+
    # min max values
    draw_label(as.character(min_val),
               x = 0.02,
               y = 0.54,
               hjust = 0,
               vjust = 1,
               lineheight = 0.75,
               color = "gray50",
               size = 9) +
    draw_label(as.character(max_val),
               x = 1,
               y = 0.54,
               hjust = 1,
               vjust = 1,
               lineheight = 0.75,
               color = "gray50",
               size = 9) 
  #429x176
  ggsave(out_file, final_legend, height = 176, width = 429, units = "px", dpi = 300, bg = "transparent")
}