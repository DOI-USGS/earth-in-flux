
#' @description establish ggplot code to create threat maps
#' @param in_dat dataframe with mean weighted threat scores by threat type and HYBAS_ID
#' @param threat_category list of target threat categories
#' @param threat_pal dataframe with color palettes and file name templates by threat type
#' @param hybas_habitat_types shape file with HYBAS IDs and their habitat types
#' @param proj character string with map projection definition
general_threat_map <- function(in_dat, threat_category, threat_pal, hybas_habitat_types, proj){
  filtered_df <- in_dat |> 
    left_join(hybas_habitat_types) |> 
    st_as_sf() |> 
    dplyr::filter(ThreatCategory == threat_category) |> 
    # remove visual bug with robinson projection
    st_wrap_dateline()
  
  proj_df <- st_transform(filtered_df, crs = st_crs(proj))
  
  threat_map <- ggplot()+
    geom_sf(data = proj_df, aes(geometry = Shape, fill = MeanWeightedThreatMetric), color = NA)+
    scale_fill_gradientn(
      colors = colorRampPalette(c(rev(unlist(threat_pal))))(100),
      limits = c(0, max(proj_df$MeanWeightedThreatMetric, na.rm = T)),
      na.value = "gray80",
      breaks = c(0 + max(proj_df$MeanWeightedThreatMetric, na.rm = T)/10, 
                 max(proj_df$MeanWeightedThreatMetric, na.rm = T) - max(proj_df$MeanWeightedThreatMetric, na.rm = T)/10)
    )+
   # scale_color_gradientn(
  #    colors = colorRampPalette(c(rev(unlist(threat_pal))))(100), 
  #    na.value= "gray80"
  #  )+
   # guides(color = "none")+
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

#' @description filter data to target major threat and apply to established ggplot function
#' @param in_dat dataframe with mean weighted threat scores by threat type and HYBAS_ID
#' @param threat_category list of target threat categories
#' @param threat_pal dataframe with color palettes and file name templates by threat type
#' @param hybas_habitat_types shape file with HYBAS IDs and their habitat types
#' @param proj character string with map projection definition
threat_map <- function(in_dat, threat_category, threat_pal, hybas_habitat_types, proj){
  
    pal <- threat_pal |> 
      filter(MajorCat == threat_category) |> 
      select(pal) |> 
      unique()
  
  final_plot <- general_threat_map(in_dat = in_dat, 
                                   threat_category = threat_category, 
                                   threat_pal = pal,
                                   hybas_habitat_types = hybas_habitat_types,
                                   proj = proj)

}

#' @description filter data to target sub threat and apply to established ggplot function
#' @param in_dat dataframe with mean weighted threat scores by threat type and HYBAS_ID
#' @param threat_category list of target threat categories
#' @param threat_pal dataframe with color palettes and file name templates by threat type
#' @param hybas_habitat_types shape file with HYBAS IDs and their habitat types
#' @param proj character string with map projection definition
subThreat_map <- function(in_dat, threat_category, threat_pal, proj, hybas_habitat_types){ 
  
  pal <- threat_pal |> 
    filter(ThreatCategory == threat_category) |> 
    select(pal) |> 
    unique()
  
  final_plot <- general_threat_map(in_dat = in_dat, 
                                   threat_category = threat_category, 
                                   threat_pal = pal,
                                   hybas_habitat_types = hybas_habitat_types,
                                   proj = proj)
  
}

#' @description cowplot code to style the legend
#' @param in_dat dataframe with mean weighted threat scores by threat type and HYBAS_ID
#' @param legend_png character string of raw extracted legend file location
#' @param threat_category list of target threat categories
#' @param out_file character string of file name and location to save final legend
#' @param height png height
#' @param width png width
#' @param unit png height and width units
#' @param dpi png dpi
cowplot_legend <- function(in_dat, legend_png, threat_category, out_file, height, width, unit, dpi){
  
  threat_df <- in_dat |> 
    filter(ThreatCategory == threat_category)
  
  min_val <- round(min(threat_df$MeanWeightedThreatMetric, na.rm = T), digits = 2)
  max_val <- round(max(threat_df$MeanWeightedThreatMetric, na.rm = T), digits = 2)
  
  # Define colors
  background_color = NA
  
  # The background canvas for your viz (DO NOT EDIT)
  canvas <- grid::rectGrob(
    x = 0, y = 0, 
    width = 9, height = 9,
    gp = grid::gpar(fill = background_color, alpha = 1, col = background_color)
  )
  
  # Load raw legend png
  legend <- magick::image_read(legend_png) 
  
  final_legend <- ggdraw(ylim = c(0,1), # 0-1 scale makes it easy to place viz items on canvas
                         xlim = c(0,1)) +
    # a background  (DO NOT EDIT)
    draw_grob(canvas,
              x = 0, y = 1,
              height = 9, width = 16,
              hjust = 0, vjust = 1) +
    draw_image(legend, 
               x = 0.079,
               y = 0.35,
               width = 0.74, 
               hjust = 0, vjust = 0, 
               halign = 0, valign = 0)+
    # min max values
    draw_label(as.character(min_val),
               x = 0.02,
               y = 0.55,
               hjust = 0,
               vjust = 1,
               lineheight = 0.75,
               color = "gray50",
               size = 8) +
    draw_label(as.character(max_val),
               x = 0.99,
               y = 0.55,
               hjust = 1,
               vjust = 1,
               lineheight = 0.75,
               color = "gray50",
               size = 8) +
    # higher lower labels
    draw_label("Lower",
               x = 0.05,
               y = 0.11,
               hjust = 0,
               vjust = 0,
               lineheight = 0.75,
               color = "black",
               size = 8) +
    draw_label("Higher",
               x = 0.86,
               y = 0.11,
               hjust = 1,
               vjust = 0,
               lineheight = 0.75,
               color = "black",
               size = 8)

  ggsave(out_file, final_legend, height = height, width = width, units = unit, dpi = dpi, bg = "transparent")
}

#' @description cowplot code to style the legend
#' @param type "threat" or "subThreat"
#' @param plot ggplot map with legend to be saved
#' @param threat_category list of target threat categories
#' @param in_dat dataframe with mean weighted threat scores by threat type and HYBAS_ID
#' @param threat_pal dataframe with color palettes and file name templates by threat type
#' @param height png height
#' @param width png width
#' @param unit png height and width units
#' @param dpi png dpi
save_legend <- function(type, plot, threat_category, in_dat, threat_pal, height, width, unit, dpi){
  
  if(type == "threat"){
    name_conv <- threat_pal |> 
      filter(MajorCat == threat_category)
    
    plot_legend <- get_plot_component(plot, "guide-box-right", return_all = T)
    
    out_file <- sprintf(unique(name_conv$threat_legend_raw), str_replace_all(threat_category, " ", "_"))
    
    ggsave(out_file, 
           plot_legend, dpi = dpi, bg = "transparent")
    knitr::plot_crop(out_file)
    
    out_file_final <- sprintf(unique(name_conv$threat_legend), str_replace_all(threat_category, " ", "_"))
    
    cowplot_legend(in_dat = in_dat, legend_png = out_file, threat_category = threat_category, 
                   out_file = out_file_final, height = height, width = width, unit = unit, dpi = dpi)
    
    return(out_file_final)
    
  } else if(type == "subThreat"){
    name_conv <- threat_pal |> 
      filter(ThreatCategory == threat_category)
    
    plot_legend <- get_plot_component(plot, "guide-box-right", return_all = T)
    
    out_file <- sprintf(unique(name_conv$subThreat_legend_raw), str_replace_all(threat_category, " ", "_"))
    
    ggsave(out_file, 
           plot_legend, dpi = dpi, bg = "transparent")
    knitr::plot_crop(out_file)
    
    out_file_final <- sprintf(unique(name_conv$subThreat_legend), str_replace_all(threat_category, " ", "_"))
    
    cowplot_legend(in_dat = in_dat, legend_png = out_file, threat_category = threat_category, 
                   out_file = out_file_final, height = 176, width = 429, unit = "px", dpi = dpi)
    
    return(out_file_final)
    
  }
}

#' @description cowplot code to style the legend
#' @param type "threat" or "subThreat"
#' @param plot ggplot map to be saved
#' @param threat_category list of target threat categories
#' @param threat_pal dataframe with color palettes and file name templates by threat type
#' @param height png height
#' @param width png width
#' @param unit png height and width units
#' @param dpi png dpi
save_map <- function(type, plot, threat_category, threat_pal, height, width, dpi){
  
  if(type == "threat"){
    name_conv <- threat_pal |> 
      filter(MajorCat == threat_category)
    
    out_file <-  sprintf(unique(name_conv$threat_map), str_replace_all(threat_category, " ", "_"))
    
    ggsave(out_file, 
           plot, height = height, width = width, dpi = dpi)
    
    knitr::plot_crop(out_file)
    
  } else if(type == "subThreat"){
    name_conv <- threat_pal |> 
      filter(ThreatCategory == threat_category)
    
    out_file <- sprintf(unique(name_conv$subThreat_map), str_replace_all(threat_category, " ", "_"))
    
    ggsave(out_file, 
           plot, height = height, width = width, dpi = dpi)
    
    knitr::plot_crop(out_file)
    
  }
}

#' @description create top threats by basin global map
#' @param in_dat dataframe with mean weighted threat scores by threat type and HYBAS_ID
#' @param threat_pal dataframe with color palettes and file name templates by threat type
#' @param hybas_habitat_types shape file with HYBAS IDs and their habitat types
#' @param proj character string with map projection definition
#' @param threat_category list of target threat categories
top_threat_plot <- function(in_dat, threat_pal, hybas_habitat_types, proj, threat_category){
  
  processed_df <- in_dat |> 
    group_by(HYBAS_ID) |> 
    filter(MeanWeightedThreatMetric == max(MeanWeightedThreatMetric, na.rm = T))
  
  processed_sf <- processed_df |> 
    left_join(hybas_habitat_types) |> 
    st_as_sf() |> 
    # remove visual bug with robinson projection
    st_wrap_dateline()
  
  proj_sf <- st_transform(processed_sf, crs = st_crs(proj))
  
  # make non-target threat category values NA so they are not plotted
  if(threat_category != "none"){
    proj_sf <- proj_sf |> 
      mutate(ThreatCategory = case_when(ThreatCategory != threat_category ~ NA, .default = as.character(ThreatCategory)))
  }
  
  pal <- threat_pal |> 
    select(MajorCat, pal) |> 
    rowwise() |> 
    mutate(pal = first(pal)) |> 
    unique() |> 
    mutate(pal = case_when(pal == "#4E6D6E" ~ "#598586",
                           pal == "#7A562B" ~ "#A97639",
                           pal == "#835192" ~ "#995EAB",
                           pal == "#B74F49" ~ "#963C36",
                           pal == "#002D5E" ~ "#002D5E")) 
  
  threat_map <- ggplot()+
    geom_sf(data = proj_sf, aes(geometry = Shape, fill = ThreatCategory), color = NA)+
    scale_fill_manual(values = pal$pal, breaks = pal$MajorCat, na.value = "gray80")+
    guides(fill = guide_legend(nrow = 2,)) +
    theme_void()+
    theme(
      legend.ticks = element_blank(),
      legend.title = element_text(face = "bold"),
      legend.title.position = "top",
      legend.text = element_text(size = 11),
      legend.direction = "horizontal"
    )
}

#' @description save the legend for top threat by basin global map
#' @param plot ggplot map to be saved
#' @param dpi png dpi
#' @param out_file file path to save the legend to
save_top_threat_legend <- function(plot, dpi, out_file){
  
  plot_legend <- get_plot_component(plot, "guide-box-right", return_all = T)
  
  ggsave(out_file, 
         plot_legend, dpi = dpi, bg = "transparent")
  
  knitr::plot_crop(out_file)
}

#' @description create and save a thumbnail for the findex global maps page
#' @param threat_category list of target threat categories
#' @param in_dat dataframe with mean weighted threat scores by threat type and HYBAS_ID
#' @param threat_pal dataframe with color palettes and file name templates by threat type
#' @param hybas_habitat_types shape file with HYBAS IDs and their habitat types
#' @param proj character string with map projection definition
#' @param height png height
#' @param width png width
#' @param unit png height and width units
#' @param dpi png dpi
top_threat_thumbnail <- function(in_dat, threat_pal, hybas_habitat_types, proj, threat_category, height, width, dpi, out_file){
  
  # these coordinates are NAD83 coordinates that correspond to "EPSG:4269", so will only work with the other data if they are also in "EPSG:4269"
  bbox <- st_bbox(c(xmin = -25, xmax = 60, ymin = -35, ymax = 40))

  # get bbox of projected area of interest
  if (!(proj == "EPSG:4269")) {
    bbox <- sf::st_as_sfc(bbox) |> 
      st_as_sf(crs = "EPSG:4269") |>
      sf::st_transform(crs = proj) |>
      sf::st_bbox()
  }
  
  # call custom plotting function
  threat_map <- top_threat_plot(in_dat = in_dat, 
                                threat_pal = threat_pal, 
                                hybas_habitat_types = hybas_habitat_types, 
                                proj = proj,
                                threat_category = threat_category)  + 
    # crop for thumbnail
    coord_sf(xlim = c(bbox['xmin'], bbox['xmax']), 
             ylim = c(bbox['ymin'], bbox['ymax']), 
             expand = F) +
    # white background and add margin
    theme(legend.position = "none",
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA),
          plot.margin = margin(t = 2,
                               r = 2,
                               b = 2,
                               l = 2,
                               unit = "cm")) 
  
  # save in square ratio
  ggsave(out_file, 
         threat_map, height = height, width = width, dpi = dpi)
  
}
