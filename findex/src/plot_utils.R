threat_map <- function(in_dat, threat_category, threat_pal){

filtered_df <- st_as_sf(in_dat) |> 
  dplyr::filter(ThreatCategory == threat_category) |> 
  # remove visual bug with robinson projection
  st_wrap_dateline()

proj_df <- st_transform(filtered_df, crs = st_crs("ESRI:54030"))

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
    legend.ticks = element_blank()
  )

return(threat_map)

}

subThreat_map <- function(in_dat, threat_category, threat_pal){
  
  filtered_df <- st_as_sf(in_dat) |> 
    dplyr::filter(Threat == threat_category) |> 
    # remove visual bug with robinson projection
    st_wrap_dateline()
  
  proj_df <- st_transform(filtered_df, crs = st_crs("ESRI:54030"))
  
  if(threat_category %in% c("Dams", "Wetland drainage", "Deforestation and associated runoff",
                            "Riparian degradation", "Agricultural extraction", "Urban extraction", 
                            "Industrial extraction")){
    pal <- threat_pal$Habitat_pal
  } else if(threat_category %in% c("Agricultural effluents", "Urban wastewater", 
                                   "Industrial effluents", "Aquaculture effluents",
                                   "Pharmaceuticals", "Oil or gas exploration",
                                   "Plastics", "Mining")){
    pal <- threat_pal$Pollution_pal
  } else if(threat_category == "Overfishing"){
    pal <- threat_pal$Exploitation_pal
  } else if(threat_category == "Invasive non-native species"){
    pal <- threat_pal$Invasive_pal
  } else if(threat_category %in% c("Change in water temperature", "Drought", "Change in flooding", 
                                 "Change in wind patterns", "Change in ice cover")){
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
      legend.ticks = element_blank()
    )
  
  return(threat_map)
  
}