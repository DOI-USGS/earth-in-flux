threat_map <- function(in_dat, threat_category, threat_pal, out_file){

filtered_df <- st_as_sf(in_dat) |> 
  dplyr::filter(ThreatCategory == threat_category)

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
  geom_sf(data = filtered_df, aes(geometry = Shape, fill = MeanWeightedThreatMetric, color = MeanWeightedThreatMetric))+
  scale_fill_gradientn(
    colors = colorRampPalette(c(rev(unlist(pal))))(100),
    limits = c(0, max(filtered_df$MeanWeightedThreatMetric)),
    na.value=NA
  )+
  scale_color_gradientn(
    colors = colorRampPalette(c(rev(unlist(pal))))(100), 
    na.value=NA
  )+
  guides(color = "none")+
  guides(fill = guide_colorbar(title.position = "top",
                               direction = "vertical",
                               barwidth = 1,
                               barheight = 5))+
  theme_void()+
  theme(
    legend.position = c(0.1, 0.21),
    legend.title = element_blank()
  )

ggsave(out_file, height = 6, width = 10, dpi = 300) 

}