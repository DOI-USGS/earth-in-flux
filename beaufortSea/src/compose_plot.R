compose_timeline <- function(timeline_grob, 
                             assemblage_grob_100, 
                             assemblage_grob_500, 
                             assemblage_grob_1000,
                             assemblage_grob_1500,
                             assemblage_grob_2000,
                             color_scheme, png_out){
  

  showtext::showtext_opts(dpi = 300, regular.wt = 200, bold.wt = 700)
  showtext::showtext_auto(enable = TRUE)
  title_font_size <- 35
  explainer_font_size <- 9
  
  # Define colors
  background_color <- "#ffffff"
  font_color <- color_scheme$Cassidulina
  
  # Define margin for cowplotting
  margin_cowplot <- 0.04
  
  # The background canvas for your viz on mobile
  mobile_height <- 24
  mobile_width <- 6
  canvas_mobile <- grid::rectGrob(
    x = 0, y = 0, 
    width = mobile_width, height = mobile_height,
    gp = grid::gpar(fill = background_color, alpha = 1, col = background_color)
  )
  
  
  ggdraw(ylim = c(0, mobile_height), 
         xlim = c(0, mobile_width)) +
    # a background
    draw_grob(canvas_mobile,
              x = 0, y = 1,
              height = mobile_height, 
              width = mobile_width,
              hjust = 0, 
              vjust = 1) +
    # the main plot
    draw_plot(timeline_grob,
              x = 0,
              y = 0,
              height = mobile_height,
              width = 3.6) +
    draw_text("Year of the sample,\nbeginning 2000 years ago\nat 0 A.D.",
              hjust = 0, vjust = 0,
              x = 0.5,
              y = 23.5,
              size = explainer_font_size, family = annotation_font) +
    draw_text("The horizontal bars and circle sizes show the relative abundance\nof each microfossil genus",
              hjust = 0, vjust = 0,
              x = 0.65,
              y = 23,
              size = explainer_font_size, family = annotation_font) +
    # t = 0
    draw_plot(assemblage_grob_100,
              x = 3.4,
              y = 21,
              height = 3,
              width = 3) +
    draw_text("2000 years ago, the system\nwas dominated by Paracyprideis,\nCassidulina and Elphidium species.",
              hjust = 0, vjust = 0,
              x = 1.5,
              y = 22,
              size = explainer_font_size, family = supporting_font) +
    # t = 500
    draw_text("The first 1000 years demonstrated relative stability",
              hjust = 0, vjust = 0,
              x = 1,
              y = 18,
              size = explainer_font_size, family = supporting_font) +
    draw_plot(assemblage_grob_500,
              x = 3.4,
              y = 16,
              height = 3,
              width = 3) +
    # t = 1000
    draw_text("Then, from 1000-1300, agglutinated species Spiroplectammina\nbecame more abundant as a result of\nclimate instability during the Medieval Climate Anomaly",
              hjust = 0, vjust = 0,
              x = 1,
              y = 12.5,
              size = explainer_font_size, family = supporting_font) +
    draw_plot(assemblage_grob_1000,
              x = 3.1,
              y = 10.5,
              height = 3,
              width = 3) +
    # t = 1500
    draw_text("During the Little Ice Age, Paracyprideis declined\nto its lowest levels in the record",
              hjust = 0, vjust = 0,
              x = 1,
              y = 7,
              size = explainer_font_size, family = supporting_font) +
    draw_plot(assemblage_grob_1500,
              x = 3.1,
              y = 4.8,
              height = 3,
              width = 3) +
    # t = 1500
    draw_text("In the modern record, Spiroplectammina and\nKotoracythere are at all-time highs, and the\npreviously dominant Elphidium and Cassidulina are sparser ",
              hjust = 0, vjust = 0,
              x = 0.5,
              y = 0.3,
              size = explainer_font_size, family = supporting_font) +
    draw_plot(assemblage_grob_2000,
              x = 3.1,
              y = 0.2,
              height = 3,
              width = 3)
  
  ggsave(filename = png_out, 
         width = mobile_width, height = mobile_height, dpi = 300)
  
}