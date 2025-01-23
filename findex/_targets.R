library(targets)
library(tarchetypes)

options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("tidyverse", 
                            "sf",
                            "readxl",
                            "cowplot"))

source("src/data_utils.R")
source("src/plot_utils.R")

p1 <- list(
  # Download from SB using SB link with token (see project notes)
  tar_target(
    p1_threats_csv,
    "in/Global_Fisheries_Threats.csv",
    format = "file"
  ),
  tar_target(
    p1_weights_csv,
    "in/Final_Weights.csv",
    format = "file"
  ),
  # level 4 hydrobasins
  # with attributes: https://www.hydrosheds.org/hydroatlas
  # Global BasinATLAS in geodatabaseformat (2.51GB).
  tar_target(
    p1_hybas_gdb,
    "in/BasinATLAS_Data_v10.gdb/BasinATLAS_v10.gdb",
    format = "file"
  ),
  tar_target(
    p1_hybas_legend_xlsx,
    "in/BasinATLAS_Data_v10.gdb/HydroATLAS_v10_Legends.xlsx",
    format = "file"
  )
)

p2 <- list(
  tar_target(
    p2_threats,
    readr::read_csv(p1_threats_csv)
  ),
  tar_target(
    p2_weights,
    readr::read_csv(p1_weights_csv)
  ),
  tar_target(
    p2_hybas_04_sf,
    sf::read_sf(p1_hybas_gdb, layer = "BasinATLAS_v10_lev04")
  ),
  tar_target(
    p2_hybas_legend,
    readxl::read_xlsx(p1_hybas_legend_xlsx, sheet = "fmh_cl")
  ),
  tar_target(
    p2_hybas_habitat_types_sf,
    get_hybas_habitat_types(
      hybas_04_sf = p2_hybas_04_sf,
      hybas_legend = p2_hybas_legend
    )
  ),
  tar_target(
    p2_total_weighted_threats_csv,
    compute_total_weighted_threats(
      threat_data = p2_threats,
      threat_weights = p2_weights,
      hybas_habitat_types = p2_hybas_habitat_types_sf,
      outfile = "../public/findex_total_weighted_threats.csv"
    ),
    format = "file"
  ),
  #### processing for threat maps ####
  tar_target(
    p2_mean_weighted_threats, 
    compute_mean_weighted_threats(
      threat_data = p2_threats, 
      threat_weights = p2_weights, 
      hybas_habitat_types = p2_hybas_habitat_types_sf
    )
  ),
  tar_target( 
    p2_mean_weighted_subThreats, 
    compute_mean_weighted_subThreats(
      threat_data = p2_threats, 
      threat_weights = p2_weights, 
      hybas_habitat_types = p2_hybas_habitat_types_sf
      #sub_threat = subThreat_cat
    )
  ),
  tar_target(
    p2_habitat_subthreats,
    c("Dams", "Wetland drainage", "Deforestation and associated runoff",
      "Riparian degradation", "Agricultural extraction", "Urban extraction", 
      "Industrial extraction")
  ),
  tar_target(
    p2_pollution_subthreats,
    c("Agricultural effluents", "Urban wastewater", 
      "Industrial effluents", "Aquaculture effluents",
      "Pharmaceuticals", "Oil or gas exploration",
      "Plastics", "Mining")
  ),
  tar_target(
    p2_climate_subthreats,
    c("Change in water temperature", "Drought", "Change in flooding", 
      "Change in wind patterns", "Change in ice cover")
  )
)

p3 <- list(
  tar_target(p3_color_pal,
             tibble(
               Habitat_pal = list(c("#7A562B", "#C7985F", "#E1C8AA")), #A07138
               Pollution_pal = list(c("#002D5E", "#B2C0CE")),
               Exploitation_pal = list(c("#B74F49", "#E2B8B6")),
               Invasive_pal = list(c("#4E6D6E", "#C9D8D9")),
               Climate_pal = list(c("#835192", "#DDCCE2")) #9D6AAC
             )),
  #### major threat category maps ####
  tar_map(
    values = tibble::tibble(threat_cat = c("Habitat", "Exploitation", 
                                           "Invasive species", "Pollution", 
                                           "Climate and weather")),
    tar_target( 
      # change all out file paths to the earth in flux directroy: src/assets/images
      p3_threat_map,
      threat_map(in_dat = p2_mean_weighted_threats, 
                 threat_category = threat_cat, 
                 threat_pal = p3_color_pal 
      )
    ),
    tar_target( 
      p3_threat_map_png,
      {
        final_plot <- p3_threat_map + 
          theme(legend.position = "none")
        
        out_file <- paste0("../src/assets/images/", str_replace_all(threat_cat, " ", "_"), "_map.png")
        
        ggsave(out_file, 
               final_plot, height = 6, width = 10, dpi = 300)
      },
      format = "file"
    ),
    tar_target(
      p3_legend_png,
      {
        plot_legend <- get_plot_component(p3_threat_map, "guide-box-right", return_all = T)
        
        out_file <- paste0("../src/assets/images/", str_replace_all(threat_cat, " ", "_"), "_legend.png")
        
        ggsave(out_file, 
               plot_legend, dpi = 300, bg = "transparent")
        knitr::plot_crop(out_file)
      }
    )),
  #### subcategory threat maps ####
  tar_map(
    values = tibble::tibble(sub_threat_cat = c("Dams", "Wetland drainage", "Deforestation and associated runoff",
                                               "Riparian degradation", "Agricultural effluents", "Agricultural extraction",
                                               "Aquaculture effluents", "Industrial effluents", "Industrial extraction",
                                               "Urban extraction", "Urban wastewater", "Overfishing",
                                               "Invasive non-native species", "Change in water temperature", "Drought",
                                               "Change in flooding", "Pharmaceuticals", "Plastics",
                                               "Oil or gas exploration", "Mining", "Change in wind patterns",
                                               "Change in ice cover")),
    tar_target(
      p3_sub_threat_map,
      subThreat_map(in_dat = p2_mean_weighted_subThreats, 
                    threat_category = sub_threat_cat, 
                    threat_pal = p3_color_pal)
    ),
    tar_target(
      p3_sub_threat_map_png,
      {
        final_plot <- p3_sub_threat_map + 
          theme(legend.position = "none")
        
        if(sub_threat_cat %in% p2_habitat_subthreats){
          out_file <- paste0("../src/assets/images/H_", str_replace_all(sub_threat_cat, " ", "_"), "_map.png")
        } else if(sub_threat_cat %in% p2_pollution_subthreats){
          out_file <- paste0("../src/assets/images/P_", str_replace_all(sub_threat_cat, " ", "_"), "_map.png")
        } else if(sub_threat_cat == "Overfishing"){
          out_file <- paste0("../src/assets/images/E_", str_replace_all(sub_threat_cat, " ", "_"), "_map.png")
        } else if(sub_threat_cat == "Invasive non-native species"){
          out_file <- paste0("../src/assets/images/IS_", str_replace_all(sub_threat_cat, " ", "_"), "_map.png")
        } else if(sub_threat_cat %in% p2_climate_subthreats){
          out_file <- paste0("../src/assets/images/CW_", str_replace_all(sub_threat_cat, " ", "_"), "_map.png")
        }
        
        ggsave(out_file, 
               final_plot, height = 6, width = 10, dpi = 300)
      },
      format = "file"
    ),
    tar_target(
      p3_sub_threat_legend_png,
      {
        plot_legend <- get_plot_component(p3_sub_threat_map, "guide-box-right", return_all = T)
        
        if(sub_threat_cat %in% p2_habitat_subthreats){
          out_file <- paste0("../src/assets/images/H_", str_replace_all(sub_threat_cat, " ", "_"), "_legend.png")
        } else if(sub_threat_cat %in% p2_pollution_subthreats){
          out_file <- paste0("../src/assets/images/P_", str_replace_all(sub_threat_cat, " ", "_"), "_legend.png")
        } else if(sub_threat_cat == "Overfishing"){
          out_file <- paste0("../src/assets/images/E_", str_replace_all(sub_threat_cat, " ", "_"), "_legend.png")
        } else if(sub_threat_cat == "Invasive non-native species"){
          out_file <- paste0("../src/assets/images/IS_", str_replace_all(sub_threat_cat, " ", "_"), "_legend.png")
        } else if(sub_threat_cat %in% p2_climate_subthreats){
          out_file <- paste0("../src/assets/images/CW_", str_replace_all(sub_threat_cat, " ", "_"), "_legend.png")
        }
        
        ggsave(out_file, 
               plot_legend, dpi = 300, bg = "transparent")
        knitr::plot_crop(out_file)
      }
    ))
  )


c(p1, p2, p3)
