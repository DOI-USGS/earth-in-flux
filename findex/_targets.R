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
  ),
  tar_target(
    p1_proj,
    "ESRI:54030"
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
  
  ## set up threat and sub-category threat lists for visualization dynamic branching
  
  tar_target(
    p2_threat_categories,
    p2_weights |>
      pull(Threat_Category) |> 
      unique()
  ),
  tar_target(
    p2_threat_subcategories,
    p2_mean_weighted_subThreats |>
      pull(ThreatCategory) |> 
      unique()
  ),
  tar_target(
    p2_habitat_subthreats,
    p2_mean_weighted_subThreats |> 
      filter(MajorCat == "Habitat") |> 
      pull(ThreatCategory) |>
      unique()
  ),
  tar_target(
    p2_pollution_subthreats,
    p2_mean_weighted_subThreats |> 
      filter(MajorCat == "Pollution") |> 
      pull(ThreatCategory) |>
      unique()
  ),
  tar_target(
    p2_climate_subthreats,
    p2_mean_weighted_subThreats |> 
      filter(MajorCat == "Climate and weather") |> 
      pull(ThreatCategory) |>
      unique()
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
    tar_target( 
      p3_threat_map_png,
      {
        final_plot <- threat_map(in_dat = p2_mean_weighted_threats, 
                                 threat_category = p2_threat_categories, 
                                 threat_pal = p3_color_pal,
                                 proj = p1_proj) + 
          theme(legend.position = "none")
        
        out_file <- paste0("../src/assets/images/", str_replace_all(p2_threat_categories, " ", "_"), "_map.png")
        
        ggsave(out_file, 
               final_plot, height = 6, width = 10, dpi = 300)
      },
      format = "file",
      pattern = p2_threat_categories
    ),
    tar_target(
      p3_legend_png,
      {
        final_plot <- threat_map(in_dat = p2_mean_weighted_threats, 
                                 threat_category = p2_threat_categories, 
                                 threat_pal = p3_color_pal,
                                 proj = p1_proj)
        
        plot_legend <- get_plot_component(final_plot, "guide-box-right", return_all = T)
        
        out_file <- paste0("../src/assets/images/", str_replace_all(p2_threat_categories, " ", "_"), "_legend.png")
        
        ggsave(out_file, 
               plot_legend, dpi = 300, bg = "transparent")
        knitr::plot_crop(out_file)
      },
      pattern = p2_threat_categories
    ),
    tar_target( # will turn this into a function
      p3_sub_threat_map_png,
      {
        final_plot <- subThreat_map(in_dat = p2_mean_weighted_subThreats, 
                      threat_category = p2_threat_subcategories, 
                      subcat_habitat = p2_habitat_subthreats,
                      subcat_pollution = p2_pollution_subthreats,
                      subcat_climate = p2_climate_subthreats,
                      threat_pal = p3_color_pal,
                      proj = p1_proj) + 
          theme(legend.position = "none")
        
        if(p2_threat_subcategories %in% p2_habitat_subthreats){
          out_file <- paste0("../src/assets/images/H_", str_replace_all(p2_threat_subcategories, " ", "_"), "_map.png")
        } else if(p2_threat_subcategories %in% p2_pollution_subthreats){
          out_file <- paste0("../src/assets/images/P_", str_replace_all(p2_threat_subcategories, " ", "_"), "_map.png")
        } else if(p2_threat_subcategories == "Overfishing"){
          out_file <- paste0("../src/assets/images/E_", str_replace_all(p2_threat_subcategories, " ", "_"), "_map.png")
        } else if(p2_threat_subcategories == "Invasive non-native species"){
          out_file <- paste0("../src/assets/images/IS_", str_replace_all(p2_threat_subcategories, " ", "_"), "_map.png")
        } else if(p2_threat_subcategories %in% p2_climate_subthreats){
          out_file <- paste0("../src/assets/images/CW_", str_replace_all(p2_threat_subcategories, " ", "_"), "_map.png")
        }
        
        ggsave(out_file, 
               final_plot, height = 6, width = 10, dpi = 300)
      },
      format = "file",
      pattern = p2_threat_subcategories
    ),
    tar_target( # will turn this into a function
      p3_sub_threat_legend_png,
      {
        final_plot <- subThreat_map(in_dat = p2_mean_weighted_subThreats, 
                                    threat_category = p2_threat_subcategories, 
                                    subcat_habitat = p2_habitat_subthreats,
                                    subcat_pollution = p2_pollution_subthreats,
                                    subcat_climate = p2_climate_subthreats,
                                    threat_pal = p3_color_pal,
                                    proj = p1_proj)
        
        plot_legend <- get_plot_component(final_plot, "guide-box-right", return_all = T)
        
        if(p2_threat_subcategories %in% p2_habitat_subthreats){
          out_file <- paste0("out/H_", str_replace_all(p2_threat_subcategories, " ", "_"), "_legend_raw.png")
        } else if(p2_threat_subcategories %in% p2_pollution_subthreats){
          out_file <- paste0("out/P_", str_replace_all(p2_threat_subcategories, " ", "_"), "_legend_raw.png")
        } else if(p2_threat_subcategories == "Overfishing"){
          out_file <- paste0("out/E_", str_replace_all(p2_threat_subcategories, " ", "_"), "_legend_raw.png")
        } else if(p2_threat_subcategories == "Invasive non-native species"){
          out_file <- paste0("out/IS_", str_replace_all(p2_threat_subcategories, " ", "_"), "_legend_raw.png")
        } else if(p2_threat_subcategories %in% p2_climate_subthreats){
          out_file <- paste0("out/CW_", str_replace_all(p2_threat_subcategories, " ", "_"), "_legend_raw.png")
        }
        
        ggsave(out_file, 
               plot_legend, dpi = 300, bg = "transparent")
        knitr::plot_crop(out_file)
        
        if(p2_threat_subcategories %in% p2_habitat_subthreats){
          out_file_final <- paste0("../src/assets/images/H_", str_replace_all(p2_threat_subcategories, " ", "_"), "_legend.png")
        } else if(p2_threat_subcategories %in% p2_pollution_subthreats){
          out_file_final <- paste0("../src/assets/images/P_", str_replace_all(p2_threat_subcategories, " ", "_"), "_legend.png")
        } else if(p2_threat_subcategories == "Overfishing"){
          out_file_final <- paste0("../src/assets/images/E_", str_replace_all(p2_threat_subcategories, " ", "_"), "_legend.png")
        } else if(p2_threat_subcategories == "Invasive non-native species"){
          out_file_final <- paste0("../src/assets/images/IS_", str_replace_all(p2_threat_subcategories, " ", "_"), "_legend.png")
        } else if(p2_threat_subcategories %in% p2_climate_subthreats){
          out_file_final <- paste0("../src/assets/images/CW_", str_replace_all(p2_threat_subcategories, " ", "_"), "_legend.png")
        }
        
        cowplot_legend(in_dat = p2_mean_weighted_subThreats, legend_png = out_file, threat_category = p2_threat_subcategories, out_file = out_file_final)
      },
      format = "file",
      pattern = p2_threat_subcategories
    )
  )


c(p1, p2, p3)
