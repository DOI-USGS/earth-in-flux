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
  
  #### processing for threat maps ####
  
  tar_target(
    p2_weighted_threats,
    compute_weighted_threats(
      threat_data = p2_threats,
      threat_weights = p2_weights,
      hybas_habitat_types = p2_hybas_habitat_types_sf
    )
  ),
  tar_target(
    p2_total_weighted_threats_csv,
    compute_total_weighted_threats(
      in_dat = p2_weighted_threats,
      outfile = "../public/findex_total_weighted_threats.csv"
    ),
    format = "file"
  ),
  tar_target(
    p2_mean_weighted_threats, 
    compute_mean_weighted_threats(
      in_dat = p2_weighted_threats
    )
  ),
  tar_target( 
    p2_mean_weighted_subThreats, 
    compute_mean_weighted_subThreats(
      in_dat = p2_weighted_threats
    )
  ),
  
  #### threat lists for branching ####
  
  tar_target(
    p2_threat_categories,
    p2_mean_weighted_threats |>
      pull(ThreatCategory) |> 
      unique()
  ),
  tar_target(
    p2_threat_subcategories,
    p2_mean_weighted_subThreats |>
      pull(ThreatCategory) |> 
      unique()
  ),
  
  #### color ramps and file name templates ####
  
  tar_target(
    p2_viz_config,
    {
      p2_mean_weighted_subThreats |>
        select(MajorCat, ThreatCategory) |>
        mutate(
          # color ramps
          pal = case_when(
            MajorCat == "Habitat" ~ list(c("#4E6D6E", "#C9D8D9")),#), 
            MajorCat == "Fishing pressure" ~ list(c("#835192", "#DDCCE2")), 
            MajorCat == "Invasive species" ~ list(c("#B74F49", "#E2B8B6")), 
            MajorCat == "Pollution" ~ list(c("#7A562B", "#C7985F", "#E1C8AA")), 
            MajorCat == "Climate and weather" ~ list(c("#002D5E", "#B2C0CE"))
          ),
          # file name templates
          threat_legend_raw = "out/%s_legend_raw.png",
          threat_legend = "../src/assets/images/%s_legend.png",
          threat_map = "../src/assets/images/%s_map.png",
          subThreat_legend_raw = case_when(
            MajorCat %in% "Habitat" ~ "out/H_%s_legend_raw.png",
            MajorCat %in% "Fishing pressure" ~ "out/FP_%s_legend_raw.png",
            MajorCat %in% "Invasive species" ~ "out/IS_%s_legend_raw.png",
            MajorCat %in% "Pollution" ~ "out/P_%s_legend_raw.png",
            MajorCat %in% "Climate and weather" ~ "out/CW_%s_legend_raw.png"
          ),
          subThreat_legend = case_when(
            MajorCat %in% "Habitat" ~ "../src/assets/images/H_%s_legend.png",
            MajorCat %in% "Fishing pressure" ~ "../src/assets/images/FP_%s_legend.png",
            MajorCat %in% "Invasive species" ~ "../src/assets/images/IS_%s_legend.png",
            MajorCat %in% "Pollution" ~ "../src/assets/images/P_%s_legend.png",
            MajorCat %in% "Climate and weather" ~ "../src/assets/images/CW_%s_legend.png"
          ),
          subThreat_map = case_when(
            MajorCat %in% "Habitat" ~ "../src/assets/images/H_%s_map.png",
            MajorCat %in% "Fishing pressure" ~ "../src/assets/images/FP_%s_map.png",
            MajorCat %in% "Invasive species" ~ "../src/assets/images/IS_%s_map.png",
            MajorCat %in% "Pollution" ~ "../src/assets/images/P_%s_map.png",
            MajorCat %in% "Climate and weather" ~ "../src/assets/images/CW_%s_map.png"
          )
        )
    }
  )
)

p3 <- list(
  
  #### major threat maps and legends ####
  
  tar_target( 
    p3_threat_map_png,
    {
      final_plot <- threat_map(in_dat = p2_mean_weighted_threats, 
                               threat_category = p2_threat_categories, 
                               threat_pal = p2_viz_config,
                               hybas_habitat_types = p2_hybas_habitat_types_sf,
                               proj = p1_proj)  + 
        theme(legend.position = "none")
      
      save_map(type = "threat", plot = final_plot, 
               threat_category = p2_threat_categories, 
               threat_pal = p2_viz_config,
               height = 6, width = 10, dpi = 300)
    },
    format = "file",
    pattern = p2_threat_categories
  ),
  tar_target(
    p3_legend_png,
    {
      final_plot <- threat_map(in_dat = p2_mean_weighted_threats, 
                               threat_category = p2_threat_categories, 
                               threat_pal = p2_viz_config,
                               proj = p1_proj,
                               hybas_habitat_types = p2_hybas_habitat_types_sf)+
        theme(legend.text = element_blank())
      
      save_legend(type = "threat", plot = final_plot, 
                  threat_category = p2_threat_categories, 
                  in_dat = p2_mean_weighted_threats,
                  threat_pal = p2_viz_config, 
                  height = 176, width = 429, unit = "px", dpi = 300)
    },
    pattern = p2_threat_categories
  ),
  
  #### sub-category treat maps and legends ####
  
  tar_target( 
    p3_sub_threat_map_png,
    {
      final_plot <- subThreat_map(in_dat = p2_mean_weighted_subThreats, 
                                  threat_category = p2_threat_subcategories, 
                                  threat_pal = p2_viz_config,
                                  proj = p1_proj,
                                  hybas_habitat_types = p2_hybas_habitat_types_sf) + 
        theme(legend.position = "none")
      
      save_map(type = "subThreat", plot = final_plot, 
               threat_category = p2_threat_subcategories, 
               threat_pal = p2_viz_config,
               height = 6, width = 10, dpi = 300)
    },
    format = "file",
    pattern = p2_threat_subcategories
  ),
  tar_target( 
    p3_sub_threat_legend_png,
    {
      final_plot <- subThreat_map(in_dat = p2_mean_weighted_subThreats, 
                                  threat_category = p2_threat_subcategories, 
                                  threat_pal = p2_viz_config,
                                  proj = p1_proj,
                                  hybas_habitat_types = p2_hybas_habitat_types_sf)+
        theme(legend.text = element_blank())
      
      save_legend(type = "subThreat", plot = final_plot, 
                  threat_category = p2_threat_subcategories, 
                  in_dat = p2_mean_weighted_subThreats,
                  threat_pal = p2_viz_config,
                  height = 176, width = 429, unit = "px", dpi = 300)
    },
    format = "file",
    pattern = p2_threat_subcategories
  ),
  # top threat in each basin globally
  tar_target(
    p3_top_threat_map_png,
    {
      final_plot <- top_threat_plot(in_dat = p2_mean_weighted_threats, 
                                    threat_pal = p2_viz_config, 
                                    hybas_habitat_types = p2_hybas_habitat_types_sf, 
                                    proj = p1_proj,
                                    threat_category = p2_threat_categories)  + 
        theme(legend.position = "none")
      
      # change to actual directory once design is finalized --------------------
      ggsave(sprintf("../src/assets/images/%s_threat_by_basin.png", str_replace_all(p2_threat_categories, " ", "_")), 
             final_plot, height = 6, width = 10, dpi = 300)
      
      # change to actual directory once design is finalized --------------------
      knitr::plot_crop(sprintf("../src/assets/images/%s_threat_by_basin.png", str_replace_all(p2_threat_categories, " ", "_")))
    },
    format = "file",
    pattern = p2_threat_categories
  ),
  tar_target(
    p3_all_top_threat_map_png,
    {
      final_plot <- top_threat_plot(in_dat = p2_mean_weighted_threats, 
                                    threat_pal = p2_viz_config, 
                                    hybas_habitat_types = p2_hybas_habitat_types_sf, 
                                    proj = p1_proj,
                                    threat_category = "none")  + 
        theme(legend.position = "none")
      
      # change to actual directory once design is finalized --------------------
      ggsave("../src/assets/images/all_threat_by_basin.png", 
             final_plot, height = 6, width = 10, dpi = 300)
      
      # change to actual directory once design is finalized --------------------
      knitr::plot_crop("../src/assets/images/all_threat_by_basin.png")
    },
    format = "file"
  ),
  tar_target(
    p3_top_threat_thumbnail,
    top_threat_thumbnail(in_dat = p2_mean_weighted_threats, 
                         threat_pal = p2_viz_config, 
                         hybas_habitat_types = p2_hybas_habitat_types_sf, 
                         proj = p1_proj,
                         threat_category = "none",
                         height = 6,
                         width = 6,
                         dpi = 300,
                         out_file = "../src/assets/images/threat_by_basin_thumbnail.png"),
    format = "file"
  )
#  tar_target(
#    p3_top_threat_legend_png,
#    {
#      final_plot <- top_threat_plot(in_dat = p2_mean_weighted_threats, 
#                                    threat_pal = p2_viz_config, 
#                                    hybas_habitat_types = p2_hybas_habitat_types_sf, 
#                                    proj = p1_proj,
#                                    threat_category = "none")
#      
#      save_top_threat_legend(plot = final_plot, 
#                             dpi = 300, 
#                             # change to actual directory once design is finalized
#                             out_file = "../src/assets/images/threat_by_basin_legend.png")
#    },
#    format = "file"
#  )
)


c(p1, p2, p3)
