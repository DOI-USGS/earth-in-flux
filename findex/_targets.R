library(targets)
library(tarchetypes)

options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("tidyverse", 
                            "sf",
                            "readxl"))

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
    p2_mean_weighted_threats, # take _csv ending off if not saving to a csv
    compute_mean_weighted_threats(
      threat_data = p2_threats, 
      threat_weights = p2_weights, 
      hybas_habitat_types = p2_hybas_habitat_types_sf#,
      #outfile = "out/findex_mean_weighted_threats.csv"
    )#,
    #format = "file"
  )
)

p3 <- list(
  tar_target(p3_color_pal,
             tibble(
               Habitat_pal = list(c("#A07138", "#C08B4B", "#CDA371", "#E1C8AA")),
               Pollution_pal = list(c("#002D5E", "#B2C0CE")),
               Exploitation_pal = list(c("#B74F49", "#cd8480", "#E2B8B6")),
               Invasive_pal = list(c("#4E6D6E", "#6f9899", "#9CB8B9", "#C9D8D9")),
               Climate_pal = list(c("#9D6AAC", "#bd9bc7", "#DDCCE2"))
             )),
  tar_map(
    values = tibble::tibble(threat_cat = c("Habitat", "Exploitation", 
                                           "Invasive species", "Pollution", 
                                           "Climate and weather")),
  tar_target(
    p3_threat_map,
    threat_map(in_dat = p2_mean_weighted_threats, 
                   threat_category = threat_cat, 
                   threat_pal = p3_color_pal, 
                   out_file = paste0("out/", str_replace_all(threat_cat, " ", "_"), "_map.png")),
    format = "file"
  ))
)

c(p1, p2, p3)
