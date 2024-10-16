library(targets)

options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("tidyverse", 
                            "sf",
                            "readxl"))

source("src/data_utils.R")

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
  )
)

c(p1, p2)
