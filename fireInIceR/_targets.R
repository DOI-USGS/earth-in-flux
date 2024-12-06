library(targets)

options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("tidyverse", 
                            "readxl", "readr"))

source('src/data_utils.R')

p0 <- list(
  tar_target(
    p0_sb_id,
    '5e472c3ee4b0ff554f6837bc'
  )
)

p1 <- list(
  tar_target(
    p1_data_xlsx,
    {
      data_file <- 'Table 4 Constituent Data.xlsx'
      out_file <- file.path('in', data_file)
      sbtools::item_file_download(sb_id = p0_sb_id, 
                                  names = data_file, 
                                  destinations = out_file)
      return(out_file)
    },
    format = 'file'
  )
)

p2 <- list(
  tar_target(
    p2_data,
    clean_input_data(data_file = p1_data_xlsx)
  ),
  tar_target(
    p2_focal_core,
    list(year = 2016,
         site = 4)
  ),
  tar_target(
    p2_data_core,
    extract_data_for_site(
      data = p2_data,
      year = p2_focal_core[['year']],
      site = p2_focal_core[['site']]
    )
  ),
  tar_target(
    p2_core_particulates_csv,
    export_particulates_data(
      data = p2_data_core,
      outfile = 
        sprintf("../public/fii_core%sparticulates.csv", p2_focal_core[['site']])
    ),
    format = 'file'
  ),
  tar_target(
    p2_core_sugars,
    process_sugars_data(
      data = p2_data_core
    )
  ),
  tar_target(
    p2_core_sugars_csv,
    export_sugars_data(
      data = p2_core_sugars,
      outfile = 
        sprintf("../public/fii_core%ssugars.csv", p2_focal_core[['site']])
    ),
    format = 'file'
  ),
  tar_target(
    p2_core_biomass_csv,
    export_biomass_data(
      sugars_data = p2_core_sugars,
      exclude_grass = FALSE,
      outfile = 
        sprintf("../public/fii_core%sbiomass.csv", p2_focal_core[['site']])
    ),
    format = 'file'
  )
)

c(p0, p1, p2)
