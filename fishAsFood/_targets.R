library(targets)

options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("tidyverse", 
                            "readxl",
                            "jsonlite",
                            "rnaturalearth"))

source('src/data_utils.R')
source('src/plot_utils.R')

p0 <- list(
  tar_target(
    p0_sb_id,
    '644ae0e0d34e45f6ddccf773'
  )
)

p1 <- list(
  tar_target(
    p1_data_xlsx,
    {
      data_file <- 'Rec fish food_20230509_for USGS data release.xlsx'
      out_file <- file.path('in', data_file)
      sbtools::item_file_download(sb_id = p0_sb_id, 
                                  names = data_file, 
                                  destinations = out_file)
      return(out_file)
    },
    format = 'file'
  ),
  tar_target(
    p2_metadata_xlsx,
    file.path('in', 'recfishfoods_families_thermal.xlsx'),
    format = 'file'
  )
)

p2 <- list(
  tar_target(
    p2_data,
    clean_input_data(data_file = p1_data_xlsx)
  ),
  tar_target(
    p2_price_json,
    build_nested_json(
      data = p2_data, 
      focal_columns = c('name' = 'admin', 'value' = 'total_value_species'), 
      out_file = '../public/total_price.json'
    ),
    format = 'file'
  ),
  tar_target(
    p2_country_climate_summary,
    build_country_climate_summary(
      data = p2_data,
      metadata_xlsx = p2_metadata_xlsx
    )
  ),
  tar_target(
    p2_country_climate_csv,
    {
      outfile <- '../public/fish_as_food_country_climate.csv'
      readr::write_csv(p2_country_climate_summary, outfile)
      return(outfile)
    },
    format = 'file'
  ),
  tar_target(
    p2_harvest_csv,
    build_harvest_csv(
      data = p2_data,
      n_top_families = 5,
      out_file = '../public/fish_as_food_harvest.csv'
    ),
    format = 'file'
  )
)

p3 <- list(
  tar_target(
    p3_continent_map_png,
    generate_continent_map(
      country_data = p2_country_climate_summary,
      out_file = '../src/assets/images/fish_as_food_continent_map.png'
    )
  )
)

c(p0, p1, p2, p3)
