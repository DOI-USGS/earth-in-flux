clean_input_data <- function(data_file) {
  data <- read_xlsx(data_file, na="NA")
  
  names(data)<-make.names(names(data),unique = TRUE)
  
  data |>
    mutate(
      family = ifelse(family == 'Salmonindae', 'Salmonidae', family),
      species_common = case_when(
        species_common == 'Bey?ehir bleak' ~ 'Beyşehir bleak',
        species_taxa_sci == 'Brachyplatystoma rousseauxii' ~ 'Gilded catfish',
        species_common == 'Speckled pavon, speckled peacock bass (percid)' ~ 
          'Speckled pavon',
        TRUE ~ species_common
      )
    )
}

build_nested_json <- function(data, focal_columns, out_file) {
  data_list <- data |>
    filter(!is.na(focal_columns[['value']])) |>
    group_by(name = family) |>
    group_modify(~ {
      species_list <- .x |> 
        group_by(name = species_common) |>
        group_modify(~ {
          country_list <- .x |>
            select(!!focal_columns) |>
            split(nrow(.x))
          return(tibble('children' = country_list))
        })
      return(tibble('children' = list(species_list)))
    })
  
  json_data <- jsonlite::toJSON(data_list, auto_unbox = TRUE, pretty = TRUE)
  json_data <- paste0('{\n "name": "fish", "children":', json_data, '\n}')
  write(json_data, file = out_file)
  
  return(out_file)
}

build_country_climate_summary <- function(data) {
  # get country-level summary of harvest, consumption, and climate vulnerability
  data_country <- data |>
    filter(total_rec_harvest_kg > 0) |>
    group_by(admin) |>
    summarize(n_fishers = first(n_fishers),
              total_consumable_harv_kg = sum(total_consumable_harv_kg),
              MCDM_VUL_2075_45 = mean(MCDM_VUL_2075_45, na.rm = TRUE)) |>
    mutate(consum_kg_fisher = total_consumable_harv_kg / n_fishers)
  
  # pull in natural earth data to get region and economic classifications
  rnaturalearth::ne_countries(returnclass = 'sf') %>%
    dplyr::select(admin, continent) |>
    sf::st_drop_geometry() |>
    dplyr::right_join(data_country)
}

build_harvest_csv <- function(data, out_file, n_top_families) {
  data_subset <- data |>
    filter(!is.na(species_common)) |>
    filter(uncertainty_classification < 4)
  
  # Determine harvest by family and order most to least
  family_summary <- data_subset |> 
    group_by(family) |>
    summarise(value = sum(total_biomass_harv_kg, na.rm = TRUE)) |>
    arrange(desc(value))
  
  # select subset of highest harvest families
  focal_families <- family_summary |>
    head(n_top_families) |>
    pull(family)
  
  # filter to subset of families
  data_subset <- data_subset |>
    filter(family %in% focal_families)
  
  # aggregate harvest to family level
  family_species <- 
    data_subset |> 
    group_by(family, species_common) |>
    summarise(value = sum(total_biomass_harv_kg, na.rm = TRUE)) |>
    mutate(percent_harvest = value/sum(value) * 100) |>
    # filter to species representing > 1% of the harvest within each family
    filter(percent_harvest > 1) |>
    select(source = family, target = species_common, value)
  
  # Pull country specific harvest for species including in family level aggregation
  species_country <- data_subset |> 
    filter(species_common %in% unique(pull(family_species, target))) |>
    select(source = species_common, target = admin, value = total_biomass_harv_kg)
  
  bind_rows(family_species, species_country) |>
    readr::write_csv(out_file)
  
  return(out_file)
}