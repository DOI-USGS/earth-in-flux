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

build_climate_csv <- function(data, metadata_file, out_file) {
  # read in metadata
  metadata <- read_xlsx(metadata_file)
  
  # filter out missing and uncertain data
  data_subset <- data |>
    filter(!is.na(species_common), !is.na(MCDM_VUL_2030_45)) |>
    left_join(metadata) |>
    filter(uncertainty_classification < 4)
  
  # get subset of better represented families with >2 species
  focal_families <- data_subset |>
    group_by(family, species_common) |>
    count() |>
    group_by(family) |>
    count() |>
    filter(n > 2) |>
    pull(family)
  
  # filter to subset of families
  data_subset <- data_subset |>
    filter(family %in% focal_families)
  
  # aggregate to species level
  data_species <- data_subset |>
    group_by(family, species_common, thermal_guild, popular_y) |>
    summarize(across(starts_with('MCDM'), 
                     ~ first(.x)),
              across(starts_with('weighted_MCDM'), 
                     ~ mean(.x, na.rm = TRUE)))
  
  # sort aggregated data
  subset_species <- pull(data_species, species_common) |> sort(decreasing = TRUE)
  data_species <- mutate(data_species, species_common = 
                           factor(species_common, levels = subset_species))
  
  # aggregate to family level
  data_family <- data_subset |>
    group_by(family, thermal_guild, popular_y) |>
    summarize(across(starts_with('MCDM'), 
                     ~ first(.x)),
              across(starts_with('weighted_MCDM'), 
                     ~ mean(.x, na.rm = TRUE)))
  
  data_species |>
    select(family, thermal_guild, species = species_common, 
           cvi_2030 = weighted_MCDM_VUL_2030_45, 
           cvi_2075 = weighted_MCDM_VUL_2075_45) |>
    left_join(select(data_family, family, thermal_guild, 
                     cvi_2030_family = weighted_MCDM_VUL_2030_45,
                     cvi_2075_family = weighted_MCDM_VUL_2075_45))|>
    readr::write_csv(out_file)
  
  return(out_file)
}

build_country_climate_csv <- function(data, out_file) {
  # get country-level summary of harvest, consumption, and climate vulnerability
  data_country <- data |>
    filter(total_rec_harvest_kg > 0) |>
    group_by(admin) |>
    summarize(population = first(population),
              n_fishers = first(n_fishers),
              total_rec_harvest_kg = first(total_rec_harvest_kg),
              total_consumable_harv_kg = sum(total_consumable_harv_kg),
              MCDM_VUL_2075_45 = mean(MCDM_VUL_2075_45, na.rm = TRUE)) |>
    mutate(consum_kg_person = total_consumable_harv_kg / population,
           consum_kg_fisher = total_consumable_harv_kg / n_fishers)
  
  # pull in natural earth data to get region and economic classifications
  rnaturalearth::ne_countries(returnclass = 'sf') %>%
    dplyr::select(admin, continent, region_un, subregion, region_wb, economy, 
                  income_grp) |>
    sf::st_drop_geometry() |>
    dplyr::right_join(data_country) |>
    readr::write_csv(out_file)
  
  return(out_file)
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