clean_input_data <- function(data_file) {
  data <- read_xlsx(data_file, na = c("<MDL", '--'))
  names(data) <- make.names(names(data),unique = TRUE)
  
  data <- data |>
    replace_na(list(Levoglucosan.in.picogram.per.milliliter..unfiltered. = 0, 
                    Mannosan.in.picogram.per.milliliter..unfiltered. = 0,
                    Galactosan.in.picogram.per.milliliter..unfiltered. = 0)) |>
    mutate(L_M = Levoglucosan.in.picogram.per.milliliter..unfiltered./
             Mannosan.in.picogram.per.milliliter..unfiltered.,
           L_MG = Levoglucosan.in.picogram.per.milliliter..unfiltered./
             (Mannosan.in.picogram.per.milliliter..unfiltered.  + 
                Galactosan.in.picogram.per.milliliter..unfiltered.),
           year = lubridate::year(Sample.date),
           L_MG_softwood = (L_MG >= 0.4 & L_MG <= 6.1),
           L_MG_hardwood = (L_MG >= 1.5 & L_MG <= 17.6)) 

}

extract_data_for_site <- function(data, year, site) {
  munged_data <- data |>
    select(id = Ice.Core.ID, depth_cm = Top.depth.of.ice.core.sample.in.centimeters,
           everything()) |>
    mutate(year = as.integer(stringr::str_sub(id, 1, 4)),
           site = stringr::str_sub(id, -1, -1), .before = 1) 
  
  requested_id <- sprintf('%s Core %s', year, site)
  
  if (!requested_id %in% munged_data[['id']]) {
    stop(message(sprintf('There are no data for site %s in %s', site, year)))
  }
  
  munged_data |>
    filter(year == !!year, site == !!site)
}

export_particulates_data <- function(data, outfile) {
  data |>
    select(year, site, depth_cm, contains('particles')) |>
    pivot_longer(cols = contains('particles'), names_to = "particle_size_microns",
                 values_to = 'number_per_mL', names_pattern = "(\\d+)",
                 names_transform = list(particle_size_microns = as.integer)) |>
    group_by(depth_cm) |>
    summarize(total = sum(number_per_mL, na.rm = TRUE)) |>
    readr::write_csv(outfile)
  
  return(outfile)
}

process_sugars_data <- function(data) {
  data |>
    select(year, site, depth_cm, contains('picogram')) |>
    pivot_longer(cols = contains('picogram'), names_to = "sugar",
                 values_to = 'picogram_per_mL', names_pattern = "([A-Za-z]+).*")  |>
    pivot_wider(id_cols = c(depth_cm, year, site),
                names_from = sugar, values_from = picogram_per_mL) |>
    # filter out high Mannosan value, per paper (https://iopscience.iop.org/article/10.1088/1748-9326/ab8fd2/pdf)
    mutate(Mannosan = ifelse(Mannosan < 600, Mannosan, 0)) |>
    pivot_longer(cols = -c(depth_cm, year, site), names_to = "sugar",
                 values_to = 'picogram_per_mL') |>
    mutate(bar_order = case_when(
      sugar == 'Mannosan' ~ 1,
      sugar == 'Galactosan' ~ 2,
      TRUE ~ 3
    ))
}

export_sugars_data <- function(data, outfile) {
  readr::write_csv(data, outfile)
  
  return(outfile)
}

export_biomass_data <- function(sugars_data, exclude_grass, outfile) {

  biomass_data <- sugars_data |>
    pivot_wider(id_cols = c(depth_cm, year, site),
                names_from = sugar, values_from = picogram_per_mL) |>
    mutate(L_MG = Levoglucosan / (Mannosan + Galactosan),
           softwood = (L_MG >= 0.4 & L_MG <= 6.1),
           hardwood = (L_MG >= 1.5 & L_MG <= 17.6),
           grass = (L_MG >= 1.7 & L_MG <= 2.5)) |>
    select(-L_MG, -Levoglucosan, -Mannosan, -Galactosan) |>
    pivot_longer(cols = c(softwood, hardwood, grass), names_to = "vegetation_type",
                 values_to = 'burned') |>
    filter(burned) |>
    mutate(vegetation_type = factor(vegetation_type, 
                                    levels = c('grass', 'softwood', 'hardwood'))) |>
    arrange(vegetation_type)
  
  if (exclude_grass) {
    biomass_data <- filter(biomass_data, !vegetation_type == 'grass')
  }
  
  readr::write_csv(biomass_data, outfile)
  
  return(outfile)
}
