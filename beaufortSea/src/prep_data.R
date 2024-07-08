join_abundance <- function(ostracode_in, foram_in, color_long){
  ## Join data
  ostracodes <- ostracode_in |>
    rename(year = `Calendar AGE`) |>
    mutate(year = round(year, 0))
  
  forams <- foram_in |>
    rename(year = 'calendar yr') |>
    mutate(year = round(year, 0))
  
  join_long <- ostracodes |>
    full_join(forams, by = "year") |>
    pivot_longer(-year, names_to = "species") |>
    mutate(type = case_when(species %in% names(ostracodes) ~ "Ostracode",
                            species %in% names(forams) ~ "Foram",
                            TRUE ~ NA)) |> 
    # add genus designation
    mutate(genus = word(species)) |>
    mutate(genus = case_when(genus == "E." ~ "Elphidium",
                             genus == "C." ~ "Cassidulina",
                             genus == "Indeterminate/IRO..." ~ "Other",
                             genus == "Other...122" ~ "Other",
                             genus == "Dentalina...39" ~ "Dentalina",
                             genus == "Polymorphinids...41" ~ "Polymorphinids",
                             genus == "Spiroplectammina_spp" ~ "Spiroplectammina",
                             TRUE ~ genus)) |>
    # designate names
    mutate(name = case_when(genus == "Spiroplectammina" ~ "Spiroplectammina",
                            genus == "Cassidulina" ~ "Cassidulina",
                            genus == "Elphidium" ~ "Elphidium",
                            genus == "Kotoracythere" ~ "Kotoracythere",
                            genus == "Paracyprideis" ~ "Paracyprideis",
                            type == "Ostracode" ~ "Ostracode",
                            type == "Foram" ~ "Foram",
                            TRUE ~ NA)) 
  
  join_to_colors <- join_long |>
    left_join(color_long, by = "name")
  
  return(join_to_colors)
}


decade_abundance <- function(data_in){
  # Not all years match exactly, so going average by every 50 years and species
  join_decades <- data_in |>
    mutate(decade = year - year %% 100) |>
    group_by(decade, species, type, genus, name, hexcode) |>
    summarize(mean_abundance = mean(value, na.rm = TRUE))
  
  # Not all 50-year periods have both species types, so need to add up total abundance by
  # period to get actual percent abundance 
  total_abundance_by_decade <- join_decades |>
    group_by(decade) |>
    summarize(total_abundance = sum(mean_abundance, na.rm = TRUE)) 
  
  # Put both values together
  pct_abundance <- join_decades |>
    left_join(total_abundance_by_decade, by = "decade") |>
    mutate(pct_abundance = (mean_abundance/total_abundance)*100)
  
  return(pct_abundance)
}