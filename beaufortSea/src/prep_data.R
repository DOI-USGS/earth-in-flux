# Function to clean up the three foram datasets and merge them
merge_foram_data <- function(GGC30_in, JPC32_in, MC29_in, age_data){
  # remove empty columns
  GGC30_temp <- GGC30_in |> dplyr::select(!contains("...")) |>
    mutate(core = "GGC30")
  
  JPC32_temp <- JPC32_in |> dplyr::select(!contains("...")) |>
    mutate(core = "JPC 32")
  
  # fix typo in name
  MC29_temp <- MC29_in |>
    rename(`E. bartleti` = `E. bartletti`) |>
    mutate(core = "MC-29A")
  
  # join rows for the three cores
  join <- bind_rows(GGC30_temp, JPC32_temp, MC29_temp) |>
    # rename depths
    rename(top = `cm (top)`,
           bottom = `cm (base)`) |>
    # remove unnecessary columns 
    select(-`wet weight`, -`dry washed weight`) |>
    # create identifier of core and depth
    mutate(id = paste0(core, "-", top))
  
  # Create identifier of core and depth and rename vars
  age_data <- age_data |>
    rename(top = `cm (top)`,
           bottom = `cm (base)`,
           year = `calendar yr`) |>
    mutate(id = paste0(core, "-", top)) |>
    select(year, id)
  
  # Merge age data with raw counts
  merged <- join |> left_join(age_data, by = "id")
  
  return(merged)
}

# clean foram data (from raw counts to % abundance)
clean_foram_data <- function(raw_in){
  # fix taxonomic name from I. teretis to Cassidulina teretis re: "T1 Taxonomic Note" in
  # the sciencebase data release table
  renamed <- raw_in |>
    rename(`Cassidulina teretis` = `I. teretis`) |>
    # Plus, spell out all abbreviated genera
    rename(`Cassidulina reniforme` = `C. reniforme`,
           `Elphidium incertum`  = `E. incertum`,
           `Elphidium excavatum clavatum`  = `E. excavatum clavatum`,
           `Elphidium bartletti` = `E. bartleti`) 
  
  # prep to calculate % abundance
  counts <- renamed |>
    select(-top, -bottom, -core, -id) |>
    # NAs should be zeros
    mutate(across(!year, ~replace_na(.x, 0))) |>
    # calculate total count 
    mutate(sum = rowSums(across(!year)))
  
  # calculate percent abundance
  percents <- counts |>
    mutate(across(!year & !sum, ~ (.x/sum)*100))
  
  # take out -6 year
  out_df <- percents |>
    filter(year >= 0)
  
  return(out_df)
  
}


join_abundance <- function(ostracode_in, foram_in, color_long, focal_species){
  ## Join data
  ostracodes <- ostracode_in |>
    rename(year = `Calendar AGE`) |>
    mutate(year = round(year, 0))
  
  forams <- foram_in |>
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
    left_join(color_long, by = "name") |>
    left_join(focal_species, by = "species")
  
  return(join_to_colors)
}


decade_abundance <- function(data_in){
  # Not all years match exactly, so going average by every 50 years and species
  join_decades <- data_in |>
    mutate(decade = year - year %% 100) |>
    group_by(decade, species, type, genus, name, hexcode, 
             species_name, focal_L, epithet, image_name) |>
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

generate_beeswarm_data <- function(data, outfile) {

  # filter data
  data_subset <- data |>
    ungroup() |>
    filter(!is.na(pct_abundance), pct_abundance > 0)
  
  # generate unique ID for each species
  unique_species <- unique(data_subset$species)
  join_df <- tibble(
    species = unique_species
  ) |>
    mutate(species_id = stringr::str_glue('species_{row_number()}'))
  
  # export for use in site
  data_subset <- data_subset |> 
    select(decade, species, type, genus, name, hexcode, species_name, pct_abundance) |>
    group_by(species) |>
    complete(nesting(type, genus, name, hexcode, species_name), 
             decade = seq(0, 2000, 100), 
             fill = list(pct_abundance = 0)) |>
    left_join(join_df) |>
    arrange(decade) 
  
  # add order for bar chart
  data_out <- data_subset |>
    mutate(bar_order = case_when(name == "Spiroplectammina" ~ 6,
                                 name == "Kotoracythere" ~ 5,
                                 name %in% c("Ostracode", "Foram") ~ 4,
                                 name == "Paracyprideis" ~ 3,
                                 name == "Cassidulina" ~ 2,
                                 name == "Elphidium" ~ 1)) |>
    readr::write_csv(outfile)
  
  return(outfile)
}