clean_input_data <- function(data_file) {
  data <- read_xlsx(data_file, na="NA")
  
  names(data)<-make.names(names(data),unique = TRUE)
  
  data |>
    mutate(
      family = ifelse(family == 'Salmonindae', 'Salmonidae', family),
      species_common = ifelse(species_common == 'Bey?ehir bleak', 
                              'Beyşehir bleak', 
                              species_common)
      )
}

build_nested_json <- function(data, focal_columns, out_file) {
  data_list <- data |>
    filter(!is.na(!!rlang::sym(focal_columns[['value']]))) |>
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