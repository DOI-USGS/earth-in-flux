#' Initialize ScienceBase session and download files (large zipped files only)
#'
#' @param sb_id chr; ScienceBase ID
#' @param unzip_file_to_check chr; names of unzipped child targets to check against
#' @param names chr; names of files to download from ScienceBase
#' @param destinations  chr; write path location for downloaded files
#' @param overwrite_destinationL logical, do you want to overwrite file?
#' @param renviron_file chr; path to .Renviron file where credentials are cached
#' @param ... additional arguments passed to `sbtools::item_file_download()`
#'
#' @return chr; path to downloaded files
#' 
sb_initialize_and_download_zipped <- function(sb_id, 
                                              unzip_file_to_check, 
                                              names, 
                                              destination_zip,
                                              download_dir, 
                                              overwrite_fileL,
                                              renviron_file = ".Renviron", ...) {
  
  # Does destination file already exist?
  does_file_exist <- file.exists(unzip_file_to_check)
  
  # If the destination file does exist, check dates before downloading
  if(does_file_exist){
    
    # Initialize ScienceBase session
    sb_login_cached(renviron_file = renviron_file)
    
    item_metadata <- sbtools::item_get(sb_id)$files
    
    file_index <- purrr::map(item_metadata, `[`, "name") |>
      unlist(use.names = FALSE) |>
      purrr::map_lgl(~ .x == names) |>
      which()
    
    file_metadata <- item_metadata[[file_index]]
    
    # If no hash is available check dates:
    # Get upload date from SB
    ul_time <- file_metadata$dateUploaded |>
      as.POSIXct(format = "%Y-%m-%dT%T", tz = "UTC")
    
    # Get download date
    dl_time <- file.info(unzip_file_to_check)$mtime |>
      `attr<-`("tzone", "UTC")
    
    download_needed <- ul_time > dl_time
    
    if(download_needed){
      # Download SB files
      sbtools::item_file_download(
        sb_id = sb_id,
        names = names,
        destinations = destination_zip,
        overwrite_file = TRUE,
        ...
      )
      
      unzip(zipfile = destination_zip, 
            overwrite = overwrite_fileL,
            exdir = download_dir)
      
      file.remove(destination_zip)
    }
  } else {
    # file doesn't yet exist, so download
    sbtools::item_file_download(
      sb_id = sb_id,
      names = names,
      destinations = destination_zip,
      overwrite_file = overwrite_fileL,
      ...
    )
    
    unzip(zipfile = destination_zip, 
          overwrite = overwrite_fileL,
          exdir = download_dir)
    
    file.remove(destination_zip)
  }
  
  return(unzip_file_to_check)
  
}


#' Initialize ScienceBase session and download files (large zipped files only)
#'
#' @param sb_id chr; ScienceBase ID
#' @param unzip_file_to_check chr; names of unzipped child targets to check against
#' @param names chr; names of files to download from ScienceBase
#' @param destinations  chr; write path location for downloaded files
#' @param overwrite_destinationL logical, do you want to overwrite file?
#' @param renviron_file chr; path to .Renviron file where credentials are cached
#' @param ... additional arguments passed to `sbtools::item_file_download()`
#'
#' @return chr; path to downloaded files
#' 
sb_initialize_and_download <- function(sb_id, 
                                       names, 
                                       destinations, 
                                       overwrite_fileL,
                                       renviron_file = ".Renviron", ...) {
  
  # Does destination file already exist?
  does_file_exist <- file.exists(destinations)
  
  # Initialize ScienceBase session
  sb_login_cached(renviron_file = renviron_file)
  
  if(does_file_exist){
    
    
    item_metadata <- sbtools::item_get(sb_id)$files
    
    file_index <- purrr::map(item_metadata, `[`, "name") |>
      unlist(use.names = FALSE) |>
      purrr::map_lgl(~ .x == names) |>
      which()
    
    file_metadata <- item_metadata[[file_index]]
    
    # If no hash is available check dates:
    # Get upload date from SB
    ul_time <- file_metadata$dateUploaded |>
      as.POSIXct(format = "%Y-%m-%dT%T", tz = "UTC")
    
    # Get download date
    dl_time <- file.info(destinations)$mtime |>
      `attr<-`("tzone", "UTC")
    
    download_needed <- ul_time > dl_time
    
    if(download_needed){
      
      
      
      # Download SB files
      sbtools::item_file_download(
        sb_id = sb_id,
        names = names,
        destinations = destinations,
        overwrite_file = TRUE
      )
      
    }
  } else {
    
    # file doesn't yet exist, so download
    sbtools::item_file_download(
      sb_id = sb_id,
      names = names,
      destinations = destinations,
      overwrite_file = overwrite_fileL
    )
  }
  
  return(destinations)
  
}


#' Login to ScienceBase using cached credentials
#'
#' @param renviron_file chr; path to .Renviron file
#'
#' @return `TRUE` if logged in. Error if not.
#' 
sb_login_cached <- function(renviron_file) {
  # If logged in, return TRUE and skip the rest
  if(sbtools::is_logged_in()) {
    return(TRUE)
  }
  
  # Try a token refresh
  tryCatch(
    sbtools:::token_refresh(),
    warning = function(x) {},
    error = function(x) FALSE
  )
  
  if(sbtools::is_logged_in()) {
    return(TRUE)
  }
  
  # If .Renviron file does not exist, re-initialize
  if(!file.exists(renviron_file)) {
    cli::cli_abort(c(
      "Could not find the specified file: {.file {renviron_file}}.",
      "i" = "Follow the instructions in {.file README.md} to initalize and cache ScienceBase login credentials."
    ))
  }
  
  # Read .Renviron file
  existing <- readLines(renviron_file)
  sb_token_idx <- which(stringr::str_detect(existing, "^sb_token="))
  sb_username_idx <- which(stringr::str_detect(existing, "^sb_username="))
  
  # If SB credentials not found, throw error
  if(any(length(sb_token_idx) == 0, length(sb_username_idx) == 0)) {
    cli::cli_abort(c(
      "Could not find the username or token in the specified file: {.file {renviron_file}}.",
      "i" = "Follow the instructions in {.file README.md} to re-initalize and cache ScienceBase login credentials."
    ))
  }
  
  # Get ScienceBase credentials
  sb_token <- stringr::str_remove(existing[sb_token_idx], "^sb_token=")
  sb_username <- stringr::str_remove(existing[sb_username_idx], "^sb_username=")
  
  # Initialize ScienceBase session with cached credentials
  sbtools::initialize_sciencebase_session(
    username = sb_username,
    token_text = sb_token
  )
  
  if(sbtools::is_logged_in()) {
    return(TRUE)
  } else {
    cli::cli_abort(c(
      "Could not login to ScienceBase using cached credentials.",
      "i" = "Follow the instructions in {.file README.md} to re-initalize and cache ScienceBase login credentials."
    ))
  }
}




read_ostracode <- function(xlsx_in){
  ostracode_raw <- read_excel(path = xlsx_in, 
                                  sheet = "T2 Ost MC29bin5G30bin10J32bin20",
                                  skip = 1, 
                                  range = cell_cols(c(7, 70:131)))
  # clean up column names
  ostracode_data <- ostracode_raw[,c(1,65:123)]
  names(ostracode_data) <- substr(names(ostracode_data), 1, 20)
  
  return(ostracode_data)
}

read_age_model <- function(xlsx_in){
  ostracode_raw <- read_excel(path = xlsx_in, 
                              sheet = "T2 Ost MC29bin5G30bin10J32bin20",
                              skip = 1, 
                              range = cell_cols(c(7, 70:131)))
  # clean up column names
  ostracode_data <- ostracode_raw[,c(1,65:123)]
  names(ostracode_data) <- substr(names(ostracode_data), 1, 20)
  
  return(ostracode_data)
}