
#' @description create shape file with HYBAS IDs and their habitat types
#' @param hybas_04_sf shape file of HYBAS IDs
#' @param hybas_legend dataframe with deifning habitat types
get_hybas_habitat_types <- function (hybas_04_sf, hybas_legend) {
  hybas_04_sf |>
    select(HYBAS_ID,fmh_cl_smj) |>
    left_join(hybas_legend, by=c("fmh_cl_smj" = "MHT_ID")) |>
    select(HYBAS_ID, Habitat = MHT_Name)
}

#' @description calculate weighted threat scores
#' @param threat_data dataframe with threat scores
#' @param threat_weights dataframe with weights for each threat type
#' @param hybas_habitat_types shape file with HYBAS IDs and their habitat types
compute_weighted_threats <- function(threat_data, threat_weights,
                                     hybas_habitat_types){
  processed_df <- threat_data |>
    select(HYBAS_ID, ends_with("_LS")) |>
    pivot_longer(cols = ends_with("LS"), names_to = c("ThreatCode", NA), 
                 names_sep = "_", values_to = "ThreatMetric") |>
    left_join(sf::st_drop_geometry(hybas_habitat_types)) |>
    # # Gretchen said there's not much data for Xeric - we could filter these out
    # filter(!Habitat %in% c("Xeric freshwaters and endorheic basins",
    #                        "Greenland",
    #                        "No Data")) |>
    left_join(select(threat_weights, ThreatCode = Threat_Code, Threat, 
                     ThreatCategory = Threat_Category, everything())) |>
    mutate(weightedThreatMetric = ThreatMetric * Final_Weight) |> 
    mutate(ThreatCategory = case_when(ThreatCategory == "Exploitation" ~ "Fishing pressure", 
                                      .default = as.character(ThreatCategory)))
}

#' @description sum weighted threat scores by threat type
#' @param in_dat dataframe with weighted threat scores
#' @param outfile file name and directory for output
compute_total_weighted_threats <- function(in_dat, outfile) {
  in_dat |>
    group_by(Threat, ThreatCategory) |> 
    summarize(TotalWeightedThreatMetric = sum(weightedThreatMetric, na.rm = TRUE)) |>
    select(ThreatCategory, Threat, TotalWeightedThreatMetric) |>
    arrange(desc(TotalWeightedThreatMetric)) |>
    readr::write_csv(outfile)
  
  return(outfile)
}

#' @description mean weighted threat scores by major threat type and HYBAS_ID
#' @param in_dat dataframe with weighted threat scores
compute_mean_weighted_threats <- function(in_dat){ 
  
  processed_df <- in_dat |> 
    group_by(HYBAS_ID, ThreatCategory) |>
    mutate(MeanWeightedThreatMetric = mean(weightedThreatMetric, na.rm = TRUE)) |>
    ungroup() |> 
    select(HYBAS_ID, ThreatCategory, MeanWeightedThreatMetric) |> 
    unique() |> 
    arrange(desc(MeanWeightedThreatMetric))
}

#' @description mean weighted threat scores by sub threat type and HYBAS_ID
#' @param in_dat dataframe with weighted threat scores
compute_mean_weighted_subThreats <- function(in_dat){
  
  processed_df <- in_dat |>
    group_by(HYBAS_ID, ThreatCategory, Threat) |>
    mutate(MeanWeightedThreatMetric = mean(weightedThreatMetric, na.rm = TRUE)) |>
    ungroup() |> 
    select(HYBAS_ID, ThreatCategory, Threat, MeanWeightedThreatMetric) |> 
    rename(MajorCat = ThreatCategory, ThreatCategory = Threat) |> 
    unique() |> 
    arrange(desc(MeanWeightedThreatMetric))
}

#' @description compute top weighted threat category for each HYBAS_ID
#' @param in_dat dataframe with weighted threat scores
#' @param keep_basins_without_scores boolean - should basins where the sum of
#' weighted threats is 0 be kept in this analysis
#' @param category_levels if `keep_basins_without_scores` is TRUE, must provide
#' category levels to use in breaking ties. Otherwise unused
compute_top_weighted_threats <- function(in_dat, keep_basins_without_scores = FALSE, 
                                         category_levels = NULL) {
  
  weighted_sums <- in_dat |>
    # compute sum of weighted threats by threat category
    group_by(HYBAS_ID, ThreatCategory) |> 
    summarize(TotalWeightedThreatMetric = sum(weightedThreatMetric, na.rm = TRUE))
  
  
  # Get threat category with max sum of weighted threats for each basin
  if (keep_basins_without_scores) {
    if (is.null(category_levels)) {
      stop(message("Must provide category levels if keep_basins_without_scores is TRUE"))
    }
    
    unique_categories = unique(pull(weighted_sums, ThreatCategory))
    if (!(all(unique_categories %in% category_levels))) {
      stop(message("All categories must be included in category_levels"))
    }
    
    weighted_sums <- weighted_sums |>
      group_by(HYBAS_ID) |> 
      # if keeping basins where total = 0, set ranking for ties
      mutate(ThreatCategory = factor(ThreatCategory, 
                                     levels = category_levels))|>
      arrange(ThreatCategory)
  } else {
    if (!is.null(category_levels)) {
      warning(message("category_levels will not be used when keep_basins_without_scores is FALSE"))
    }
    
    weighted_sums <- weighted_sums |>
      # filter out basins where total for all threat categories is 0
      filter(sum(TotalWeightedThreatMetric) > 0)
  }
  
  # Get top threat for each HYBAS ID
  top_threats <- weighted_sums |>
    group_by(HYBAS_ID) |> 
    slice_max(
      n = 1,
      order_by = TotalWeightedThreatMetric,
      with_ties = FALSE,
      na_rm = TRUE
    )
}

