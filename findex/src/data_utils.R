get_hybas_habitat_types <- function (hybas_04_sf, hybas_legend) {
  hybas_04_sf |>
    select(HYBAS_ID,fmh_cl_smj) |>
    left_join(hybas_legend, by=c("fmh_cl_smj" = "MHT_ID")) |>
    select(HYBAS_ID, Habitat = MHT_Name)
}

compute_total_weighted_threats <- function(threat_data, threat_weights,
                                           hybas_habitat_types, outfile) {
  threat_data |>
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
    group_by(Threat, ThreatCategory) |>
    summarize(TotalWeightedThreatMetric = sum(weightedThreatMetric, na.rm = TRUE)) |>
    select(ThreatCategory, Threat, TotalWeightedThreatMetric) |>
    arrange(desc(TotalWeightedThreatMetric)) |>
    readr::write_csv(outfile)
  
  return(outfile)
}

compute_mean_weighted_threats <- function(threat_data, threat_weights, 
                                     hybas_habitat_types){ #outfile
  
  processed_df <- threat_data |> 
    select(HYBAS_ID, ends_with("_LS")) |>
    pivot_longer(cols = ends_with("LS"), names_to = c("ThreatCode", NA), 
                 names_sep = "_", values_to = "ThreatMetric") |>
    left_join(hybas_habitat_types) |>
    # # Gretchen said there's not much data for Xeric - we could filter these out
    # filter(!Habitat %in% c("Xeric freshwaters and endorheic basins",
    #                        "Greenland",
    #                        "No Data")) |>
    left_join(select(threat_weights, ThreatCode = Threat_Code, Threat, 
                     ThreatCategory = Threat_Category, everything())) |>
    mutate(weightedThreatMetric = ThreatMetric * Final_Weight) |>
    group_by(HYBAS_ID, ThreatCategory) |>
    mutate(MeanWeightedThreatMetric = mean(weightedThreatMetric, na.rm = TRUE)) |>
    ungroup() |> 
    select(HYBAS_ID, ThreatCategory, MeanWeightedThreatMetric, Shape) |>
    unique() |> 
    arrange(desc(MeanWeightedThreatMetric)) #|>
    #readr::write_csv(outfile)
  
  #return(outfile)
}

# same as comput_mean_weighted_threats except it includes subcategories
compute_mean_weighted_subThreats <- function(threat_data, threat_weights, 
                                             hybas_habitat_types){
  
  processed_df <- threat_data |> 
    select(HYBAS_ID, ends_with("_LS")) |>
    pivot_longer(cols = ends_with("LS"), names_to = c("ThreatCode", NA), 
                 names_sep = "_", values_to = "ThreatMetric") |>
    left_join(hybas_habitat_types) |>
    # # Gretchen said there's not much data for Xeric - we could filter these out
    # filter(!Habitat %in% c("Xeric freshwaters and endorheic basins",
    #                        "Greenland",
    #                        "No Data")) |>
    left_join(select(threat_weights, ThreatCode = Threat_Code, Threat, 
                     ThreatCategory = Threat_Category, everything())) |>
    mutate(weightedThreatMetric = ThreatMetric * Final_Weight) |>
    group_by(HYBAS_ID, ThreatCategory, Threat) |>
    mutate(MeanWeightedThreatMetric = mean(weightedThreatMetric, na.rm = TRUE)) |>
    ungroup() |> 
    select(HYBAS_ID, ThreatCategory, Threat, MeanWeightedThreatMetric, Shape) |>
    unique() |> 
    arrange(desc(MeanWeightedThreatMetric))
}
