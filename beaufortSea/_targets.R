library(targets)
library(tarchetypes)
library(tidyverse)

options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("tidyverse", 
                            "readxl",
                            "packcircles",
                            "cowplot",
                            "showtext",
                            "sysfonts"))

source("src/fetch_data.R")
source("src/prep_data.R")
source("src/plot_data.R")
source("src/compose_plot.R")

# Global settings
color_scheme = tibble(
  Cassidulina = "#3c475a", # declines, Foram
  Elphidium = "#697a93", # declines, Foram
  Paracyprideis = "#b4cfd1", # commonly dominant, except LIA
  Ostracode = "#e7f0e7",
  Foram = "#e7f0e7",
  Spiroplectammina = "#dd605a", # increases, Foram
  Kotoracythere = "#c49051" # increases, Ostracode
)
# Mutate color scheme to long
color_long <- pivot_longer(color_scheme, cols = everything(), values_to = "hexcode") |>
  mutate(hexcode = factor(hexcode, levels = c("#3c475a", "#697a93", "#b4cfd1","#e7f0e7", "#c49051", "#dd605a")))

list(
  
  ################ FETCH DATA #########################
  
  # Read in ostracode data
  tar_target(p1_ostracode_raw_csv,
             "in/Data release HLY1302 IP135359.xlsx",
             format = "file"),
  tar_target(p1_ostracode_raw_df,
             read_ostracode(csv_in = p1_ostracode_raw_csv)),
  
  
  # Read in foraminifera data
  tar_target(p1_foram_raw_csv,
             "in/HLY1302 faunal bin foraminifera.xlsx",
             format = "file"),
  tar_target(p1_foram_raw_df,
             read_foram(csv_in = p1_foram_raw_csv)),
  
  ################ PROCESS DATA #########################
  
  # Join abundance data in long format
  tar_target(p2_join_abundance_long,
             join_abundance(ostracode_in = p1_ostracode_raw_df,
                            foram_in = p1_foram_raw_df,
                            color_long = color_long)),
  
  # Because foram and ostracode years don't align perfectly, summarize by decade
  tar_target(p2_decade_abundance_long,
             decade_abundance(data_in = p2_join_abundance_long)),
  
  ################ PLOT DATA #########################
  
  # Create assemblage plots as grobs at select years
  tar_map(
    values = tibble(years = c(100, 500, 1000, 1500, 2000)),
    tar_target(p3_assemblage_plots,
               plot_assemblages(data_in = p2_decade_abundance_long, year = years))),
  
  # Create timeline plot as grob
  tar_target(p3_timeline_plot,
             plot_timeline(data_in = p2_decade_abundance_long)
  ),
  
  
  ################ COMPOSE PLOT FOR WEBSITE #############
  
  tar_target(p3_visualization_png,
             compose_timeline(timeline_grob = p3_timeline_plot,
                              assemblage_grob_100 = p3_assemblage_plots_100,
                              assemblage_grob_500 = p3_assemblage_plots_500,
                              assemblage_grob_1000 = p3_assemblage_plots_1000,
                              color_scheme = color_scheme,
                              png_out = "out/BeaufortSeaTimeline.png"))
  
)

