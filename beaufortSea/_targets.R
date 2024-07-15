library(targets)
library(tarchetypes)
library(tidyverse)

options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("tidyverse", 
                            "readxl",
                            "packcircles",
                            "ggimage",
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
  Spiroplectammina = "#dd605a", # increases, Foram #"#dd605a"
  Kotoracythere = "#c49051" # increases, Ostracode
)
# Mutate color scheme to long
color_long <- pivot_longer(color_scheme, cols = everything(), values_to = "hexcode") |>
  mutate(hexcode = factor(hexcode, levels = c("#3c475a", "#697a93", "#b4cfd1","#e7f0e7", "#c49051", "#dd605a")))

# Focal species
focal_species <- tibble(
  species = c("E. excavatum clavatu", "Spiroplectammina bif", "Paracyprideis pseudo",
              "Kotoracythere arctob", "C. reniforme...25"),
  species_name = c("E. excavatum", "S. biformis", "P. pseudopunctillata",
                   "K. arctoborealis", "C. reniforme"),
  epithet = c("excavatum", "biformis", "pseudopunctillata",
              "arctoborealis", "reniforme"),
  focal_L = rep(TRUE, 5),
  image_name = c("F_Elphidium.png", "F_Spiroplectammina.png", "O_Paracyprideis.png",
                 "O_Kotoracythere.png", "F_Cassidulina.png"))


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
                            color_long = color_long,
                            focal_species = focal_species)),
  
  # Because foram and ostracode years don't align perfectly, summarize by decade
  tar_target(p2_decade_abundance_long,
             decade_abundance(data_in = p2_join_abundance_long)),
  
  ################ PLOT DATA #########################
  
  # Create assemblage plots as grobs at select years
  tar_map(
    values = tibble(years = c(100, 500, 1000, 1500, 2000)),
    tar_target(p3_assemblage_plots,
               plot_assemblages(data_in = p2_decade_abundance_long, year = years)),
    tar_target(
      p3_assemblage_plot_pngs,
      save_plot(plot_grob = p3_assemblage_plots,
                save_name = sprintf("out/assemblage_%s.png", years),
                width = 900, height = 900),
      format = "file"
    )),
  tar_target(
    p3_viz_thumbnail_png,
    ggsave(p3_assemblage_plots_2000, file = "out/BeaufortSeaTimeline_thumbnail.PNG",
           width = 3, height = 3)
  ),
  
  # Create timeline plot as grob
  tar_target(p3_timeline_plot,
             plot_timeline(data_in = p2_decade_abundance_long)
  ),
  
  ## Species trend plots
  tar_map(
    values = tibble(species = focal_species$epithet),
    tar_target(
      p3_species_trend_plots,
      plot_species_trend(data_in = p2_decade_abundance_long,
                         species_name = species)
    ),
    tar_target(
      p3_species_trend_pngs,
      save_plot(plot_grob = p3_species_trend_plots,
                save_name = sprintf("out/species_trend_%s.png", species),
                width = 1600, height = 900),
      format = "file"
    )
  ),
  
  
  ################ COMPOSE PLOT FOR WEBSITE #############
  
  tar_target(p3_visualization_png,
             compose_timeline(timeline_grob = p3_timeline_plot,
                              assemblage_grob_100 = p3_assemblage_plots_100,
                              assemblage_grob_500 = p3_assemblage_plots_500,
                              assemblage_grob_1000 = p3_assemblage_plots_1000,
                              assemblage_grob_1500 = p3_assemblage_plots_1500,
                              assemblage_grob_2000 = p3_assemblage_plots_2000,
                              color_scheme = color_scheme,
                              png_out = "out/BeaufortSeaTimeline.png"))
  
)

