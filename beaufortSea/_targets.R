library(targets)
library(tarchetypes)
library(tidyverse)

###############################################################################
#########       TO RUN PIPELINE FOR THE FIRST TIME
#########         Download data from Sciencebase, using the code in 00_config.R
#########         Then, run `targets::tar_make()`

###############################################################################
##########      RUN PIPELINE
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


# Global settings
color_scheme = tibble(
  Cassidulina = "#3c475a", # declines, Foram
  Elphidium = "#697a93", # declines, Foram
  Paracyprideis = "#729C9D", # commonly dominant, except LIA
  Ostracode = "#e7f0e7",
  Foram = "#e7f0e7",
  Spiroplectammina = "#dd605a", # increases, Foram #"#dd605a"
  Kotoracythere = "#c49051" # increases, Ostracode
)
# Mutate color scheme to long
color_long <- pivot_longer(color_scheme, cols = everything(), values_to = "hexcode") |>
  mutate(hexcode = factor(hexcode, levels = c("#3c475a", "#697a93", "#729C9D","#e7f0e7", "#c49051", "#dd605a")))

# Load some custom fonts and set some custom settings
title_font <- "Lora"
sysfonts::font_add_google(title_font)
supporting_font <- "Source Sans 3"
sysfonts::font_add_google(supporting_font)
annotation_font <- "Caveat"
sysfonts::font_add_google(annotation_font)

# Focal species
focal_species <- tribble(
  ~species, ~species_name, ~epithet, ~focal_L, ~image_name, ~image_size, ~image_y, ~ylim, ~fig_number,
  "Elphidium excavatum clavatum","E. excavatum","excavatum",TRUE,"F_Elphidium.png",0.2,48,58,"2b",
  "Spiroplectammina biformis","S. biformis","biformis",TRUE,"F_Spiroplectammina.png",0.15,55,70,"2c",
  "Paracyprideis pseudo","P. pseudopunctillata","pseudopunctillata",TRUE,"O_Paracyprideis.png",0.3,100,110,"3b", 
  "Kotoracythere arctob","K. arctoborealis","arctoborealis",TRUE,"O_Kotoracythere.png",0.3,25,30,"3a",
  "Cassidulina reniforme","C. reniforme","reniforme",TRUE,"F_Cassidulina.png",0.2,80,90,"2a"
)


list(
  
  ################ FETCH DATA #########################
  
  # Read in ostracode data from ScienceBase ID 6197a410d34eb622f692ce0e
  tar_target(p1_ostracode_raw_xlsx, 
             sb_initialize_and_download(
               sb_id ="6197a410d34eb622f692ce0e",  
               names = "Data release HLY1302 IP135359.xlsx",
               destinations = "in/Data release HLY1302 IP135359.xlsx",
               overwrite_fileL = FALSE
             ),
             format = "file"),
  # Read ostracode data
  tar_target(p1_ostracode_raw_df,
             read_ostracode(xlsx_in = p1_ostracode_raw_xlsx)),
  # read in sheet with age model information
  tar_target(p1_age_model_data_df,
             readr::read_csv(file = "in/temp_crosswalk.csv",
                             show_col_types = FALSE)),
  
  
  # Read in foraminifera data (downloaded from supplementary information)
  tar_target(p1_foram_raw_xlsx,
             "in/41063_2018_58_MOESM1_ESM.xlsx",
             format = "file"),
  # one page from excel for each of the three cores
  tar_target(p1_foram_HLY1302_MC29_raw_df,
             read_excel(path = p1_foram_raw_xlsx, 
                        sheet = "HLY1302 MC29 Foram",
                        range = "C1:AJ46")),
  tar_target(p1_foram_HLY1302_GGC30_raw_df,
             read_excel(path = p1_foram_raw_xlsx, 
                        sheet = "HLY1302 GGC30 Foram ",
                        range = "C1:AM57")),
  tar_target(p1_foram_HLY1302_JPC32_raw_df,
             read_excel(path = p1_foram_raw_xlsx, 
                        sheet = "HLY1302 JPC32 Foram",
                        range = "C1:AJ142")),
  
  ################ PROCESS DATA #########################
  # Merge foram data (raw counts)
  tar_target(p2_foram_counts_df,
             merge_foram_data(GGC30_in = p1_foram_HLY1302_GGC30_raw_df, 
                              JPC32_in = p1_foram_HLY1302_JPC32_raw_df,
                              MC29_in = p1_foram_HLY1302_MC29_raw_df,
                              age_data = p1_age_model_data_df)),
  # Clean foram data
  tar_target(p2_foram_df,
             clean_foram_data(raw_in = p2_foram_counts_df)),
  
  

  
  # Join abundance data in long format
  tar_target(p2_join_abundance_long,
             join_abundance(ostracode_in = p1_ostracode_raw_df,
                            foram_in = p2_foram_df,
                            color_long = color_long,
                            focal_species = focal_species)),
  
  # Because foram and ostracode years don't align perfectly, summarize by decade
  tar_target(p2_decade_abundance_long,
             decade_abundance(data_in = p2_join_abundance_long)),
  
  tar_target(p2_decade_abundance_csv,
             generate_beeswarm_data(data = p2_decade_abundance_long,
                                    outfile = '../public/beaufort_species_abundance.csv'),
             format = 'file'),
  
  ################ PLOT DATA #########################
  
  ## Species trend plots
  tar_map(
    values = tibble(species = focal_species$epithet,
                    fig_num = focal_species$fig_number),
    tar_target(
      p3_species_trend_plots,
      plot_species_trend(data_in = p2_join_abundance_long,
                         species_name = species)
    ),
    tar_target(
      p3_species_trend_pngs,
      save_plot(plot_grob = p3_species_trend_plots,
                save_name = sprintf("images/BeaufortSeaSpecies_%s.png", fig_num),
                width = 1600, height = 900),
      format = "file"
    )
  )
)

