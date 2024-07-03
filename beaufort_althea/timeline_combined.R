
library(readxl)
library(tidyverse)

# Assemblage bee swarm data - ostrocodes
raw_pct_abundance <- read_excel(path = "beaufort_althea/in/Data release HLY1302 IP135359.xlsx", 
                                sheet = "T2 Ost MC29bin5G30bin10J32bin20",
                                skip = 1, 
                                range = cell_cols(c(7, 70:131)))


abundance_data <- raw_pct_abundance[,c(1,65:123)] 
names(abundance_data) <- substr(names(abundance_data), 1, 20)

# Assemblage bee swarm data - forams
forams_raw <- read_excel(path = "beaufort_althea/in/HLY1302 faunal bin foraminifera.xlsx", 
                         sheet = "HLY1302 FORAM counts as valu",
                         range = cell_cols(c(7, 70:131)))

foram_data <- forams_raw[c(1:25, 27:120),c(2,24:42)] 
names(foram_data) <- substr(names(foram_data), 1, 20)

## Combine the two Spiroplectammina spp
foram_data <- foram_data |>
  rename(s_bif = `Spiroplectammina bif`,
         s_ear = `Spiroplectammina ear`) |>
  mutate(Spiroplectammina_spp = s_bif + s_ear) |>
  select(-s_bif, -s_ear)

## Join data
ostracodes <- abundance_data |>
  rename(year = `Calendar AGE`) |>
  mutate(year = round(year, 0))

forams <- foram_data |>
  rename(year = 'calendar yr') |>
  mutate(year = round(year, 0))

join_long <- ostracodes |>
  full_join(forams, by = "year") |>
  pivot_longer(-year, names_to = "species") |>
  mutate(type = case_when(species %in% names(ostracodes) ~ "Ostracode",
                          species %in% names(forams) ~ "Foram",
                          TRUE ~ NA))

# Not all years match exactly, so going average by every 50 years and species
join_decades <- join_long |>
  mutate(decade = year - year %% 50) |>
  group_by(decade, species, type) |>
  summarize(mean_abundance = mean(value, na.rm = TRUE))

# Not all 50-year periods have both species types, so need to add up total abundance by
# period to get actual percent abundance 
total_abundance_by_decade <- join_decades |>
  group_by(decade) |>
  summarize(total_abundance = sum(mean_abundance, na.rm = TRUE)) 

pct_abundance <- join_decades |>
  left_join(total_abundance_by_decade, by = "decade") |>
  mutate(pct_abundance = (mean_abundance/total_abundance)*100)


# Add in hex codes for colors
plot_data <- pct_abundance |>
  mutate(color_hex = case_when(species == "Kotoracythere arctob" ~ "#c49051",
                               species == "Spiroplectammina_spp" ~ "#dd605a",
                               species == "E. excavatum clavatu" ~ "#697a93",
                               species == "Paracyprideis pseudo" ~ "#8c939e",#3c475a
                               species == "C. reniforme...25" ~ "#3c475a",
                               type == "Ostracode" ~ "#B6B2A6",
                               TRUE ~ "#F6E6CB"),
         timeperiod_group = case_when(species %in% c("Kotoracythere arctob", 
                                                     "Spiroplectammina_spp") ~ "Present",
                                      species %in% c("E. excavatum clavatu",
                                                     "Paracyprideis pseudo",
                                                     "C. reniforme...25") ~ "Past",
                                      TRUE ~ "Other"),
         species_group = case_when(species %in% c("Kotoracythere arctob", 
                                                  "Paracyprideis pseudo") ~ "Ostracode",
                                   species %in% c("E. excavatum clavatu", 
                                                  "Spiroplectammina_spp",
                                                  "C. reniforme...25") ~ "Foram",
                                   TRUE ~ "Other"),
         name = case_when(species == "Kotoracythere arctob" ~ "K. arctoborealis",
                          species == "E. excavatum clavatu" ~ "E. excavatum",
                          species == "Spiroplectammina_spp" ~ "S. spp",
                          species == "Paracyprideis pseudo" ~ "P. pseudopunctillata",
                          species == "C. reniforme...25" ~ "C. reniforme",
                          type == "Ostracode" ~ "Ostracode",
                          type == "Foram" ~ "Foram",
                          TRUE ~ NA))


ggplot(data = plot_data, aes(x = decade, y = pct_abundance, fill = timeperiod_group)) +
  geom_bar(position = "stack", stat = "identity",
           width = 15, color = NA)

# Drop NAs and any species with 0% abundance
plot_data_noNA <- plot_data |>
  filter(pct_abundance > 0) #|>
  #drop_na()

# Number by alphabetic species names
plot_alphabetic <- plot_data |>
  arrange(species) |>
  group_by(decade) |>
  mutate(order = row_number())

plot_alpha_noNA <- plot_alphabetic |>
  filter(pct_abundance > 0) |>
  drop_na()

# Plot
ggplot(data = plot_alpha_noNA, aes(y = decade, x = order)) +
  geom_point(aes(size = pct_abundance, color = color_hex)) +
  scale_color_identity()
