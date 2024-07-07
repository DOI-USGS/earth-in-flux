
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

join_genera <- join_long |>
  mutate(genus = word(species)) |>
  mutate(genus = case_when(genus == "E." ~ "Elphidium",
                           genus == "C." ~ "Cassidulina",
                           genus == "Indeterminate/IRO..." ~ "Other",
                           genus == "Other...122" ~ "Other",
                           genus == "Dentalina...39" ~ "Dentalina",
                           genus == "Polymorphinids...41" ~ "Polymorphinids",
                           genus == "Spiroplectammina_spp" ~ "Spiroplectammina",
                           TRUE ~ genus))

#sum_by_genera <- join_genera |>
#  group_by(genus, year, type) |>
#  summarize(genus_abundance = sum(value, na.rm = TRUE))

# Not all years match exactly, so going average by every 50 years and species
join_decades <- join_genera |>
  mutate(decade = year - year %% 100) |>
  group_by(decade, species, type, genus) |>
  summarize(mean_abundance = mean(value, na.rm = TRUE))

# Not all 50-year periods have both species types, so need to add up total abundance by
# period to get actual percent abundance 
total_abundance_by_decade <- join_decades |>
  group_by(decade) |>
  summarize(total_abundance = sum(mean_abundance, na.rm = TRUE)) 

pct_abundance <- join_decades |>
  left_join(total_abundance_by_decade, by = "decade") |>
  mutate(pct_abundance = (mean_abundance/total_abundance)*100)

# list of genera that have no abundance ever
#names_no_abundance <- plot_data |> 
#  group_by(genus) |> 
#  summarize(total = sum(pct_abundance)) |>
#  filter(total == 0)

# Add in hex codes for colors
plot_data <- pct_abundance |>
  mutate(color_hex = case_when(genus == "Cassidulina" ~ "#3c475a", #"#c49051",
                               genus == "Spiroplectammina" ~ "#dd605a",
                               genus == "Elphidium" ~ "#697a93",
                               genus == "Kotoracythere" ~ "#c49051",
                               type == "Ostracode" ~ "#e7f0e7",
                               TRUE ~ "#e7f0e7"),
         name = case_when(genus == "Spiroplectammina" ~ "S. spp",
                          genus == "Cassidulina" ~ "C. spp",
                          genus == "Elphidium" ~ "E. spp",
                          genus == "Kotoracythere" ~ "K. spp",
                          type == "Ostracode" ~ "Ostracode",
                          type == "Foram" ~ "Foram",
                          TRUE ~ NA))

plot_data <- plot_data |>
  mutate(color_factor = factor(color_hex, ordered = FALSE,
                               levels = c("#697a93", "#3c475a", "grey90", "#dd605a", "#c49051")))


ggplot(data = plot_data, aes(y = decade, x = pct_abundance, fill = color_factor)) +
  geom_bar(position = "stack", stat = "identity", orientation = "y",
           width = 15, color = NA) +
  scale_fill_identity() +
  scale_y_reverse() +
  theme_minimal() +
  ylab("Year") +
  xlab("Relative Abundance")


ggsave("beaufort_althea/out/timeline_bar_genera.png",
       height = 20, width = 4)






############ old versions

# 
# 
# 
# # Drop NAs and any species with 0% abundance
# plot_data_noNA <- plot_data |>
#   filter(pct_abundance > 0) #|>
#   #drop_na()
# 
# # Number by alphabetic species names
# plot_alphabetic <- plot_data |>
#   arrange(genus) |>
#   group_by(decade) |>
#   mutate(order = row_number())
# 
# plot_alpha_noNA <- plot_alphabetic |>
#   filter(pct_abundance > 0) |>
#   drop_na()
# 
# # Plot
# ggplot(data = plot_alpha_noNA, aes(y = decade, x = order)) +
#   geom_point(aes(size = pct_abundance, color = color_hex)) +
#   scale_color_identity() +
#   #scale_size_continuous(limits = c(1, 5)) +
#   theme_bw() +
#   theme(legend.position = "none")
# ggsave("beaufort_althea/out/timeline_assemblage_genera.png",
#        height = 20, width = 4)
# 
# 
# 
# ################# FORAMS ONLY ###################
# 
# forams <- foram_data |>
#   rename(year = 'calendar yr') |>
#   mutate(year = round(year, 0))
# 
# join_long <- forams |>
#   pivot_longer(-year, names_to = "species") |>
#   mutate(type = case_when(species %in% names(forams) ~ "Foram",
#                           TRUE ~ NA)) 
# 
# join_genera <- join_long |>
#   mutate(genus = word(species)) |>
#   mutate(genus = case_when(genus == "E." ~ "Elphidium",
#                            genus == "C." ~ "Cassidulina",
#                            genus == "Indeterminate/IRO..." ~ "Other",
#                            genus == "Other...122" ~ "Other",
#                            genus == "Dentalina...39" ~ "Dentalina",
#                            genus == "Polymorphinids...41" ~ "Polymorphinids",
#                            TRUE ~ genus))
# 
# sum_by_genera <- join_genera |>
#   group_by(genus, year, type) |>
#   summarize(genus_abundance = sum(value, na.rm = TRUE))
# 
# # Not all years match exactly, so going average by every 50 years and species
# join_decades <- sum_by_genera |>
#   mutate(decade = year - year %% 50) |>
#   group_by(decade, genus, type) |>
#   summarize(mean_abundance = mean(genus_abundance, na.rm = TRUE))
# 
# # Not all 50-year periods have both species types, so need to add up total abundance by
# # period to get actual percent abundance 
# total_abundance_by_decade <- join_decades |>
#   group_by(decade) |>
#   summarize(total_abundance = sum(mean_abundance, na.rm = TRUE)) 
# 
# pct_abundance <- join_decades |>
#   left_join(total_abundance_by_decade, by = "decade") |>
#   mutate(pct_abundance = (mean_abundance/total_abundance)*100)
# 
# # list of genera that have no abundance ever
# names_no_abundance <- plot_data |> 
#   group_by(genus) |> 
#   summarize(total = sum(pct_abundance)) |>
#   filter(total == 0)
# 
# # Add in hex codes for colors
# plot_data <- pct_abundance |>
#   filter(! genus %in% names_no_abundance$genus) |>
#   mutate(color_hex = case_when(genus == "Kotoracythere" ~ "#c49051",
#                                genus == "Spiroplectammina" ~ "#dd605a",
#                                genus == "Elphidium" ~ "#697a93",
#                                genus == "Paracyprideis" ~ "#8c939e",#3c475a
#                                genus == "Cassidulina" ~ "#3c475a",
#                                type == "Ostracode" ~ "#B6B2A6",
#                                TRUE ~ "#F6E6CB"),
#          name = case_when(genus == "Kotoracythere" ~ genus,
#                           genus == "Elphidium" ~ "Other Foram",
#                           genus == "Spiroplectammina" ~ genus,
#                           genus == "Paracyprideis" ~ genus,
#                           genus == "Cassidulina" ~ genus,
#                           type == "Ostracode" ~ "Ostracode",
#                           type == "Foram" ~ "Other Foram",
#                           TRUE ~ NA))
# 
# ggplot(data = plot_data, aes(x = decade, y = pct_abundance, fill = name)) +
#   geom_bar(position = "stack", stat = "identity",
#            width = 15, color = NA) 
# 
# ggplot(data = plot_data |>
#          filter(genus == "Cassidulina"), aes(x = decade, y = pct_abundance, fill = name)) +
#   geom_bar(position = "stack", stat = "identity",
#            width = 15, color = NA) +
#   geom_smooth(method = "lm")
# ggplot(data = plot_data |>
#          filter(genus == "Spiroplectammina"), aes(x = decade, y = pct_abundance, fill = name)) +
#   geom_bar(position = "stack", stat = "identity",
#            width = 15, color = NA) +
#   geom_smooth(method = "lm")

