
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
  mutate(genus = word(species)) |>
  mutate(type = case_when(species %in% names(ostracodes) ~ "Ostracode",
                          species %in% names(forams) ~ "Foram",
                          TRUE ~ NA)) |>
  mutate(genus = case_when(genus == "E." ~ "Elphidium",
                           genus == "C." ~ "Cassidulina",
                           genus == "Indeterminate/IRO..." ~ "Other",
                           genus == "Other...122" ~ "Other",
                           genus == "Dentalina...39" ~ "Dentalina",
                           genus == "Polymorphinids...41" ~ "Polymorphinids",
                           genus == "Spiroplectammina_spp" ~ "Spiroplectammina",
                           TRUE ~ genus))

# Not all years match exactly, so going average by every 100 years and species
join_centuries <- join_long |>
  mutate(century = year - year %% 100) |>
  group_by(century, species, type, genus) |>
  summarize(mean_abundance = mean(value, na.rm = TRUE))

# Add in hex codes for colors
plot_data <- join_centuries |>
  mutate(color_hex = case_when(genus == "Cassidulina" ~ "#3c475a", #"#c49051",
                               genus == "Spiroplectammina" ~ "#dd605a",
                               genus == "Elphidium" ~ "#697a93",
                               genus == "Kotoracythere" ~ "#c49051",
                               #genus == "Paracyprideis" ~ "blue",
                               #genus == "Semicytherura" ~ "purple",
                               #species == "E. excavatum clavatu" ~ "#697a93",
                               #species == "Paracyprideis pseudo" ~ "#8c939e",#3c475a
                               type == "Ostracode" ~ "#B6B2A6",
                               TRUE ~ "#F6E6CB"),
         timeperiod_group = case_when(genus == "Spiroplectammina" ~ "Present",
                                      genus == "Cassidulina" ~ "Past",
                                      TRUE ~ "Other"),
         name = case_when(genus == "Spiroplectammina" ~ "S. spp",
                          genus == "Cassidulina" ~ "C. spp",
                          genus == "Elphidium" ~ "E. spp",
                          genus == "Kotoracythere" ~ "K. spp",
                          #genus == "Paracyprideis" ~ "Para",
                          #genus == "Semicytherura" ~ "Semicyth",
                          type == "Ostracode" ~ "Ostracode",
                          type == "Foram" ~ "Foram",
                          TRUE ~ NA))

# circles by century
library(packcircles)
library(ggimage)

y100 <- plot_data |> filter(century == 100)
y100_pack <- circleProgressiveLayout(y100$mean_abundance, sizetype = "area")
y100_join <- cbind(y100, y100_pack)
y100_plot <- circleLayoutVertices(y100_pack, npoints = 50)
y100_plot$species <- rep(y100$species, each = 51)
y100_plot$color <- rep(y100$color_hex, each = 51)
y100_plot$timeperiod_group <- rep(y100$timeperiod_group, each = 51)
y100_plot$name <- rep(y100$name, each = 51)
y100_plot$genus <- rep(y100$genus, each = 51)

# Make the plot
ggplot() + 
  # Make the bubbles
  geom_polygon(data = y100_plot, 
               aes(x, y, group = id, fill = color), 
               colour = "grey80", linewidth = 0.3) +
  # Scatter plot 
  # ggimage::geom_image(
  #   data = y100_join |> filter(name == "S. spp"),
  #   aes(x, y),
  #   size = 0.01,
  #   image = "beaufort_althea/S_app.png",
  #   # Aspect ratio
  #   asp = 1.9) +
  # Add text in the center of each bubble + control its size
  geom_text(data = y100_join, aes(x, y, size = mean_abundance, label = name)) +
  scale_size_continuous(range = c(1,6)) +
  scale_fill_identity() +
  # General theme:
  theme_void() + 
  theme(legend.position="none") +
  coord_equal()

ggsave("beaufort_althea/out/y0100.png")


y500 <- plot_data |> filter(century == 500)
y500_pack <- circleProgressiveLayout(y500$mean_abundance, sizetype = "area")
y500_join <- cbind(y500, y500_pack)
y500_plot <- circleLayoutVertices(y500_pack, npoints = 50)
y500_plot$species <- rep(y500$species, each = 51)
y500_plot$color <- rep(y500$color_hex, each = 51)
y500_plot$timeperiod_group <- rep(y500$timeperiod_group, each = 51)
y500_plot$name <- rep(y500$name, each = 51)
y500_plot$genus <- rep(y500$genus, each = 51)


# Make the plot
ggplot() + 
  # Make the bubbles
  geom_polygon(data = y500_plot, 
               aes(x, y, group = id, fill = color), 
               colour = "grey80", linewidth = 0.3) +
  # Add text in the center of each bubble + control its size
  geom_text(data = y500_join, aes(x, y, size = mean_abundance, label = name)) +
  scale_size_continuous(range = c(1,6)) +
  scale_fill_identity() +
  # General theme:
  theme_void() + 
  theme(legend.position="none") +
  coord_equal()

ggsave("beaufort_althea/out/y0500.png")

y1000 <- plot_data |> filter(century == 1000)
y1000_pack <- circleProgressiveLayout(y1000$mean_abundance, sizetype = "area")
y1000_join <- cbind(y1000, y1000_pack)
y1000_plot <- circleLayoutVertices(y1000_pack, npoints = 50)
y1000_plot$species <- rep(y1000$species, each = 51)
y1000_plot$color <- rep(y1000$color_hex, each = 51)
y1000_plot$timeperiod_group <- rep(y1000$timeperiod_group, each = 51)
y1000_plot$name <- rep(y1000$name, each = 51)
y1000_plot$genus <- rep(y1000$genus, each = 51)


# Make the plot
ggplot() + 
  # Make the bubbles
  geom_polygon(data = y1000_plot, 
               aes(x, y, group = id, fill = color), 
               colour = "grey80", linewidth = 0.3) +
  # Add text in the center of each bubble + control its size
  geom_text(data = y1000_join, aes(x, y, size = mean_abundance, label = name)) +
  scale_size_continuous(range = c(1,4)) +
  scale_fill_identity() +
  # General theme:
  theme_void() + 
  theme(legend.position="none") +
  coord_equal()

ggsave("beaufort_althea/out/y1000.png")

y1500 <- plot_data |> filter(century == 1500)
y1500_pack <- circleProgressiveLayout(y1500$mean_abundance, sizetype = "area")
y1500_join <- cbind(y1500, y1500_pack)
y1500_plot <- circleLayoutVertices(y1500_pack, npoints = 50)
y1500_plot$species <- rep(y1500$species, each = 51)
y1500_plot$color <- rep(y1500$color_hex, each = 51)
y1500_plot$timeperiod_group <- rep(y1500$timeperiod_group, each = 51)
y1500_plot$name <- rep(y1500$name, each = 51)
y1500_plot$genus <- rep(y1500$genus, each = 51)




# Make the plot
ggplot() + 
  # Make the bubbles
  geom_polygon(data = y1500_plot, 
               aes(x, y, group = id, fill = color), 
               colour = "grey80", linewidth = 0.3) +
  # Add text in the center of each bubble + control its size
  geom_text(data = y1500_join, aes(x, y, size = mean_abundance, label = name)) +
  scale_size_continuous(range = c(1,4)) +
  scale_fill_identity() +
  # General theme:
  theme_void() + 
  theme(legend.position="none") +
  coord_equal()

ggsave("beaufort_althea/out/y1500.png")

y1900 <- plot_data |> filter(century == 1900)
y1900_pack <- circleProgressiveLayout(y1900$mean_abundance, sizetype = "area")
y1900_join <- cbind(y1900, y1900_pack)
y1900_plot <- circleLayoutVertices(y1900_pack, npoints = 50)
y1900_plot$species <- rep(y1900$species, each = 51)
y1900_plot$color <- rep(y1900$color_hex, each = 51)
y1900_plot$timeperiod_group <- rep(y1900$timeperiod_group, each = 51)
y1900_plot$name <- rep(y1900$name, each = 51)
y1900_plot$genus <- rep(y1900$genus, each = 51)



# Make the plot
ggplot() + 
  # Make the bubbles
  geom_polygon(data = y1900_plot, 
               aes(x, y, group = id, fill = color), 
               colour = "grey80", linewidth = 0.3) +
  # Add text in the center of each bubble + control its size
  geom_text(data = y1900_join, aes(x, y, 
                                   size = mean_abundance, 
                                   label = name)) +
  scale_size_continuous(range = c(1,4)) +
  scale_fill_identity() +
  # General theme:
  theme_void() + 
  theme(legend.position="none") +
  coord_equal()

ggsave("beaufort_althea/out/y1900.png")

y2000 <- plot_data |> filter(century == 2000)
y2000_pack <- circleProgressiveLayout(y2000$mean_abundance, sizetype = "area")
y2000_join <- cbind(y2000, y2000_pack)
y2000_plot <- circleLayoutVertices(y2000_pack, npoints = 50)
y2000_plot$species <- rep(y2000$species, each = 51)
y2000_plot$color <- rep(y2000$color_hex, each = 51)
y2000_plot$timeperiod_group <- rep(y2000$timeperiod_group, each = 51)
y2000_plot$name <- rep(y2000$name, each = 51)
y2000_plot$genus <- rep(y2000$genus, each = 51)


# Make the plot
ggplot() + 
  # Make the bubbles
  geom_polygon(data = y2000_plot, 
               aes(x, y, group = id, fill = color), 
               colour = "grey80", linewidth = 0.3) +
  # Add text in the center of each bubble + control its size
  geom_text(data = y2000_join, aes(x, y, 
                                   size = mean_abundance, 
                                   label = name)) +
  scale_size_continuous(range = c(1,4)) +
  scale_fill_identity() +
  # General theme:
  theme_void() + 
  theme(legend.position="none") +
  coord_equal()
ggsave("beaufort_althea/out/y2000.png")

