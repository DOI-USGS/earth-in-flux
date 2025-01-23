### Findex Visualizations ###
### Shar Siddiqui         ###
### 1/22/2025             ###

## Load Libraries
library(sf)
library(tidyverse)
library(ggplot2)

sf_use_s2(FALSE) ## for degenerate vertice error message

## Load data
hybas_04_sf <- read_sf('./in/BasinAtlas_v10.gdb','BasinAtlas_v10_lev04')
weights <- read.csv('./in/Final_Weights.csv')
threat_data <- read.csv('./in/Global_Fisheries_Threats.csv') ##Q: in which script where this variables defined for the other script?

## Data Processing 
## grab only columns with '.LS' and HYBAS_ID
threat_data <- threat_data[,grep('_ID|_LS',names(threat_data))]
## transform from short to long
threat_data <- threat_data |> pivot_longer(!HYBAS_ID,names_to = 'Threat_Code',values_to = 'value') 
## append column with long form name and weight information
threat_data$Threat_Code <- gsub("_LS","",threat_data$Threat_Code) 
threat_data <- left_join(threat_data,weights, by = "Threat_Code")

## append threat data to basin df
hybas_04_threats <- left_join(hybas_04_sf,threat_data,by="HYBAS_ID")

## Plot 1: Visualize Threats
ggplot(threat_data)+
  geom_boxplot(mapping=aes(x=Threat_Code,y=value)) +
  xlab('Threat') + ylab('Score') +
  theme_bw()

## Plot 2: Replicate Map
summary <- hybas_04_threats %>%
  group_by(HYBAS_ID) %>%
  summarise(meanTS = mean(value,na.rm=T)) %>%
  ungroup()
map <- ggplot() +
  geom_sf(summary,mapping=aes(fill=meanTS),color='transparent') +
  scale_fill_gradient(low="tan", high="tan4")+
  labs(fill='Mean Threat Score')+
  theme_bw()

ggsave('map.png',map,width=8.27,height = 3.44,unit='in')

## Plot 4: Plot basins by number of threats 

## "We find that nearly 98% of all watersheds globally face multiple 
## threats, and nearly half face more than 16 different threats, likely 
## acting synergistically." - Abstract Q: Can't find how this was defined

## Going to define facing a threat as being 0

## Plot 4a: Map of basins filled by number of threats
## Add column 'N_threats'
num_threats <- threat_data %>%
  filter(value > 0) %>%
  mutate(threat_present = 1) %>%
  group_by(HYBAS_ID) %>%
  summarise(N_threats = sum(threat_present)) %>% ## calculate number of threats not equal to zero
  ungroup()

hybas_04_sf <- left_join(hybas_04_sf,num_threats %>% select('HYBAS_ID','N_threats'),relationship='one-to-one')
hybas_04_sf$N_threats[which(is.na(hybas_04_sf$N_threats))] <- 0
n_map <- ggplot() +
  geom_sf(hybas_04_sf,mapping=aes(fill=N_threats),color='transparent') +
  scale_fill_gradient(low="white", high="darkblue")+
  labs(fill='Number of Threats')+
  theme_bw()

ggsave('number_threats_map.png',n_map,width=8.27,height = 3.44,unit='in')

## Plot 4b: Grid of number of threats by continent/habitat type
## Note: First two digits of HYBAS ID denote region

hybas_04_sf$cont_code <- substr(hybas_04_sf$HYBAS_ID,1,2)
hybas_04_sf$Continent <- NA
hybas_04_sf$Continent[hybas_04_sf$cont_code == 10] <- 'Africa'
hybas_04_sf$Continent[hybas_04_sf$cont_code == 20] <- 'Europe'
hybas_04_sf$Continent[hybas_04_sf$cont_code == 30] <- 'Asia'
hybas_04_sf$Continent[hybas_04_sf$cont_code == 40] <- 'Asia'
hybas_04_sf$Continent[hybas_04_sf$cont_code == 50] <- 'Oceania'
hybas_04_sf$Continent[hybas_04_sf$cont_code == 60] <- 'South America'
hybas_04_sf$Continent[hybas_04_sf$cont_code == 70] <- 'North America'
hybas_04_sf$Continent[hybas_04_sf$cont_code == 80] <- 'North America'
hybas_04_sf$Continent[hybas_04_sf$cont_code == 90] <- 'North America'

install.packages("ggridges")
library(ggridges)
ggplot() +
  geom_density_ridges2(hybas_04_sf, mapping=aes(x=N_threats,y=Continent),
                       fill = "lightblue", alpha = 0.5) +
  xlab('Distribution of Number of Threats per Basin') +
  theme_bw()

## Plot 4c: Relationship between Number of Threats and Mean Threat Score
summary <- left_join(summary, num_threats %>% select('HYBAS_ID','N_threats'),relationship='one-to-one')
ggplot() +
  geom_point(summary,mapping=aes(x=N_threats,y=meanTS))+
  geom_smooth(summary,mapping=aes(x=N_threats,y=meanTS))+
  xlab('Number of Threats') + ylab('Mean Threat Score')+
  theme_bw() #(For each basin, across all Threat Types)
cor(summary$meanTS,summary$N_threats,use='complete.obs')^2

ggplot() +
  geom_boxplot(summary,mapping=aes(x=N_threats, group=as.character(N_threats),y=meanTS,fill=N_threats))+
  scale_fill_gradient(low="white", high="darkblue")+
  xlab('Number of Threats') + ylab('Mean Threat Score')+
  theme_bw() #(For each basin, across all Threat Types)

table(hybas_04_sf$N_threats)
1342-23

## What are the most common co-occurring threats?

