# clean files
# btw, I used chatGPT a lot to write this code
rm(list=ls())

# libraries
library(sf)
library(ggplot2)
library(dplyr)
library(readr)
library(rmapshaper)

# bring csv
data <- read_csv("../output/results_inversion.csv", col_names=FALSE)
colnames(data) <- c("CD_MICRO", "Tbar", "Bbar")

# bring shapefile
shapes <- st_read("../input/shapefile/BR_Microrregioes_2021.shp")
shapes <- ms_simplify(shapes, keep = 0.05)

# ensure both datasets are in characters
shapes$CD_MICRO <- as.character(shapes$CD_MICRO)
data$CD_MICRO <- as.character(data$CD_MICRO)

# join
map_data <- shapes %>% left_join(data, by = "CD_MICRO")

###### TBAR
# Create quintile bins
map_data <- map_data %>%
  mutate(Tbar_quintile = cut(Tbar, 
                             breaks = quantile(Tbar, probs = seq(0, 1, 0.2), na.rm = TRUE), 
                             include.lowest = TRUE, 
                             labels = c("Q1", "Q2", "Q3", "Q4", "Q5")))

# Plot with discrete fill
ggplot() +   
  geom_sf(data = map_data, aes(fill = Tbar_quintile), color = "black", size = 0.01) +   
  scale_fill_viridis_d(option = "plasma", na.value = "white") +  # Discrete scale
  theme_minimal() +   
  labs(fill = "Tbar Quintile", title = "Productivity Shifters")
ggsave("../output/microregion_map_Tbar.png", width = 8, height = 6, dpi = 300)

###### BBAR
# Create quintile bins
map_data <- map_data %>%
  mutate(Bbar_quintile = cut(Bbar, 
                             breaks = quantile(Bbar, probs = seq(0, 1, 0.2), na.rm = TRUE), 
                             include.lowest = TRUE, 
                             labels = c("Q1", "Q2", "Q3", "Q4", "Q5")))

# Plot with discrete fill
ggplot() +   
  geom_sf(data = map_data, aes(fill = Bbar_quintile), color = "black", size = 0.01) +   
  scale_fill_viridis_d(option = "plasma", na.value = "white") +  # Discrete scale
  theme_minimal() +   
  labs(fill = "Bbar Quintile", title = "Amenity Shifters")
ggsave("../output/microregion_map_Bbar.png", width = 8, height = 6, dpi = 300)