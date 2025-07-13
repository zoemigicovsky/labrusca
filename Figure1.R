#Figure 1 

library(tidyverse)
library(CoordinateCleaner) 
library(cowplot)
library(rnaturalearth)
library(terra) 
library(geodata) 

# Load in supplemental table with coordinate information for all samples 
TableS1 <- read_delim("TableS1.csv")

# Filter to only include geo-referenced accessions
occ_wild_conserved <- TableS1 %>%
  filter(!is.na(LON),!is.na(LON)) %>% 
  mutate(long = LON, lat = LAT, Source, .keep="none")

###### Figure 1A ####### 
# To generate Figure 1A requires first downloading the occurrence data from GBIF for Vitis labrusca. This step is included below but commented out as the occurrence data are also posted to the repository 'occurrence_labrusca_dat.csv' 

# The information associated with this download is: GBIF.org (13 July 2025) GBIF Occurrence Download https://doi.org/10.15468/dl.xu3a52


# # GBIF user info
# user= #goes here
# pwd= #goes here
# email= #goes here 
# 
# # Vitis labrusca
# taxonKey <- 5372424
# basisOfRecord <- c('PRESERVED_SPECIMEN', 'HUMAN_OBSERVATION')
# hasCoordinates <- TRUE # limit to records with coordinates
# country_codes <- c("CA", "US") # limit to Canada and USA
# 
# # Download data
# down_code = occ_download(
#   pred("taxonKey", taxonKey),
#   pred_in("basisOfRecord", basisOfRecord),
#   pred("hasCoordinate", hasCoordinates),
#   pred_in("country", country_codes),
#   format = "DWCA",
#   user=user, pwd=pwd, email=email)
# 
# # Wait for download to finish
# occ_download_wait(down_code)
# 
# download_labrusca <- occ_download_get(down_code[1], overwrite = TRUE)

# Extract and save as csv
# Only load in columns of interest and ensure the character type is correct

gbif_labrusca <- read_csv("occurrence_labrusca_dat.csv", col_types = cols_only(species = col_character(), countryCode = col_character(), decimalLatitude = col_double(), decimalLongitude = col_double(), coordinateUncertaintyInMeters = col_double(), year = col_double(), basisOfRecord = col_character(), gbifID = col_double(), stateProvince = col_character(), hasGeospatialIssues = col_character()))

# Define allowed locations, all lowercase, based on FNA

allowed_states <- tolower(c(
  "Ontario", "Ont.",
  "Alabama", "Ala.",
  "Connecticut", "Conn.",
  "Delaware", "Del.",
  "District of Columbia", "D.C.",
  "Georgia", "Ga.",
  "Illinois", "Ill.",
  "Indiana", "Ind.",
  "Kentucky", "Ky.",
  "Maine",
  "Maryland", "Md.",
  "Massachusetts", "Mass.",
  "Michigan", "Mich.",
  "Mississippi", "Miss.",
  "New Hampshire", "N.H.",
  "New Jersey", "N.J.",
  "New York", "N.Y.",
  "North Carolina", "N.C.",
  "Ohio",
  "Pennsylvania", "Pa.",
  "Rhode Island", "R.I.",
  "South Carolina", "S.C.",
  "Tennessee", "Tenn.",
  "Vermont", "Vt.",
  "Virginia", "Va.",
  "Wisconsin", "Wis."
))

# Standardize the stateProvince column to lowercase and filter
filtered_data <- gbif_labrusca %>%
  mutate(stateProvince_lower = tolower(stateProvince)) %>%
  filter(stateProvince_lower %in% allowed_states)

occ_lab_clean <- filtered_data %>% 
  filter(hasGeospatialIssues == FALSE) %>% # remove records with geospatial issues
  dplyr::select(species, countryCode, decimalLatitude, decimalLongitude, coordinateUncertaintyInMeters, year, basisOfRecord, gbifID) %>% #grab necessary columns 
  filter(!is.na(decimalLongitude)) %>% # should have coords but you never know
  filter(decimalLongitude != 0) %>% 
  filter(coordinateUncertaintyInMeters < 30000 | is.na(coordinateUncertaintyInMeters)) %>%
  CoordinateCleaner::cc_cen(buffer = 2000) %>%
  CoordinateCleaner::cc_inst(buffer = 200) %>%
  CoordinateCleaner::cc_sea()%>% 
  distinct(decimalLatitude, decimalLongitude, gbifID, .keep_all = T) 

# combine gbif records with accessions sampled 
occ_gbif <- occ_lab_clean %>%
  transmute(long = decimalLongitude,
            lat = decimalLatitude,
            Source = "GBIF")

all_points <- bind_rows(occ_gbif, occ_wild_conserved)

# Re order levels so that GBIF plots first 

all_points <- all_points %>%
  mutate(Source = fct_relevel(Source, "GBIF", "Conserved", "Wild"))

# Define colors, shapes, sizes by Source
fill_vals <- c("GBIF" = "red", "Wild" = "white", "Conserved" = "black")
color_vals <- c("GBIF" = "red", "Wild" = "black", "Conserved" = "black")  
shape_vals <- c("GBIF" = 21, "Wild" = 21, "Conserved" = 24)
size_vals <- c("GBIF" = 2, "Wild" = 3, "Conserved" = 3)

# Load world land for background
land <- ne_countries(scale = "medium", returnclass = "sf")
lakes <- ne_download(scale = "medium", type = "lakes", category = "physical", returnclass = "sf")
usa_states <- map_data("state")

# USDA in Geneva coordinates
geneva_coords <- data.frame(
  long = -77.2198236,
  lat = 42.8637196
)

# Calculate plotting extent with margin
margin <- 2
lon_range <- range(all_points$long, na.rm = TRUE) + c(-margin, margin)
lat_range <- range(all_points$lat, na.rm = TRUE) + c(-margin, margin)

occ_plot <- ggplot() +
  # Background layers
  geom_sf(data = land, fill = "grey90", color = NA) +
  geom_sf(data = lakes, fill = "white", color = NA) +
  geom_path(data = usa_states, aes(x = long, y = lat, group = group),
            color = "grey60", linewidth = 0.3, alpha = 0.8) +
  
  # Points layer with mapped aesthetics
  geom_point(data = all_points,
             aes(x = long, y = lat,
                 fill = Source,
                 shape = Source,
                 size = Source,
                 color = Source),
             stroke = 0.6,
             alpha = 0.6) +
  
  # Geneva star
  geom_point(data = geneva_coords,
             aes(x = long, y = lat),
             shape = 23,
             size = 5.5,
             fill = "yellow",
             color = "black",
             stroke = 1) +
  
  # Scales for legend and aesthetics
  scale_fill_manual(name = "Source", values = fill_vals) +
  scale_color_manual(name = "Source", values = color_vals) +
  scale_shape_manual(name = "Source", values = shape_vals) +
  scale_size_manual(name = "Source", values = size_vals) +
  
  # Coordinate limits
  coord_sf(xlim = lon_range, ylim = lat_range, expand = FALSE) +
  
  # Labels and theme
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal(base_size = 11) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    legend.key = element_blank()
  )

occ_plot


###### Figure 1B ####### 

# Load the world map
climate <- geodata::worldclim_global(var = "bio", res = 2.5, path = tempdir())
bio1 <- climate[[1]]  # Annual mean temperature
ne_extent <- ext(-87, -66, 33, 45)
bio1_crop <- crop(bio1, ne_extent)
bio1_df <- as.data.frame(bio1_crop, xy = TRUE)
colnames(bio1_df)[3] <- "bio1"

# Filter for eastern USA map
north_east_usa_map <- map_data("state") %>%
  filter(long >= -87 & long <= -66,
         lat >= 33 & lat <= 45)

# Coordinates for USDA Genebank in Geneva, NY
geneva_coords <- data.frame(
  long = -77.2198236,
  lat = 42.8637196
)

# Load lakes as sf polygons
lakes <- ne_download(scale = "medium", type = "lakes", category = "physical", returnclass = "sf")

# US state borders
usa_states <- map_data("state") %>%
  filter(long >= -87 & long <= -66,
         lat >= 33 & lat <= 45)

# Geneva coordinates
geneva_coords <- data.frame(
  long = -77.2198236,
  lat = 42.8637196
)

# Longitude and latitude limits based on raster extent + margin
lon_range <- range(bio1_df$x, na.rm = TRUE) + c(-0.5, 0.5)
lat_range <- range(bio1_df$y, na.rm = TRUE) + c(-0.5, 0.5)

temp_plot <- ggplot() +
  # Temperature raster base
  geom_tile(data = bio1_df, aes(x = x, y = y, fill = bio1), alpha = 1) +
  scale_fill_viridis_c(option = "B", name = "Annual Mean\nTemperature\n(°C)") +
  
  ggnewscale::new_scale_fill() +
  
  # White lakes on top
  geom_sf(data = lakes, fill = "white", color = NA) +
  
  # US state borders
  geom_path(data = usa_states,
            aes(x = long, y = lat, group = group),
            color = "grey40", linewidth = 0.3, alpha = 0.6) +
  
  # Points with population type fill/shape
  geom_point(data = occ_wild_conserved,
             aes(x = long, y = lat, fill = Source, shape = Source),
             color = "black",
             stroke = 0.6,
             size = 2.5,
             alpha = 0.7) +
  scale_shape_manual(name = "Source", values = c("Wild" = 21, "Conserved" = 24)) +
  scale_fill_manual(name = "Source", values = c("Wild" = "white", "Conserved" = "black")) +
  
  # Geneva star
  geom_point(data = geneva_coords,
             aes(x = long, y = lat),
             shape = 23,
             size = 5.5,
             fill = "yellow",
             color = "black",
             stroke = 1) +
  
  # Coordinate limits
  coord_sf(xlim = lon_range, ylim = lat_range, expand = FALSE) +
  
  # Theme and labels
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal(base_size = 11) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8)
  )

temp_plot

###### Figure 1C ####### 

# Load the world map
climate <- geodata::worldclim_global(var = "bio", res = 2.5, path = tempdir())
bio12 <- climate[[12]]  # Annual precipitation (mm)
ne_extent <- ext(-87, -66, 33, 45)
bio12_crop <- crop(bio12, ne_extent)
bio12_df <- as.data.frame(bio12_crop, xy = TRUE)
colnames(bio12_df)[3] <- "bio12"

# Convert precipitation values from mm to cm
bio12_df$bio12_cm <- bio12_df$bio12 / 10

# Filter for eastern USA map
north_east_usa_map <- map_data("state") %>%
  filter(long >= -87 & long <= -66,
         lat >= 33 & lat <= 45)

# Coordinates for USDA Genebank in Geneva, NY
geneva_coords <- data.frame(
  long = -77.2198236,
  lat = 42.8637196
)

# Load lakes as sf polygons
lakes <- ne_download(scale = "medium", type = "lakes", category = "physical", returnclass = "sf")

# US state borders filtered to your region
usa_states <- map_data("state") %>%
  filter(long >= -87 & long <= -66,
         lat >= 33 & lat <= 45)

# Longitude and latitude limits based on raster extent + margin
lon_range <- range(bio12_df$x, na.rm = TRUE) + c(-0.5, 0.5)
lat_range <- range(bio12_df$y, na.rm = TRUE) + c(-0.5, 0.5)

precip_plot <- ggplot() +
  # Precipitation raster base (converted to cm)
  geom_tile(data = bio12_df, aes(x = x, y = y, fill = bio12_cm), alpha = 1) +
  scale_fill_viridis_c(name = "Annual\nPrecipitation\n(cm)", direction = -1) +
  
  ggnewscale::new_scale_fill() +
  
  # White lakes on top
  geom_sf(data = lakes, fill = "white", color = NA) +
  
  # US state borders
  geom_path(data = usa_states,
            aes(x = long, y = lat, group = group),
            color = "grey40", linewidth = 0.3, alpha = 0.6) +
  
  # Points with population type fill/shape
  geom_point(data = occ_wild_conserved,
             aes(x = long, y = lat, fill = Source, shape = Source),
             color = "black",
             stroke = 0.6,
             size = 2.5,
             alpha = 0.7) +
  scale_shape_manual(name = "Source", values = c("Wild" = 21, "Conserved" = 24)) +
  scale_fill_manual(name = "Source", values = c("Wild" = "white", "Conserved" = "black")) +
  
  # Geneva star
  geom_point(data = geneva_coords,
             aes(x = long, y = lat),
             shape = 23,
             size = 5.5,
             fill = "yellow",
             color = "black",
             stroke = 1) +
  
  # Coordinate limits
  coord_sf(xlim = lon_range, ylim = lat_range, expand = FALSE) +
  
  # Theme and labels
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal(base_size = 11) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8)
  )

precip_plot

###### Make Figure 1B/1C boxplots ####### 

# Convert to spatial points
coords <- occ_wild_conserved %>% dplyr::select(long, lat) 
points_sp <- terra::vect(coords, geom = c("long", "lat"), crs = crs(bio1))

# Extract temperature and precipitation values
occ_wild_conserved$bio1 <- terra::extract(bio1, points_sp)[,2] 
occ_wild_conserved$bio12 <- terra::extract(bio12, points_sp)[,2]/10

# Temperature boxplot
temp_box <- ggplot(occ_wild_conserved, aes(x = Source, y = bio1, fill = Source)) +
  geom_jitter(aes(x = Source, y = bio1), width = 0.15, color = "black", size = 0.8) +
  geom_boxplot(alpha = 0.7, color = "black", width = 0.4, outlier.shape = NA) +  # hides default outliers
  scale_fill_manual(values = c("Wild" = "white", "Conserved" = "darkgrey")) +
  labs(y = "Temperature (°C)", x = NULL) +
  theme_bw() +
  theme(plot.background = element_blank(),    # Remove the background around the plot (outside the panel)
        panel.background = element_rect(fill = "white", color = "black"), 
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.position = "none",
        text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        legend.title = element_text(size = 7, color = "black", face="bold"),
        legend.text = element_text(size = 5.5, color = "black"),
        panel.grid = element_blank())

# Precipitation boxplot
precip_box <- ggplot(occ_wild_conserved, aes(x = Source, y = bio12, fill = Source)) +
  geom_jitter(aes(x = Source, y = bio12), width = 0.15, color = "black", size = 0.8) +
  geom_boxplot(alpha = 0.7, color = "black", width = 0.4) +  # narrower
  scale_fill_manual(values = c("Wild" = "white", "Conserved" = "darkgrey")) +
  labs(y = "Precipitation (cm)",x = NULL) +
  theme_bw() +
  theme(plot.background = element_blank(),    # Remove the background around the plot (outside the panel)
        panel.background = element_rect(fill = "white", color = "black"), 
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.position = "none",
        text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        legend.title = element_text(size = 7, color = "black", face="bold"),
        legend.text = element_text(size = 5.5, color = "black"),
        panel.grid = element_blank())

###### Figure 1 (Merged) ####### 

library(cowplot)

# Adjusted boxplot positions for vertical layout (lower in each panel)
temp_combined <- ggdraw() +
  draw_plot(temp_plot, 0, 0, 1, 1) +
  draw_plot(temp_box, 0.51, 0.12, 0.27, 0.47)  

precip_combined <- ggdraw() +
  draw_plot(precip_plot, 0, 0, 1, 1) +
  draw_plot(precip_box, 0.51, 0.12, 0.27, 0.47) 

# Stack the plots vertically
full_figure <- plot_grid(occ_plot, temp_combined, precip_combined, nrow = 3,labels = c("A", "B", "C"))

# Save the figure
ggsave("Figure1.pdf", full_figure, width = 7, height = 13.5)

###### Figure 1 Stats ####### 

# What proportion of the temp and precip variables did conserved vs wild cover?

# Calculate total ranges
temp_range <- range(occ_wild_conserved$bio1, na.rm = TRUE)
precip_range <- range(occ_wild_conserved$bio12, na.rm = TRUE)

# Calculate proportion of total range captured by each group
occ_wild_conserved %>%
  group_by(Source) %>%
  summarise(
    temp_min = min(bio1, na.rm = TRUE),
    temp_max = max(bio1, na.rm = TRUE),
    temp_prop = (temp_max - temp_min) / (temp_range[2] - temp_range[1]),
    precip_min = min(bio12, na.rm = TRUE),
    precip_max = max(bio12, na.rm = TRUE),
    precip_prop = (precip_max - precip_min) / (precip_range[2] - precip_range[1])
  )
