### Interactive map for the manuscript:
### Seroprevalence of Leptospira antibodies in dogs and cats attending to municipal 
### spay/neuter campaigns, Santa Fe, Argentina.
### Spatial analysis
### Author: Tamara Ricardo
### Last update:
# Fri Jan 12 14:55:41 2024 ------------------------------


# LOAD PACKAGES -----------------------------------------------------------
pacman::p_load(
  # Map tools
  sf, 
  tmap,
  tmaptools, 
  maptiles,
  # Analysis
  spdep,
  # GGplot tools
  GGally,
  # Colorblind-friendly palettes
  scico,
  # Data management
  janitor,
  tidyverse)


# Load and edit spatial layers --------------------------------------------
### Location of informal settlements and slums
inf_shp <- st_read("SHP/raw/Registro Nacional de Barrios Populares - Santa Fe.shp", 
                   stringsAsFactors = T) %>% 
  
  ## Clean variable names
  clean_names() %>% 
  
  ## Validate geometry
  st_make_valid()


### IMUSA markers
imu_loc <- st_read("SHP/IMUSA_SF.geojson")


### Socioeconomic indicators by census tract
var_ct <- st_read("SHP/VAR_CT_SF.geojson") %>% 
  
  ## Transform variables
  mutate(
    # Data to percentages
    across(.cols = c("h_nbi", "pc_hog", "h_hacinami", "h_agua_red","h_agua_viv", 
                     "h_cloaca", "h_hoyo", "h_paviment", "h_residuos"), 
           .fns = list(pct = ~ round(.x*100/n_hogares, 2))),
    
    # Categorize incidence of chronic poverty
    pc_incid_h = cut(pc_hog_pct, breaks = c(-Inf, 1, 5, 10, Inf),
                     labels = c("Very low (<1%)", "Low (1-5%)", 
                                "Moderate (5-10%)", "High/Very high (>10%)"))
  )


### Sampling data (points)
imu_shp <- readxl::read_excel("data_IMUSA_clean.xlsx") %>% 
  
  ## Relevel origin of samples
  mutate(origen_muestra = if_else(grepl("Móvil", origen_muestra), 
                                  "IMUSA mobile truck location", 
                                  "IMUSA fixed units")) %>%
  
  ## Create spatial points by census tract
  select(id_muestra, origen_muestra, animal, MAT_res, lon, lat) %>% 
  
  ## Dataframe to shapefile
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% 
  
  ## Join socioeconomic indicators (census tract)
  st_join(var_ct, join = st_nearest_feature) %>% 
  
  ## Drop unused factor levels
  mutate(across(.cols = where(is.factor), .fns = fct_drop))


# Generate new spatial layers ---------------------------------------------
### Aggregate seropositives by census tract
imu_p_ct <- imu_shp %>% 
  
  ## Summarise data by census tract
  group_by(across(c(redcode, adm_dis, pc_incid_h, p_casab_ho, 
                  h_nbi_pct:h_residuos_pct))) %>%
  summarise(
    sampled = n(),
    POS = sum(MAT_res, na.rm = T)) %>% 
  
  ## Simplify geometry
  st_cast(to = "POINT") %>% 
  st_make_valid()


### Aggregate samples by animal, origin and census tract
imu_shp <- imu_shp %>% 
  count(redcode, origen_muestra, animal) %>% 
  ## Validate geometry
  st_make_valid()


### Create layer for admin districts
adm_dis <- var_ct %>%
  ## Group data by administrative district
  group_by(adm_dis) %>% 
  summarise(n_hogares = sum(n_hogares, na.rm = F)) %>%
  
  ## Validate geometry
  st_make_valid()


# Spatial analysis (census tract) -----------------------------------------
### Create spatial polygons
imu_poly_ct <- var_ct %>% 
  inner_join(imu_p_ct %>% st_drop_geometry() %>% 
               select(redcode, sampled, POS)) %>% 
  # Replace missing values with 0
  mutate_at(c("sampled","POS"), .funs = ~ replace_na(.x, 0))

### Nearest neighbor analysis
nb <- st_centroid(imu_poly_ct) %>% 
  knearneigh(k = 3) %>% 
  knn2nb()

# Neighbors matrix
nbw <- nb2listw(nb, style = "W")

# Local Moran's I (census tract) ------------------------------------------
lmoran <- localmoran(imu_poly_ct$POS, listw = nbw, alternative = "two.sided")

### Identify clusters
mp <- imu_poly_ct$POS %>% 
  scale() %>% 
  as.vector() %>% 
  moran.plot(., nbw) 

### Add Moran's data to shapefile
imu_poly_Ii <- imu_poly_ct %>% 
  ## Local Moran's I
  bind_cols(lmoran) %>% 
  ## Cluster data
  bind_cols(mp %>% select(x, wx)) %>% 
  # Clean variable names
  rename(lmp = "Pr(z != E(Ii))") %>% 
  # Create new variables
  mutate(
    # Categorized Z-score
    Zii_cat = cut(Z.Ii, breaks = c(-Inf, -1.96, 1.96, Inf),
                  labels = c("Negative SAC", "No SAC", "Positive SAC")),
    # Cluster quadrants
    quadrant = case_when(
      (x>=0 & wx>=0) & lmp<=.05 ~ "High-High",
      (x<=0 & wx<=0) & lmp<=.05 ~ "Low-Low",
      (x>=0 & wx<=0) & lmp<=.05 ~ "High-Low",
      (x<=0 & wx>=0) & lmp<=.05 ~ "Low-High",
      lmp>.05 ~ "Non-significant")
  )

# Generate map ------------------------------------------------------------
tmap_mode("view")

### Interactive map
int_map <- tm_basemap(server = "OpenStreetMap") +
  
  ## Plot socioeconomic variables by census tract
  # Unsatisfied basic needs
  tm_shape(shp = var_ct) +
  tm_polygons(col = "h_nbi_pct",
              style = "pretty",
              alpha = .9,
              border.alpha = .45,
              palette = scico(n = 6, palette = "lipari", direction = -1),
              title = "Unsatisfied basic needs",
              popup.vars = F,
              group = "Unsatisfied basic needs") +
  
  # Lack of piped water
  tm_shape(shp = var_ct) +
  tm_polygons(col = "h_agua_viv_pct",
              style = "pretty",
              alpha = .9,
              border.alpha = .45,
              palette = scico(n = 6, palette = "lipari", direction = -1),
              title = "Absence of water pipes (%)",
              popup.vars = F,
              group = "Water pipes") +
  
  # Lack of sewage system
  tm_shape(shp = var_ct) +
  tm_polygons(col = "h_cloaca_pct",
              style = "pretty",
              alpha = .9,
              border.alpha = .45,
              palette = scico(n = 6, palette = "lipari", direction = -1),
              title = "Absence of sewage system (%)",
              popup.vars = F,
              group = "Sewage system") +
  
  # Regular garbage collection
  tm_shape(shp = var_ct) +
  tm_polygons(col = "h_residuos_pct",
              style = "pretty",
              alpha = .9,
              border.alpha = .45,
              palette = scico(n = 6, palette = "lipari", direction = -1),
              title = "Regular garbage collection (%)",
              popup.vars = F,
              group = "Garbage collection") +
  
  # Paved roads
  tm_shape(shp = var_ct) +
  tm_polygons(col = "h_paviment_pct",
              style = "pretty",
              alpha = .9,
              border.alpha = .45,
              palette = scico(n = 6, palette = "lipari", direction = -1),
              title = "Paved roads (%)",
              popup.vars = F,
              group = "Paved roads") +
  
  # Incidence of chronic poverty
  tm_shape(shp = var_ct) +
  tm_polygons(col = "pc_incid_h",
              style = "pretty",
              alpha = .9,
              border.alpha = .45,
              palette = scico(n = 6, palette = "lipari", direction = -1),
              title = "Incidence of chronic poverty",
              popup.vars = c("Adm. district" = "adm_dis", 
                             "Number of households" = "n_hogares",
                             "Overcrowding (%)" = "h_hacinami_pct", 
                             "Absence of water pipes (%)" = "h_agua_viv_pct", 
                             "Absence of sewage (%)" = "h_cloaca_pct", 
                             "Pit or cesspool drainage (%)" = "h_hoyo_pct",
                             "Paved roads (%)" = "h_paviment_pct",
                             "Regular garbage collection (%)" = "h_residuos_pct",
                             "Unsatisfied Basic Needs (%)" = "h_nbi_pct"),
              id = "redcode",
              group = "Incidence of chronic poverty") +
  
  ## Plot variables from spatial analysis
  # Spatial clusters
  tm_shape(shp = imu_poly_Ii) +
  tm_polygons(col = "quadrant",
              style = "pretty",
              alpha = .9,
              palette = scico(n = 5, palette = "acton", direction = 1),
              title = "Spatial clusters",
              popup.vars = c("Local Moran's I" = "Ii",
                             "Spatial autocorrelation" = "Zii_cat",
                             "quadrant", 
                             "Seropositive animals" = "POS"),
              id = "redcode",
              group = "Spatial clusters") +
  
  ## Add number of sampled animals per species
  tm_shape(shp = imu_shp) +
  tm_bubbles(col = "animal",
             size = "n",
             palette = "Set2",
             alpha = .8, 
             jitter = .1,
             popup.vars = F,
             group = "Sampled",
             labels = c("Cat", "Dog")) +

  ## Add number of seropositive animals per census tract
  tm_shape(shp = imu_p_ct) +
  tm_bubbles(size = "POS",
             col = "#FC8D62",
             alpha = .8,
             style = "cat",
             group = "Seropositives") +
  
  ## Informal settlements and slums
  tm_shape(shp = inf_shp) +
  tm_fill(col = "#B32F9C", 
          alpha = .9,
          popup.vars = c("Name" = "nombre_bar"),
          group = "Informal settlements/slums") +
  
  ## Plot IMUSA locations
  tm_shape(shp = imu_loc) +
  tm_markers(border.col = "white",
             border.lwd = 0, 
             icon.scale = 1.1, 
             popup.vars = F,
             id = "type",
             group = "IMUSA markers")


# Configure and export map ------------------------------------------------
int_map1 <- int_map %>% tmap_leaflet() %>% 
  ## Hide layers
  leaflet::hideGroup(group = c("Unsatisfied basic needs", 
                               "Water pipes", 
                               "Sewage system",
                               "Garbage collection", 
                               "Paved roads", 
                               "Sampled"))
## Export
library(htmlwidgets) 

# Save the map as an HTML widget 
saveWidget(int_map1, file = "FileS4_interactive_map.html") 
