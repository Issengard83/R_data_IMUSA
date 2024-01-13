### Data analysis for the manuscript:
### Seroprevalence of Leptospira antibodies in dogs and cats attending to municipal 
### spay/neuter campaigns, Santa Fe, Argentina.
### Spatial analysis
### Author: Tamara Ricardo
### Last update:
# Wed Jan 10 16:08:13 2024 ------------------------------


# LOAD PACKAGES -----------------------------------------------------------
# devtools::install_github('Chrisjb/basemapR')

pacman::p_load(
  # Map tools
  sf, 
  tmap,
  tmaptools,
  maptiles,
  basemapR,
  # Analysis
  gtsummary,
  glmmTMB,
  performance,
  spdep,
  # GGplot tools
  GGally,
  # Colorblind-friendly palettes
  scico,
  # Data management
  flextable,
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
                                "Moderate (5-10%)", "High/Very high (>10%)")),
    # # Characters as factor
    # across(.cols = where(is.character), .fns = as.factor)
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


# Aggregate seropositives by census tract ----------------------------------
imu_p_ct <- imu_shp %>% 
  
  ## Summarise data by census tract
  summarise(
    sampled = n(),
    POS = if_else(MAT_res=="POS", 1, 0) %>% sum(na.rm = T),
    .by = c(redcode, adm_dis, pc_incid_h, p_casab_ho, ends_with("_pct"))) %>% 
  
  ## Simplify geometry
  st_cast(to = "POINT") %>% 
  st_make_valid()

# Base map ----------------------------------------------------------------
base_osm <- get_tiles(provider = "OpenStreetMap", 
                      x = st_bbox(var_ct) %>% expand_bbox(X = 1500, Y = 1500),
                      zoom = 11, crop = T)

### Base map template
base_map <- 
  tm_shape(shp = base_osm) +
  tm_rgb() +

  ## North arrow
  tm_compass(size = 1.25, text.size = .9) +
  
  ## Scale bar
  tm_scale_bar() +
  # tm_scalebar() + # v.4.0 

  ## Credits 
  tm_credits("© OpenStreetMap contributors and \n POBLACIONES (CONICET, Universidad Católica Argentina)",
             position = c(.01, 0), size = .45) +
  
  ## Layout
  tm_layout(inner.margins = rep(0, 4), 
            outer.margins = rep(.01, 4),
            legend.outside = T,
            legend.outside.position = "right",
            legend.stack = "vertical",
            legend.title.size = .8,
            legend.text.size = .7,
            panel.label.size = .7)


# Figure 1 ----------------------------------------------------------------
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

### Generate map
fig1 <- base_map +
  
  ## Plot districts
  tm_shape(shp = adm_dis) +
  tm_polygons(col = "adm_dis",
              alpha = .6,
              palette = scico(n = 8, palette = "batlow"), 
              title = "Adm. districts",
              legend.is.portrait = F) +
  
  ## Informal settlements and slums
  tm_shape(shp = inf_shp) +
  tm_fill(col = "#B32F9C", alpha = .9) +
  
  ## Plot samples by origin
  tm_shape(shp = imu_shp) +
  tm_bubbles(col = "animal",
             size = "n",
             palette = "Set2",
             alpha = .8, 
             jitter = .1,
             title.col = "Animal",
             title.size ="Number of samples",
             legend.col.show = F,
             legend.size.show = T,
             legend.col.is.portrait = F,
             legend.size.is.portrait = F) +
  
  ## Facets by origin
  tm_facets(by = "origen_muestra", free.coords = F, ncol = 1) +
  
  ## Plot IMUSA locations
  tm_shape(shp = imu_loc) +
  tm_markers(border.col = "white", 
             border.lwd = 0, 
             icon.scale = 1.1) +
  
  ## Facets by type of marker
  tm_facets(by = "type") +
  tm_layout(legend.outside.position = "bottom", 
            legend.stack = "vertical", 
            panel.labels = c("A", "B"))

 
### Save map
tmap_save(fig1, "FIGS/fig1.png", dpi = 300)


### Clean working environment
rm(imu_loc, imu_shp, base_osm, fig1)

# Poisson analysis --------------------------------------------------------
### Generate data
imu_ct <- imu_p_ct %>% 
  st_drop_geometry() 

### Association with response variable
imu_ct %>% 
  select(POS, adm_dis, pc_incid_h, pc_hog_pct, h_nbi_pct, p_casab_ho, 
         h_hacinami_pct, h_cloaca_pct, h_hoyo_pct, h_agua_viv_pct, h_paviment_pct, h_residuos_pct) %>% 
  tbl_uvregression(y = POS, 
                   method = glmmTMB::glmmTMB,
                   formula = "{y} ~ {x} + (1|adm_dis)",
                   method.args = list(family = poisson),
                   exponentiate = T,
                   label = list(
                     pc_incid_h ~ "Incidence of chronic poverty",
                     pc_hog_pct ~ "Chronic poverty (%)",
                     h_nbi_pct ~ "NBI (%)",
                     p_casab_ho ~ "Poor housing (%)",
                     h_hacinami_pct ~ "Overcrowding (%)",
                     h_cloaca_pct ~ "Absence of sewage (%)",
                     h_hoyo_pct ~ "Pit or cesspool drainage (%)",
                     h_agua_viv_pct ~ "Absence of water pipes (%)",
                     h_paviment_pct ~ "Paved roads (%)",
                     h_residuos_pct ~ "Regular garbage collection (%)")) %>% 
  bold_labels() %>% 
  bold_p() 

### Correlation of numeric variables
ggpairs(data = imu_ct, 
        columns = c("h_nbi_pct", "p_casab_ho", "h_residuos_pct"))

## Univariate models
fit1 = glmmTMB(POS ~ pc_incid_h + (1|adm_dis), family = poisson,
               data = imu_ct)

fit2 = glmmTMB(POS ~ h_nbi_pct + (1|adm_dis), family = poisson,
               data = imu_ct)


fit3 = glmmTMB(POS ~ p_casab_ho + (1|adm_dis), family = poisson,
               data = imu_ct)

fit4 = glmmTMB(POS ~ h_residuos_pct + (1|adm_dis), family = poisson,
               data = imu_ct)

## Multivariate model
fit5 = glmmTMB(POS ~ p_casab_ho + h_residuos_pct + (1|adm_dis), family = poisson,
               data = imu_ct)

### Check performance
compare_performance(fit1, fit2, fit3, fit4, fit5, metrics = "common", rank = T) %>% 
  plot()

### R-squared
r2(fit1)

# ### Table regression coefficients
# tbl_regression(fit1, exp = T) %>% 
#   bold_p() %>% 
#   as_flex_table() %>% 
#     save_as_docx(path = "tab2.docx")

### Check residuals
DHARMa::testResiduals(fit1)

check_overdispersion(fit1)

check_zeroinflation(fit1)

### Clean working environment
rm(fit1, fit2, fit3, fit4, fit5, imu_ct)

# Figure 2 ----------------------------------------------------------------
# ### Update base map template
# base_map <- tm_shape(shp = base_osm) + 
#   tm_rgb() +
#   # tm_basemap(server = "OpenStreetMap", bbox) + # v.4.0
#   ## North arrow
#   tm_compass() +
#   ## Scale bar
#   tm_scale_bar() +
#   # tm_scalebar() + # v.4.0 
#   ## Credits 
#   tm_credits("© OpenStreetMap contributors and \n POBLACIONES (CONICET, Universidad Católica Argentina)",
#              position = c(.01, 0), size = .45) +
#   ## Layout
#   tm_layout(inner.margins = rep(0, 4), 
#             outer.margins = rep(.01, 4),
#             legend.outside = T,
#             legend.outside.position = "right",
#             legend.stack = "vertical",
#             legend.title.size = .8,
#             legend.text.size = .7,
#             panel.label.size = .7)

### Generate map
fig2 <- base_map +
  ## Plot census tracts
  tm_shape(shp = var_ct) +
  tm_polygons(col = c("pc_incid_h", "h_nbi_pct", "p_casab_ho", "h_hacinami_pct"),
          style = "pretty",
          alpha = .9,
          border.alpha = .45,
          palette = scico(n = 6, palette = "lipari", direction = -1),
          title = c("Category", rep("Percentage (%)", 5)),
          legend.is.portrait = T) +
  tm_facets(free.scales = T, ncol = 2) + 
  ## Add admin districts
  tm_shape(shp = adm_dis) +
  tm_borders() +
  ## Informal settlements and slums
  tm_shape(shp = inf_shp) +
  tm_fill(col = "#B32F9C", alpha = .9) +
  ## Add number of seropositive animals per census tract
  tm_shape(shp = imu_p_ct) +
  tm_bubbles(size = "POS",
             col = "#FC8D62",
             alpha = .8,
             style = "cat",
             legend.size.show = F) +
  ## Facet labels
  tm_layout(panel.labels = c("Incidence of chronic poverty",
                             "Unsatisfied basic needs (NBI)",
                             "Poor Housing",
                             "Overcrowding"),
            legend.outside = F,
            legend.position = c("right","top"),
            legend.bg.color = "white",
            legend.bg.alpha = .5)

# ### Save map
# tmap_save(fig2, "FIGS/fig2.png", width = 16, height = 12, units = "cm", dpi = 300)

# Figure 3 ----------------------------------------------------------------
### Generate map
fig3 <- base_map +
  
  ## Plot census tracts
  tm_shape(shp = var_ct) +
  tm_polygons(col = c("h_agua_viv_pct","h_cloaca_pct","h_hoyo_pct",
                  "h_paviment_pct", "h_residuos_pct"),
          style = "pretty",
          alpha = .75,
          border.alpha = .45,
          palette = scico(n = 6, palette = "lipari", direction = -1),
          title = rep("Percentage (%)", 6),
          legend.is.portrait = T) +
  tm_facets(free.scales = T, ncol = 2) + 
  
  ## Add admin districts
  tm_shape(shp = adm_dis) +
  tm_borders() +
  ## Informal settlements and slums
  tm_shape(shp = inf_shp) +
  tm_fill(col = "#B32F9C", alpha = .9) +
  ## Add number of seropositive animals per census tract
  tm_shape(shp = imu_p_ct) +
  tm_bubbles(size = "POS",
             col = "#FC8D62",
             alpha = .8,
             style = "cat",
             legend.size.show = F) +
  
  ## Facet labels
  tm_layout(panel.labels = c("Absence of water pipes",
                             "Absence of sewage",
                             "Pit or cesspool drainage",
                             "Paved roads",
                             "Regular garbage collection"),
            legend.outside = F,
            legend.position = c("right","top"),
            legend.bg.color = "white",
            legend.bg.alpha = .5)


# ### Save map
# tmap_save(fig3, "FIGS/fig3.png", width = 16, height = 16, units = "cm", dpi = 300)

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

# # Plot neighbors
# plot(st_geometry(var_ct), border = "lightgray")
# plot.nb(nb, st_geometry(imu_poly_ct), add = T, col = "purple")

# Neighbors matrix
nbw <- nb2listw(nb, style = "W")

# Global Moran's I (census tract) -----------------------------------------
gmoran <- moran.test(imu_poly_ct$POS, listw = nbw, alternative = "greater")

gmoran$p.value

# Plot results
moran.plot(imu_poly_ct$POS, nbw)

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

# Figure 4 ----------------------------------------------------------------
fig4 <- base_map +
  ## Plot census tracts
  tm_shape(shp = var_ct) +
  tm_fill() +
  ## Plot variables
  tm_shape(shp = imu_poly_Ii) +
  tm_fill(col = c("POS","Ii","Zii_cat","quadrant"),
          # style = "pretty",
          alpha = .9,
          # palette = "viridis",
          palette = scico(n = 5, palette = "lipari", direction = -1),
          legend.is.portrait = T) +
  tm_facets(free.scales = T, ncol = 2) + 
  ## Plot admin districts
  tm_shape(shp = adm_dis) +
  tm_borders() + 
  ## Facet labels
  tm_layout(panel.labels = c("Seropositives",
                             "Local Moran's I",
                             "Spatial autocorrelation",
                             "Spatial clusters"),
            legend.outside = F,
            legend.position = c("right","top"),
            legend.bg.color = "white",
            legend.bg.alpha = .5) +
  ## add borders
    tm_shape(shp = var_ct) +
    tm_borders(alpha = .5)

### Save map
tmap_save(fig4, "FIGS/fig4.png", width = 16, height = 12, units = "cm", dpi = 300)

### Clean working environment
rm(list = setdiff(ls(), c("inf_shp", "var_ct", "imu_poly_Ii", "imu_p_ct")))
