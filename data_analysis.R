### Data analysis for the manuscript:
### Seroprevalence of Leptospira antibodies in dogs and cats attending to municipal 
### spay/neuter campaigns, Santa Fe, Argentina.
### Logistic regression analysis
### Author: Tamara Ricardo
### Last update:
# Tue Jan  2 13:10:17 2024 ------------------------------


# LOAD PACKAGES -----------------------------------------------------------
pacman::p_load(
  # Descriptive statistics
  gtsummary,
  # Model fit
  glmmTMB,
  # LASSO
  glmnet,
  caret,
  # Residuals and inference
  performance,
  DHARMa,
  # Data management
  rio,
  janitor,
  flextable,
  tidyverse)

# Load and clean data -----------------------------------------------------
imu_clean <- import("data_IMUSA_clean.xlsx") %>% 
  
  ### Modify variables
  mutate(
    edad_cat = fct_relevel(edad_cat, "cachorro", after = 0),
    
    bcs = as.numeric(bcs),
    
    donde_duerme_limp = fct_lump_min(donde_duerme_limp,  min = 20, 
                                     other_level = "ocasionalmente") %>% 
      fct_relevel("diariamente", after = 0),
    
    vio_roedores_frec = if_else(is.na(vio_roedores_frec) & vio_roedores=="no",
                                "nunca", vio_roedores_frec),
    
    across(where(is.character), .fns = ~ as.factor(.x))) %>% 
  
  ### Create variable for overall number of pets in the household
  rowwise() %>% 
  mutate(n_pets = sum(c_across(n_perros:n_gatos), na.rm = T))


# Create dataset for dog samples ------------------------------------------
imu_dog <- imu_clean %>% 
  ### Select only dog samples
  filter(animal=="perro") %>%
  
  ### Select relevant variables
  select(MAT_res, adm_dis, edad_yrs, edad_cat, sexo, raza_cat, bcs_cat,
         vac_alguna, animal_fun:n_gatos, n_pets, enferm_reciente: donde_duerme_limp) %>% 
  
  ### Remove NAs
  drop_na()


# Explore data ------------------------------------------------------------
### Number of samples per administrative district
imu_clean %>%
  tbl_summary(include = adm_dis, by = MAT_res,
              sort = list(everything() ~ "frequency")) %>% 
  add_p()

# Table 1: descriptive statistics -----------------------------------------
tab1 <- imu_clean %>% 
  select(animal:donde_duerme_limp, n_pets, -bcs) %>% 
  tbl_summary(by = animal,
              missing = "no",
              label = list(
                edad_yrs = "age (years)",
                edad_cat = "age group",
                tuvo_cria = "pregnancies",
                tuvo_abortos = "abortions/stillbirths",
                raza_cat = "breed",
                bcs_cat = "BCS",
                vac_alguna = "vaccines (any)",
                vac_rabia = "vaccines (rabies)",
                vac_rabia_st = "vaccines (rabies at sterilization)",
                vac_sext_canina = "vaccines (hexavalent)",
                vac_ns = "vaccines (DK/DR)",
                animal_fun = "role of the animal",
                animal_alim = "feeding",
                sale_calle = "street access",
                sale_suelto = "unsupervised street access",
                con_agua_barro = "contact w. water/mud",
                con_basurales = "contact w. garbage dumps",
                con_perros = "contact w. dogs",
                con_gatos = "contact w. cats",
                con_anim_otros = "contact w. other animals",
                con_anim_propios = "contact w. household pets",
                con_anim_vecinos = "contact w. neighbor pets",
                con_anim_callejeros = "contact w. stray pets",
                n_perros = "number of dogs",
                n_gatos = "number of cats",
                n_pets = "number of dogs and cats",
                enferm_reciente = "recent illness",
                caza_animales = "hunts animals (any)",
                caza_roedores = "hunts rodents",
                caza_silvestres = "hunts wild animals",
                lugar_mascota = "housing",
                donde_duerme_limp = "cleaning frequency",
                vio_roedores_frec = "frequency of rodent sight")
              ) %>% 
  
  ## Table format
  add_overall() %>% 
  italicize_labels() %>% 
  
  ## Add significance
  add_p() %>% 
  bold_p()

## Save to word
tab1 %>% as_flex_table() %>% 
  font(fontname = "calibri", part = "all") %>% 
  fontsize(size = 12, part = "all") %>% 
  line_spacing(space = 1.5, part = "all") %>% 
  save_as_docx(path = "table1.docx")


# Table 2: univariate GLMMs for dogs --------------------------------------
tab2 <- imu_dog %>% 
  select(-adm_dis) %>% 
  tbl_uvregression(
    y = MAT_res, 
    method = glmmTMB::glmmTMB,
    # formula = "{y} ~ {x} + (1|adm_dis)",
    formula = "{y} ~ {x}",
    method.args = list(family = binomial),
    exponentiate = T,
    show_single_row = vac_any| 
      starts_with("street")|
      starts_with("con")|
      starts_with("hunt")|
      mud_grass_yard
    ) %>% 
  bold_p()


# LASSO selection ---------------------------------------------------------
### Control structure
ctrl <- trainControl(method = "repeatedcv", 
             number = 5, 
             repeats = 10)

### LASSO
lasso_fit <- train(form = MAT_res ~ age_yrs + sex + breed_cat + BCS_cat +
                     vac_any + animal_feed + street_access_uns + con_garbage_dumps +
                     con_mud_water + con_dogs + con_cats + con_anim_others + 
                     con_anim_house + con_anim_neighbors + con_anim_stray +
                     n_dogs_t + n_cats_t + n_pets_t + hunt_animals + hunt_rodents +
                     hous_type + hous_fr_clean + mud_grass_yard,  
                   data = imu_dog,
                   family = binomial,
                   method = "glmnet", 
                   trControl = ctrl)
                  
### Predictors
pred <- predictors(lasso_fit) %>% 
  str_extract(string = ., pattern = paste(colnames(imu_dog), collapse = "|")) %>% 
  unique()

### Regression table
## Update variables
imu_dog <- imu_dog %>% select(MAT_res, pred) 

### Geenerate table
reg_table <- glm(MAT_res ~ ., data = imu_dog, family = binomial) %>% 
  tbl_regression(exponentiate = T) %>% 
  bold_p()


# Logistic regression models ----------------------------------------------
### Multivariate model
fit = glmmTMB(MAT_res ~ street_access + con_garbage_dumps + rodent_sight_fr +
                (1|adm_dis), 
              family = binomial, data = imu_dog)

### Variable selection
drop1(fit)

fit1 = update(fit, ~.-con_garbage_dumps)
drop1(fit1)

fit2 = update(fit1, ~.-rodent_sight_fr)

# Compare models
compare_performance(fit, fit1, fit2, metrics = "AIC", rank = T)

# Check model residuals
testResiduals(fit2)

### Coefficients
tbl_regression(fit2, exponentiate = T)

r2(fit2)