
# 1. Chargement des données -----------------------------------------------

# On charge la table de données préparées pour la modélisation CPUE
data_modelo = readRDS("data/outputs/df/data_modelo-crozet.rds")


# 2. Préparation des données ----------------------------------------------
# Selon la valeur de "repartition", on prépare les données différemment :
# - Cas repartition = TRUE : normalisation par ID, calcul du nb_ham et du poids par section
# - Cas repartition = FALSE : calcul direct de la CPUE, renommage des variables

if (repartition) {
  spat_peche_format <- data_modelo %>%
    filter(!is.na(prof_bathy), !is.na(slope)) %>%
    mutate(prof_bathy = abs(prof_bathy)) %>%
    group_by(StationID) %>%
    mutate(
      nb_ham = HookHaul / max(ID, na.rm = TRUE),
      pds_section = ScanWeightKg / max(ID, na.rm = TRUE)
    ) %>%
    group_by(StationID) %>%
    mutate(
      CPUE  = (pds_section * 1000) / nb_ham,
      tps_h = round(tps_trempe_hm_num, digits = 0)
    )
  
} else {
  spat_peche_format <- data_modelo %>%
    filter(!is.na(prof_bathy), !is.na(slope)) %>%
    mutate(prof_bathy = abs(prof_bathy)) %>%
    group_by(StationID) %>%
    mutate(
      CPUE  = (ScanWeightKg * 1000) / HookHaul,
      tps_h = round(tps_trempe_hm_num, digits = 0)
    ) %>%
    rename(nb_ham = HookHaul)
}

# On associe les captures aux points médians de chaque ligne de pêche
df_filtre <- spat_peche_format %>%
  dplyr::filter(if_all(all_of(filtres_actifs), ~ .x == 1))


# 3. Transformation -------------------------------------------------------
# Mise à l’échelle des covariables continues (centrer-réduire)
# Ajout des colonnes UTM pour la modélisation spatiale
dat <- df_filtre %>%
  ungroup() %>%
  mutate(
    depth_scaled = (prof_bathy - mean(prof_bathy)) / sd(prof_bathy),
    MDepth = mean(prof_bathy),
    SDDepth = sd(prof_bathy),
    slope_scaled = (slope - mean(slope)) / sd(slope)
  )
dat <- sdmTMB::add_utm_columns(dat, ll_names = c("x", "y"), units = "km")  # Conversion coordonnées



# 4. Mesh -----------------------------------------------------------------

# Création du maillage spatial pour le modèle (choix du nombre de nœuds)
noeuds <- 120


# Supprime les lignes avec CPUE manquante
dat_clean <- dat %>% filter(!is.na(CPUE), !is.na(prof_bathy), !is.na(Year))

dat_sf <- st_as_sf(dat_clean, coords = c("x","y"), crs = 4326)

dat_sf <- st_transform(dat_sf, 32743)

coords <- st_coordinates(dat_sf)

dat_clean <- dat_clean %>%
  mutate(
    X = coords[,1] / 1000,
    Y = coords[,2] / 1000
  )

dat_pres <- dat_clean %>%
  mutate(pres = as.integer(CPUE > 0))

dat_pos <- dat_pres %>% filter(CPUE > 0, prof_bathy >= 500)

q95 <- quantile(dat_pos$CPUE, 0.95, na.rm = TRUE)

dat_pos <- dat_pos %>%
  mutate(CPUE = pmin(CPUE, q95))

dat_pos <- dat_pos %>%
  mutate(
    NB_MEAN_CACHALOT = case_when(
      !is.na(nb_min_cachalot) & !is.na(nb_max_cachalot) ~ (nb_min_cachalot + nb_max_cachalot) / 2,
      is.na(nb_min_cachalot) & !is.na(nb_max_cachalot) ~ nb_max_cachalot,
      !is.na(nb_min_cachalot) & is.na(nb_max_cachalot) ~ nb_min_cachalot,
      TRUE ~ NA_real_
    )
  )

dat_pos <- dat_pos %>%
  mutate(
    NB_MEAN_ORQUE = case_when(
      !is.na(nb_min_orque) & !is.na(nb_max_orque) ~ (nb_min_orque + nb_max_orque) / 2,
      is.na(nb_min_orque) & !is.na(nb_max_orque) ~ nb_max_orque,
      !is.na(nb_min_orque) & is.na(nb_max_orque) ~ nb_min_orque,
      TRUE ~ NA_real_
    )
  )

dat_pos <- dat_pos %>%
  mutate(
    
    # indicateur d'observation
    cachalot_obs = case_when(
      pres_cachalot %in% c("présence", "absence") ~ 1,
      pres_cachalot %in% c("non obsé", NA) ~ 0
    ),
    
    orque_obs = case_when(
      pres_orque %in% c("présence", "absence") ~ 1,
      pres_orque %in% c("non obsé", NA) ~ 0
    ),
    
    # mise à zéro des valeurs si pas observé
    NB_MEAN_CACHALOT = ifelse(cachalot_obs == 0, 0, NB_MEAN_CACHALOT),
    NB_MEAN_ORQUE    = ifelse(orque_obs == 0, 0, NB_MEAN_ORQUE)
  )


dat_pos <- dat_pos %>%
  filter(
    !TAAFseason %in%  2026
  )



dat_pos <- dat_pos %>%
  filter(
    !is.na(NB_MEAN_CACHALOT),
    !is.na(NB_MEAN_ORQUE)
  )

dat_pos$Vessel <- droplevels(dat_pos$Vessel)


dat_pos$HookLength <- as.numeric(sub("_.*", "", dat_pos$HookSize))


dat_pos <- dat_pos %>%
  filter(HookLength != 7)


mean_value <- mean(dat_pos$HookLength, na.rm = TRUE)


dat_pos$HookLength[is.na(dat_pos$HookLength)] <- mean_value






Study_domain_Hole <- readRDS("Data/Ker_polygon.rds")%>%
  mutate(geometry = st_geometry(.) / 1000)

# Update CRS to indicate units in kilometers (custom CRS)

st_crs(Study_domain_Hole) <- "+proj=utm +zone=42 +datum=WGS84 +units=km +no_defs"



dat_pos_std <- dat_pos %>%
  mutate(
    TAAFseason_sc = as.numeric(scale(TAAFseason)),
    prof_bathy_sc = as.numeric(scale(prof_bathy)),
    slope_sc = as.numeric(scale(slope)),
    HookLength_sc = as.numeric(scale(HookLength)),
    NB_MEAN_CACHALOT_sc = as.numeric(scale(NB_MEAN_CACHALOT)),
    NB_MEAN_ORQUE_sc = as.numeric(scale(NB_MEAN_ORQUE))
  )







# 
# # ------------------------------------------------------------------------------#
# ####  COMPUTE AND PLOT MORAN's INDEX = SPATIAL AUTOCORRELATION  #### 
# # ------------------------------------------------------------------------------#
# # source(here('Function', 'function_model', 'compute_plot_Moran.R'))
# 
# # ------------------------------------------------------------------------------#
# ####  PLOT RESIDUALS OF MODEL FITTED WITHOUT RANDOM FIELD #### 
# # ------------------------------------------------------------------------------#
# source(here('Function', 'function_model', 'compute_plot_residuals_withoutRF.R'))

# ------------------------------------------------------------------------------#
####  CREATE THE MESH CONFIGURATIONS #### 
# ------------------------------------------------------------------------------#
source(here('Function', 'function_model', 'MSA_add_boundary_to_mesh.R'))
source(here('Function', 'function_model', 'create_mesh_all_configs.R'))


params <- list(data_m = dat_pos,
               cutoff_values = c(10, 11, 15),
               max_edge_in  = c( 5, 20, 20),
               max_edge_out = rep(100, 6),
               utm_crs_used = st_crs(32743),
               land_domain = Study_domain_Hole)


meshes <- do.call(create_meshes_configurations, params)


# ------------------------------------------------------------------------------#
####  BUILD & FIT CANDITATE MODELS  #### 
# ------------------------------------------------------------------------------#
# Vector of response variables to be tested sequentially
responses <- c("ScanWeightKg")

# Fixed-effects
fixed_effect <- "TAAFseason + Month + Vessel + prof_bathy_sc + slope_sc + HookLength_sc + cachalot_obs + orque_obs + NB_MEAN_CACHALOT_sc * NB_MEAN_ORQUE_sc"

# candidate_terms <- c(
#   "s(TAAFseason_sc)",
#   "s(Month, bs = 'cc')",      # on ne standardise pas un mois cyclique
#   "Vessel",
#   "s(prof_bathy_sc)",
#   "s(slope_sc)",
#   "s(HookLength_sc)",
#   "cachalot_obs",
#   "orque_obs",
#   "s(NB_MEAN_CACHALOT_sc, NB_MEAN_ORQUE_sc)"
# )
# 
# selection <- forward_selection_fixed_effects(
#   data = dat_pos_std,
#   response = "ScanWeightKg",
#   mesh = meshes$cutoff_50$mesh,
#   family = sdmTMB::tweedie(link = "log"),
#   candidate_terms = candidate_terms
# )
# 
# selection$selected_terms
# selection$final_formula
# selection$history
# selection$best_model

source(here('Function', 'function_model', 'fit_candidate_models_all_mesh_configs.R'))



fitted_models <- fit_candidate_models(
  data_CGFS = dat_pos_std,
  meshes_configurations = meshes
)

#ajout poisson

# ------------------------------------------------------------------------------#
#### SANITY FILTER FITTED MODELS ####
# ------------------------------------------------------------------------------#

source(here('Function', 'function_model', 'sanity_filter_models_all_mesh_configs.R'))

sanity_models <- list()
# Apply sanity checks to West English Channel models
sanity_models <- sanity_filter_models(fitted_models)



# ------------------------------------------------------------------------------#
#### EXTRACT CONVERGED MODELS and THEIR PARAMETER ESTIMATES ####
# ------------------------------------------------------------------------------#

source(here('Function', 'function_model', 'get_converged_models_all_mesh_configs.R'))
converged_models <- get_converged_models(sanity_models)

source(here('Function', 'function_model', 'extract_parameter_estimates_all_mesh_configs.R'))
parameters <- extract_params_converged_models(converged_models)


# ------------------------------------------------------------------------------#
####  REPEATED RANDOM CROSS-VALIDATION #### 
# ------------------------------------------------------------------------------#
source(here('Function', 'function_model', 'cross_validation_repeated_random.R'))
repeated_randomCV <- list()


repeated_randomCV <- repeated_random_CV(
  converged_models = converged_models, meshes = meshes, data_CGFS = dat_pos_std, n_repeats = 10, k_folds = 10, seed = 123)



# ------------------------------------------------------------------------------#
####  conditional AIC (cAIC)  #### 
# ------------------------------------------------------------------------------#
source(here('Function', 'function_model', 'AIC_weights_all_mesh_configs.R'))



AIC_models <- compute_cAIC_weights(converged_models)
