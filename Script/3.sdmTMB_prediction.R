
# 1. Chargement des données -----------------------------------------------

# On charge la table de données préparées pour la modélisation CPUE
data_modelo = readRDS("data/outputs/df/data_modelo.rds")


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

dat_sf <- st_transform(dat_sf, 32742)

coords <- st_coordinates(dat_sf)

dat_clean$X <- coords[,1] / 1000
dat_clean$Y <- coords[,2] / 1000

Study_domain_Hole <- readRDS("Data/Ker_polygon.rds")
Study_domain_Hole <- st_transform(Study_domain_Hole, 32742)

# Update CRS to indicate units in kilometers (custom CRS)
Study_domain_Hole <- Study_domain_Hole %>%
  mutate(geometry = st_geometry(.) / 1000)
st_crs(Study_domain_Hole) <- "+proj=utm +zone=42 +datum=WGS84 +units=km +no_defs"


mesh <- make_mesh(
  dat_clean,
  xy_cols = c("X", "Y"),
  n_knots = 200
)

plot(mesh)

# Add Barrier to account for donut effect of the island
bspde <- sdmTMBextra::add_barrier_mesh(
  mesh,
  Study_domain_Hole,
  range_fraction = 0.1,
  proj_scaling = 1,
  plot = TRUE
)


# Plotting the mesh and barrier
mesh_df_water <- bspde$mesh_sf[bspde$normal_triangles, ]
mesh_df_land  <- bspde$mesh_sf[bspde$barrier_triangles, ]

ggplot() +
  geom_sf(data = Study_domain_Hole, fill = NA, colour = "red", size = 1) +
  geom_sf(data = mesh_df_water, colour = "blue") +
  geom_sf(data = mesh_df_land, colour = "green") +
  ggtitle("Mesh avec effet barrière") +
  theme_minimal()

plot(bspde)


# Estimation du modèle GAMMA simple
m_gamma <- sdmTMB(
  formula = CPUE ~ prof_bathy,
  data = dat_clean,
  mesh = bspde,
  family = tweedie(link = "log"),
  spatial = TRUE   
)



# Pour modèle sans temps
# 5. Visualisation relation -----------------------------------------------
# Vérification du modèle + graphes de relations CPUE ~ variables explicatives

sanity(m_gamma)

AIC(m_gamma)

dat_clean$resids <- residuals(m_gamma)
hist(dat_clean$resids)

qqnorm(dat_clean$resids)
abline(a = 0, b = 1)
ggplot(dat_clean, aes(X, Y, col = resids)) + scale_colour_gradient2() +
  geom_point() + facet_wrap(~Year) + coord_fixed()

# Graphique : relation CPUE / profondeur
pred_depth <- ggeffects::ggpredict(m_gamma, terms = "prof_bathy [all]")

depth_graph = ggplot(pred_depth, aes(x = x, y = predicted)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_point(
    data = dat_clean,
    aes(x = prof_bathy, y = CPUE),
    inherit.aes = FALSE,
    color = "red",
    alpha = 0.5
  ) +
  labs(x = "Profondeur (en mètres)", y = "CPUE de la palangre (en g/hameçon)") +
  theme_bw()


if (SavePlots) {
  ggsave(
    filename = file.path("data/", output_plot, "relation_CPUE_Profondeur.png"),
    plot     = depth_graph,
    width    = 8,
    height   = 6,
    dpi      = 300
  )
}



# 6. Résidus --------------------------------------------------------------
# Extraction et visualisation des résidus du modèle
dat_clean$resids <- residuals(m_gamma, type = "mle-mvn")

# Histogramme des résidus

hist_residual <- ggplot(dat_clean, aes(x = resids)) +
  geom_histogram(
    bins = 30,
    fill = "steelblue",
    color = "black",
    alpha = 0.7
  ) +
  labs(title = "Histogramme des résidus", x = "Résidus", y = "Fréquence") +
  theme_bw()

if (SavePlots) {
  ggsave(
    filename = file.path("data/", output_plot, "hist_residuals.png"),
    plot     = hist_residual,
    width    = 8,
    height   = 6,
    dpi      = 300
  )
}

# QQplot des résidus

qqplot_resdat = dat_clean %>% ggplot(aes(sample = resids)) +
  geom_abline(slope = 1,
              intercept = 0,
              col = 'red') +
  stat_qq() +
  stat_qq_line() + theme_bw()

if (SavePlots) {
  ggsave(
    filename = file.path("data/", output_plot, "qqplot_residual.png"),
    plot     = qqplot_resdat,
    width    = 8,
    height   = 6,
    dpi      = 300
  )
}

# Simulation des résidus de type DHARMa

s <- simulate(m_gamma, nsim = number_sim, type = "mle-mvn")
dharma_qqplot = dharma_residuals(s, m_gamma, test_uniformity = TRUE)


# 7. Grilles d’extrapolation ----------------------------------------------
# Import et préparation de la grille d’extrapolation spatiale
Extrapol_grid <- readRDS(here("data/outputs/df/grid_extrapol.rds"))

# Mise à l’échelle de la profondeur et de la pente sur la grille
Mean_SD <- dat_clean %>% ungroup() %>% summarise(Mean = mean(prof_bathy), SD = sd(prof_bathy))
Extrapol_grid = Extrapol_grid %>% rename(depth = profondeur.pred)
Extrapol_grid <- Extrapol_grid %>% mutate(depth_scaled = (depth - Mean_SD$Mean) / Mean_SD$SD)

Mean_SD <- dat_clean %>% ungroup() %>% summarise(Mean = mean(slope), SD = sd(slope))
Extrapol_grid <- Extrapol_grid %>% mutate(slope_scaled = (slope - Mean_SD$Mean) / Mean_SD$SD)

# Nettoyage et sélection des colonnes nécessaires pour la prédiction
grid <- Extrapol_grid %>%
  mutate(
    nb_ham = nb_ham,
    tps_h = mean(dat$tps_h),
    x = x / 1000,
    y = y / 1000
  ) %>%
  filter(!is.na(slope_scaled)) %>%
  filter(!is.na(depth_scaled)) %>%
  filter(depth > min_depth) %>%
  rename(X = x, Y = y) %>%
  dplyr::select(depth_scaled, slope_scaled, slope, depth, nb_ham, tps_h, X, Y) %>%
  rename(prof_bathy = depth)


# 8. Prédiction -----------------------------------------------------------
# Prédictions du modèle sur la grille avec séparation des effets
predictions <- predict(m_gamma, newdata = grid, return_tmb_object = TRUE)
grid$est <- predictions$data$est
grid$est_non_rf <- predictions$data$est_non_rf
grid$est_rf <- predictions$data$est_rf
grid$omega_s  <- predictions$data$omega_s
grid$epsilon_st <- predictions$data$epsilon_st

# Transformation des coordonnées pour passer en latitude / longitude
grid$X <- grid$X * 1000
grid$Y <- grid$Y * 1000
grid_sf <- st_as_sf(grid, coords = c("X", "Y"), crs = 32742)
grid_sf_transformed <- st_transform(grid_sf, crs = 4326)
grid$Lon <- st_coordinates(grid_sf_transformed)[, 1]
grid$Lat <- st_coordinates(grid_sf_transformed)[, 2]



# 9. Génération des cartes de prédiction ----------------------------------
# Génération des cartes de prédictions :
# - prédiction globale
# - prédiction effets fixes uniquement
# - cartes des effets aléatoires (omega)
r <- raster(extent(60, 76, -54, -44), resolution = c(0.0150, 0.0150))
raster_data <- rasterize(grid[, c("Lon", "Lat")], r, field = grid$est)
raster_df <- as.data.frame(raster_data, xy = TRUE)

min = min(exp(raster_df$layer), na.rm = TRUE)
max = max(exp(raster_df$layer), na.rm = TRUE)


# Carte de prédiction globale
#load("data/KER_SHELF_map.RData")
police_size <- 14

global_predict = KER_SHELF_map +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = exp(layer))) +
  scale_fill_gradientn(
    name = "CPUE (g par ham)",
    trans = "log",
    colours = c("#440154", "#31688E", "#35B779", "#FDE725"),
    # palette fixe (ex : viridis)
    limits = c(min, max),
    na.value = NA
  ) +
  theme_bw() +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("Prédiction globale de la CPUE") +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = police_size),
    axis.title.y = element_text(size = police_size),
    axis.text.x = element_text(size = police_size),
    axis.text.y = element_text(size = police_size),
    legend.title = element_text(size = police_size),
    legend.text = element_text(size = police_size),
    strip.text = element_text(size = police_size)
  )

if (SavePlots) {
  ggsave(
    filename = file.path("data/", output_plot, "global_predict.png"),
    plot     = global_predict,
    width    = 8,
    height   = 6,
    dpi      = 300
  )
}

# Carte de prédiction effets fixes uniquement
r <- raster(extent(60, 76, -54, -44), resolution = c(0.0150, 0.0150))
raster_data <- rasterize(grid[, c("Lon", "Lat")], r, field = grid$est_non_rf)
raster_df <- as.data.frame(raster_data, xy = TRUE)

#load("data/KER_SHELF_map.RData")
police_size <- 16

fixed_effect_predict = KER_SHELF_map +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = exp(layer))) +
  scale_fill_gradientn(
    name = "CPUE (g par ham)",
    trans = "log",
    colours = c("#440154", "#31688E", "#35B779", "#FDE725"),
    # palette fixe (ex : viridis)
    limits = c(min, max),
    na.value = NA
  ) +
  theme_bw() +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("Prédiction uniquement à partir des effets fixes") +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = police_size),
    axis.title.y = element_text(size = police_size),
    axis.text.x = element_text(size = police_size),
    axis.text.y = element_text(size = police_size),
    legend.title = element_text(size = police_size),
    legend.text = element_text(size = police_size),
    strip.text = element_text(size = police_size)
  )

if (SavePlots) {
  ggsave(
    filename = file.path("data/", output_plot, "fixed_effect_predict.png"),
    plot     = fixed_effect_predict,
    width    = 8,
    height   = 6,
    dpi      = 300
  )
}


# Carte des effets spatiaux aléatoires omega
r <- raster(extent(60, 76, -54, -44), resolution = c(0.0150, 0.0150))
raster_data <- rasterize(grid[, c("Lon", "Lat")], r, field = grid$omega_s)
raster_df <- as.data.frame(raster_data, xy = TRUE)

#load("data/KER_SHELF_map.RData")
police_size <- 16

omega_effect = KER_SHELF_map +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = layer)) +
  scale_fill_viridis_c(name = "Omega_s", na.value = NA)  +
  theme_bw() + xlab("Longitude") + ylab("Latitude") +
  ggtitle("Effets spatiaux aléatoires uniquement ~ omega")

if (SavePlots) {
  ggsave(
    filename = file.path("data/", output_plot, "omega_effect.png"),
    plot     = omega_effect,
    width    = 8,
    height   = 6,
    dpi      = 300
  )
}




