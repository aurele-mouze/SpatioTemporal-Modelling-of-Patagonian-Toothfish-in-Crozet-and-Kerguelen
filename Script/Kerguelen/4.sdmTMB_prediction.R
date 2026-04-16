
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


Study_domain_Hole <- readRDS("Data/Ker_polygon.rds")%>%
  mutate(geometry = st_geometry(.) / 1000)

# Update CRS to indicate units in kilometers (custom CRS)

st_crs(Study_domain_Hole) <- "+proj=utm +zone=42 +datum=WGS84 +units=km +no_defs"


mesh_pos <- make_mesh(dat_pos, xy_cols = c("X", "Y"), n_knots = 200)


# Add Barrier to account for donut effect of the island
bspde_pos <- sdmTMBextra::add_barrier_mesh(mesh_pos, Study_domain_Hole,
                                           range_fraction = 0.1, proj_scaling = 1, plot = TRUE)

cat("Mesh : N nœuds =", mesh_pos$mesh$n, "\n\n")

# Visualisation mesh

loc  <- mesh_pos$mesh$loc[, 1:2]
loc_df <- data.frame(loc)
names(loc_df) <- c("X", "Y")
idx  <- mesh_pos$mesh$graph$tv[, 1:3]
segs <- lapply(seq_len(nrow(idx)), function(i) {
  pts <- loc[idx[i, ], ]
  st_linestring(rbind(pts, pts[1, ]))
})
mesh_lines <- st_sf(
  geometry = st_sfc(segs),
  crs = st_crs(Study_domain_Hole)  # utiliser le CRS de ton domaine
)
mesh_nodes <- st_as_sf(
  loc_df,
  coords = c("X", "Y"),
  crs = st_crs(Study_domain_Hole)
)


# Créer le plot
mesh_plot <- ggplot() +
  geom_sf(data = mesh_lines, color = "grey55", linewidth = 0.2) +
  geom_sf(data = mesh_nodes, aes(color = "Nodes"), size = 1) +  #Si on veut des points sur les noeuds
  # geom_point(data = dat_pos, aes(x = X, y = Y),  #Si on veut afficher les stations
  #            color = "orange", size = 1.2) +
  geom_sf(data = Study_domain_Hole, fill = "grey80", color = "black") +
  coord_sf(
    xlim = range(dat_pos$X) + c(-35, 35),
    ylim = range(dat_pos$Y) + c(-35, 35)
  ) +
  labs(
    title = paste0("Mesh"),
    subtitle = paste0("N nœuds = ", mesh_pos$mesh$n),
    x = "Eastings (km)", y = "Northings (km)"
  ) +
  theme_bw() +
  theme(panel.grid = element_blank())


if (SavePlots) {
  ggsave(
    filename = file.path("data/", output_plot, "mesh_plot.png"),
    plot     = mesh_plot,
    width    = 8,
    height   = 6,
    dpi      = 300
  )
}
cat("→ mesh.png sauvegardé\n\n")




# # Estimation du modèle GAMMA simple, sans effets
# m_null <- sdmTMB(
#   CPUE ~ 1,
#   data = dat_pos,
#   mesh = bspde_pos,
#   family = Gamma(link = "log"),
#   spatial = "off",
#   spatiotemporal = "off"
# )
# 
# m_spatial <- sdmTMB(
#   CPUE ~ 1,
#   data = dat_pos,
#   mesh = bspde_pos,
#   family = Gamma(link = "log"),
#   spatial = "on",
#   spatiotemporal = "off"
# )
# 
# 
m_spatemp <- sdmTMB(
  CPUE ~ 1,
  data = dat_pos,
  mesh = bspde_pos,
  family = Gamma(link = "log"),
  spatial = "on",
  time= "Year",
  spatiotemporal = "iid"
)
# 
# 
# m_0 <- sdmTMB(
#   CPUE ~ s(Year),
#   data = dat_pos,
#   mesh = bspde_pos,
#   family = Gamma(link = "log"),
#   spatial = "on",
#   time= "Year",
#   spatiotemporal = "iid"
# )

m_1 <- sdmTMB(
  ScanWeightKg ~ s(Year),
  offset = log(dat_pos$nb_ham),
  data = dat_pos,
  mesh = bspde_pos,
  family = Gamma(link = "log"),
  spatial = "on",
  time= "Year",
  spatiotemporal = "iid"
)



ll_1 <- as.numeric(logLik(m_1))
ll_spat <- as.numeric(logLik(m_spatemp))

delta_x1 <- 2 * (ll_1 - ll_spat)


#5. Contrôle des modèles ---------------------------------
# 
# sanity(m_null)
# sanity(m_spatial)
# sanity(m_spatemp)
# 
# dat_pos$resids <- residuals(m_null)
# hist_residual <- ggplot(dat_pos, aes(x = resids)) +
#   geom_histogram(
#     bins = 30,
#     fill = "steelblue",
#     color = "black",
#     alpha = 0.7
#   ) +
#   labs(title = "Histogramme des résidus", x = "Résidus", y = "Fréquence") +
#   theme_bw()
# 
# s <- simulate(m_null, nsim = number_sim, type = "mle-mvn")
# dharma_qqplot = dharma_residuals(s, m_null, test_uniformity = TRUE)
# 
# QQplot <- ggplot(dharma_qqplot, aes(x = expected, y = observed)) +
#   geom_point(alpha = 0.4) +
#   geom_abline(slope = 1, intercept = 0) +
#   theme_minimal()
# 
# if (SavePlots) {
#   ggsave(
#     filename = file.path("data/", output_plot, "residual-histogramme-m-null.png"),
#     plot     = hist_residual,
#     width    = 8,
#     height   = 6,
#     dpi      = 300
#   )
# }
# 
# if (SavePlots) {
#   ggsave(
#     filename = file.path("data/", output_plot, "qqplot_residual-m-null.png"),
#     plot     = QQplot,
#     width    = 8,
#     height   = 6,
#     dpi      = 300
#   )
# }
# 
# 
# dat_pos$resids <- residuals(m_spatial)
# hist_residual <- ggplot(dat_pos, aes(x = resids)) +
#   geom_histogram(
#     bins = 30,
#     fill = "steelblue",
#     color = "black",
#     alpha = 0.7
#   ) +
#   labs(title = "Histogramme des résidus", x = "Résidus", y = "Fréquence") +
#   theme_bw()
# 
# s <- simulate(m_spatial, nsim = number_sim, type = "mle-mvn")
# dharma_qqplot = dharma_residuals(s, m_spatial, test_uniformity = TRUE)
# 
# QQplot <- ggplot(dharma_qqplot, aes(x = expected, y = observed)) +
#   geom_point(alpha = 0.4) +
#   geom_abline(slope = 1, intercept = 0) +
#   theme_minimal()
# 
# if (SavePlots) {
#   ggsave(
#     filename = file.path("data/", output_plot, "residual-histogramme-spatial.png"),
#     plot     = hist_residual,
#     width    = 8,
#     height   = 6,
#     dpi      = 300
#   )
# }
# 
# if (SavePlots) {
#   ggsave(
#     filename = file.path("data/", output_plot, "qqplot_residual-spatial.png"),
#     plot     = QQplot,
#     width    = 8,
#     height   = 6,
#     dpi      = 300
#   )
# }
# 
# 
# 
# 
# dat_pos$resids <- residuals(m_spatemp)
# hist_residual <- ggplot(dat_pos, aes(x = resids)) +
#   geom_histogram(
#     bins = 30,
#     fill = "steelblue",
#     color = "black",
#     alpha = 0.7
#   ) +
#   labs(title = "Histogramme des résidus", x = "Résidus", y = "Fréquence") +
#   theme_bw()
# 
# s <- simulate(m_spatemp, nsim = number_sim, type = "mle-mvn")
# dharma_qqplot = dharma_residuals(s, m_spatemp, test_uniformity = TRUE)
# 
# QQplot <- ggplot(dharma_qqplot, aes(x = expected, y = observed)) +
#   geom_point(alpha = 0.4) +
#   geom_abline(slope = 1, intercept = 0) +
#   theme_minimal()
# 
# if (SavePlots) {
#   ggsave(
#     filename = file.path("data/", output_plot, "residual-histogramme-m-spatemp.png"),
#     plot     = hist_residual,
#     width    = 8,
#     height   = 6,
#     dpi      = 300
#   )
# }
# 
# if (SavePlots) {
#   ggsave(
#     filename = file.path("data/", output_plot, "qqplot_residual-m-spatemp.png"),
#     plot     = QQplot,
#     width    = 8,
#     height   = 6,
#     dpi      = 300
#   )
# }




#6. Fabrication d'une grille d'extrapolation ------------------------------------------------------
# Import et préparation de la grille d’extrapolation spatiale
Extrapol_grid <- readRDS(here("data/outputs/df/grid_extrapol.rds"))

# Mise à l’échelle de la profondeur et de la pente sur la grille
Mean_SD <- dat_pos %>% ungroup() %>% summarise(Mean = mean(prof_bathy), SD = sd(prof_bathy))
Extrapol_grid = Extrapol_grid %>% rename(depth = profondeur.pred)
Extrapol_grid <- Extrapol_grid %>% mutate(depth_scaled = (depth - Mean_SD$Mean) / Mean_SD$SD)

Mean_SD <- dat_pos %>% ungroup() %>% summarise(Mean = mean(slope), SD = sd(slope))
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

#7. Prédictions -----------------------------------
# # Vérifications rapides
# stopifnot("Year" %in% names(dat_pos))
# stopifnot(is.numeric(dat_pos$Year))
# 
# #Préparer grid avec un identifiant de cellule
# grid_with_id <- grid %>%
#   mutate(cell_id = row_number())
# 
# #Réplication temporelle : une ligne par cellule × année
# years <- sort(unique(dat_pos$Year))
# grid_all_years <- grid_with_id %>%
#   tidyr::uncount(length(years)) %>%
#   mutate(Year = rep(years, times = nrow(grid_with_id))) %>%
#   relocate(cell_id)
# 
# #Vérifier colonnes requises
# required_cols <- c("prof_bathy", "Year", "X", "Y")
# missing_cols <- setdiff(required_cols, names(grid_all_years))
# if (length(missing_cols) > 0) stop("Colonnes manquantes dans newdata : ", paste(missing_cols, collapse = ", "))
# 
# # Prédiction (utilise chunking si nécessaire)
# pred_grid_raw <- predict_in_chunks(m_null, grid_all_years, chunk_size = 20000)
# 
# # Rattacher cell_id et transformer sur l'échelle réponse
# pred_grid <- pred_grid_raw %>%
#   mutate(cell_id = grid_all_years$cell_id,
#          CPUE_pred = est) %>%
#   left_join(grid_with_id %>% dplyr::select(cell_id, nb_ham, tps_h, slope_scaled), by = "cell_id")
# 
# histo <- ggplot(pred_grid, aes(x = "", y = exp(CPUE_pred))) +
#   geom_boxplot(fill = "lightgreen", outlier.color = "red") +
#   labs(
#     x = "",
#     y = "CPUE prédite",
#     title = "Boxplot des CPUE prédites par le modèle null"
#   ) +
#   theme_minimal()
# 
# if (SavePlots) {
#   ggsave(
#     filename = file.path("data/", output_plot, "mean_histo_null.png"),
#     plot     = histo,
#     width    = 8,
#     height   = 6,
#     dpi      = 300
#   )
# }
# 
# 
# 
# 
# 
# 
# 
# 
# stopifnot("Year" %in% names(dat_pos))
# stopifnot(is.numeric(dat_pos$Year))
# 
# #Préparer grid avec un identifiant de cellule
# grid_with_id <- grid %>%
#   mutate(cell_id = row_number())
# 
# #Réplication temporelle : une ligne par cellule × année
# years <- sort(unique(dat_pos$Year))
# grid_all_years <- grid_with_id %>%
#   tidyr::uncount(length(years)) %>%
#   mutate(Year = rep(years, times = nrow(grid_with_id))) %>%
#   relocate(cell_id)
# 
# #Vérifier colonnes requises
# required_cols <- c("prof_bathy", "Year", "X", "Y")
# missing_cols <- setdiff(required_cols, names(grid_all_years))
# if (length(missing_cols) > 0) stop("Colonnes manquantes dans newdata : ", paste(missing_cols, collapse = ", "))
# 
# # Prédiction (utilise chunking si nécessaire)
# pred_grid_raw <- predict_in_chunks(m_spatial, grid_all_years, chunk_size = 20000)
# 
# # Rattacher cell_id et transformer sur l'échelle réponse
# pred_grid <- pred_grid_raw %>%
#   mutate(cell_id = grid_all_years$cell_id,
#          CPUE_pred = est) %>%
#   left_join(grid_with_id %>% dplyr::select(cell_id, nb_ham, tps_h, slope_scaled), by = "cell_id")
# 
# # Moyenne temporelle par cellule spatiale
# pred_mean_cell <- pred_grid %>%
#   group_by(cell_id, X, Y, prof_bathy) %>%
#   summarise(CPUE_mean = mean(CPUE_pred, na.rm = TRUE), .groups = "drop")
# 
# #  Carte moyenne -> sortie tout effets confondus
# mean_map <- ggplot(pred_mean_cell, aes(x = X, y = Y, fill = exp(CPUE_mean))) +
#   geom_raster() +
#   scale_fill_viridis_c(option = "viridis", name = "CPUE moyenne") +
#   coord_equal() +
#   labs(title = "log CPUE moyenne (moyenne sur les années)", x = "X (km)", y = "Y (km)") +
#   theme_minimal()
# 
# if (SavePlots) {
#   ggsave(
#     filename = file.path("data/", output_plot, "mean_map_spatial.png"),
#     plot     = mean_map,
#     width    = 8,
#     height   = 6,
#     dpi      = 300
#   )
# }
# 
# 
# 
# 
# 
# 
# 
# # Vérifications rapides
# stopifnot("Year" %in% names(dat_pos))
# stopifnot(is.numeric(dat_pos$Year))
# 
# #Préparer grid avec un identifiant de cellule
# grid_with_id <- grid %>%
#   mutate(cell_id = row_number())
# 
# grid_with_id <- dplyr::as_tibble(grid_with_id)
# 
# #Réplication temporelle : une ligne par cellule × année
# years <- sort(unique(dat_pos$Year))
# grid_all_years <- grid_with_id %>%
#   tidyr::uncount(length(years)) %>%
#   mutate(Year = rep(years, times = nrow(grid_with_id))) %>%
#   relocate(cell_id)
# 
# #Vérifier colonnes requises
# required_cols <- c("prof_bathy", "Year", "X", "Y")
# missing_cols <- setdiff(required_cols, names(grid_all_years))
# if (length(missing_cols) > 0) stop("Colonnes manquantes dans newdata : ", paste(missing_cols, collapse = ", "))
# 
# # Prédiction (utilise chunking si nécessaire)
# pred_grid_raw <- predict_in_chunks(m_spatemp, grid_all_years, chunk_size = 20000)
# 
# 
# 
# # Rattacher cell_id et transformer sur l'échelle réponse
# pred_grid <- pred_grid_raw %>%
#   mutate(cell_id = grid_all_years$cell_id,
#          CPUE_pred = est) %>%
#   left_join(grid_with_id %>% dplyr::select(cell_id, nb_ham, tps_h, slope_scaled), by = "cell_id")
# 
# # Moyenne temporelle par cellule spatiale
# pred_mean_cell <- pred_grid %>%
#   group_by(cell_id, X, Y, prof_bathy) %>%
#   summarise(CPUE_mean = mean(CPUE_pred, na.rm = TRUE), .groups = "drop")
# 
# #  Carte moyenne -> sortie tout effets confondus
# mean_map <- ggplot(pred_mean_cell, aes(x = X, y = Y, fill = exp(CPUE_mean))) +
#   geom_raster() +
#   scale_fill_viridis_c(option = "viridis", name = "CPUE moyenne") +
#   coord_equal() +
#   labs(title = "log CPUE moyenne (moyenne sur les années)", x = "X (km)", y = "Y (km)") +
#   theme_minimal()
# 
# if (SavePlots) {
#   ggsave(
#     filename = file.path("data/", output_plot, "mean_map_spatiotemp.png"),
#     plot     = mean_map,
#     width    = 8,
#     height   = 6,
#     dpi      = 300
#   )
# }
# 
# 
# annees <- sort(unique(dat_pos$Year))
# 
# liste_cartes <- list()
# 
# 
# for (an in annees) {
#   # filtrer l'année
#   pred_year <- pred_grid %>%
#     filter(Year == an)
#   
#   # carte année courante
#   carte <- ggplot(pred_year, aes(x = X, y = Y, fill = exp(CPUE_pred))) +
#     geom_tile() +
#     scale_fill_viridis_c(option = "viridis", name = "CPUE prédite") +
#     coord_equal() +
#     labs(
#       title = paste("CPUE prédite - année", an),
#       x = "X (km)",
#       y = "Y (km)"
#     ) +
#     theme_minimal()
#   
#   if (SavePlots) {
#     ggsave(
#       filename = file.path("data/", output_plot, paste0("map_CPUE_pred_", an, ".png")),
#       plot = carte,
#       width = 8,
#       height = 6,
#       dpi = 300
#     )}
# }
# 
# 
# # Effet spatial
# grid_spatial <- grid %>%
#   mutate(Year = 2015)
# pred_spatial <- predict(m_spatemp, newdata = grid_spatial)
# 
# grid_spatial$effet_spatial <- pred_spatial$omega_s
# 
# lim <- max(abs(grid_spatial$effet_spatial), na.rm = TRUE)
# 
# Spat_effect <- ggplot() +
#   
#   # fond gris sur toute l'emprise
#   annotate(
#     "rect",
#     xmin = min(grid_spatial$X),
#     xmax = max(grid_spatial$X),
#     ymin = min(grid_spatial$Y),
#     ymax = max(grid_spatial$Y),
#     fill = "grey80",
#     alpha = 0.4
#   ) +
#   
#   # raster prédit
#   geom_raster(
#     data = grid_spatial,
#     aes(X, Y, fill = effet_spatial)
#   ) +
#   
#   scale_fill_gradient2(
#     low = "blue",
#     mid = "white",
#     high = "red",
#     midpoint = 0,
#     limits = c(-lim, lim),
#     name = "Effet spatial"
#   ) +
#   
#   coord_equal() +
#   theme_minimal()
# 
# if (SavePlots) {
#   ggsave(
#     filename = file.path("data/", output_plot, "Effet_spatial.png"),
#     plot     = Spat_effect,
#     width    = 8,
#     height   = 6,
#     dpi      = 300
#   )
# }
# 
# 
# # Effet Spatio-temporel
# 
# pred_spatiotemp <- predict(m_spatemp, newdata = grid_all_years)
# 
# grid_all_years$effet_spatiotemp <- pred_spatiotemp$epsilon_st
# 
# 
# annees <- sort(unique(grid_all_years$Year))
# 
# liste_cartes_spatiotemp <- list()
# 
# lim <- max(abs(grid_all_years$effet_spatiotemp), na.rm = TRUE)
# 
# for (an in annees) {
#   
#   pred_year <- grid_all_years %>%
#     filter(Year == an)
#   
#   carte <- ggplot() +
#     
#     annotate(
#       "rect",
#       xmin = min(grid_all_years$X),
#       xmax = max(grid_all_years$X),
#       ymin = min(grid_all_years$Y),
#       ymax = max(grid_all_years$Y),
#       fill = "grey90",
#       alpha = 0.4
#     ) +
#     
#     geom_raster(
#       data = pred_year,
#       aes(X, Y, fill = effet_spatiotemp)
#     ) +
#     
#     scale_fill_gradient2(
#       low = "blue",
#       mid = "white",
#       high = "red",
#       midpoint = 0,
#       limits = c(-lim, lim),
#       name = "Effet spatio-temporel"
#     ) +
#     
#     coord_equal() +
#     
#     labs(
#       title = paste("Effet spatio-temporel - année", an),
#       x = "X (km)",
#       y = "Y (km)"
#     ) +
#     
#     theme_minimal()
#   
#   if (SavePlots) {
#     ggsave(
#       filename = file.path("data/", output_plot, paste0("map_Effet_Spatiotemp_", an, ".png")),
#       plot = carte,
#       width = 8,
#       height = 6,
#       dpi = 300
#     )}
# }
# 
# 
# pred_fix <- predict(
#   m_spatemp,
#   newdata = grid_all_years,
#   re_form = NA
# )
# 
# grid_all_years$effet_fixe <- pred_fix$est
# unique(grid_all_years$effet_fixe)
# 
# 
# 
# # Vérifications rapides
# stopifnot("Year" %in% names(dat_pos))
# stopifnot(is.numeric(dat_pos$Year))
# 
# #Préparer grid avec un identifiant de cellule
# grid_with_id <- grid %>%
#   mutate(cell_id = row_number())
# 
# grid_with_id <- dplyr::as_tibble(grid_with_id)
# 
# #Réplication temporelle : une ligne par cellule × année
# years <- sort(unique(dat_pos$Year))
# grid_all_years <- grid_with_id %>%
#   tidyr::uncount(length(years)) %>%
#   mutate(Year = rep(years, times = nrow(grid_with_id))) %>%
#   relocate(cell_id)
# 
# #Vérifier colonnes requises
# required_cols <- c("prof_bathy", "Year", "X", "Y")
# missing_cols <- setdiff(required_cols, names(grid_all_years))
# if (length(missing_cols) > 0) stop("Colonnes manquantes dans newdata : ", paste(missing_cols, collapse = ", "))
# 
# # Prédiction (utilise chunking si nécessaire)
# pred_grid_raw <- predict_in_chunks(m_0, grid_all_years, chunk_size = 20000)
# 
# 
# 
# # Rattacher cell_id et transformer sur l'échelle réponse
# pred_grid <- pred_grid_raw %>%
#   mutate(cell_id = grid_all_years$cell_id,
#          CPUE_pred = est) %>%
#   left_join(grid_with_id %>% dplyr::select(cell_id, nb_ham, tps_h, slope_scaled), by = "cell_id")
# 
# # Moyenne temporelle par cellule spatiale
# pred_mean_cell <- pred_grid %>%
#   group_by(cell_id, X, Y, prof_bathy) %>%
#   summarise(CPUE_mean = mean(CPUE_pred, na.rm = TRUE), .groups = "drop")
# 
# #  Carte moyenne -> sortie tout effets confondus
# mean_map <- ggplot(pred_mean_cell, aes(x = X, y = Y, fill = exp(CPUE_mean))) +
#   geom_raster() +
#   scale_fill_viridis_c(option = "viridis", name = "CPUE moyenne") +
#   coord_equal() +
#   labs(title = "CPUE moyenne (moyenne sur les années)", x = "X (km)", y = "Y (km)") +
#   theme_minimal()
# 
# if (SavePlots) {
#   ggsave(
#     filename = file.path("data/", output_plot, "mean_map_m_0.png"),
#     plot     = mean_map,
#     width    = 8,
#     height   = 6,
#     dpi      = 300
#   )
# }
# 
# 
# annees <- sort(unique(dat_pos$Year))
# 
# liste_cartes <- list()
# 
# 
# for (an in annees) {
#   # filtrer l'année
#   pred_year <- pred_grid %>%
#     filter(Year == an)
#   
#   # carte année courante
#   carte <- ggplot(pred_year, aes(x = X, y = Y, fill = exp(CPUE_pred))) +
#     geom_tile() +
#     scale_fill_viridis_c(option = "viridis", name = "CPUE prédite") +
#     coord_equal() +
#     labs(
#       title = paste("CPUE prédite - année", an),
#       x = "X (km)",
#       y = "Y (km)"
#     ) +
#     theme_minimal()
#   
#   if (SavePlots) {
#     ggsave(
#       filename = file.path("data/", output_plot, paste0("map_CPUE_pred_m_0_", an, ".png")),
#       plot = carte,
#       width = 8,
#       height = 6,
#       dpi = 300
#     )}
# }
# 
# 
# # Effet spatial
# grid_spatial <- grid %>%
#   mutate(Year = 2015)
# pred_spatial <- predict(m_0, newdata = grid_spatial)
# 
# grid_spatial$effet_spatial <- pred_spatial$omega_s
# 
# lim <- max(abs(grid_spatial$effet_spatial), na.rm = TRUE)
# 
# Spat_effect <- ggplot() +
#   
#   # fond gris sur toute l'emprise
#   annotate(
#     "rect",
#     xmin = min(grid_spatial$X),
#     xmax = max(grid_spatial$X),
#     ymin = min(grid_spatial$Y),
#     ymax = max(grid_spatial$Y),
#     fill = "grey80",
#     alpha = 0.4
#   ) +
#   
#   # raster prédit
#   geom_raster(
#     data = grid_spatial,
#     aes(X, Y, fill = effet_spatial)
#   ) +
#   
#   scale_fill_gradient2(
#     low = "blue",
#     mid = "white",
#     high = "red",
#     midpoint = 0,
#     limits = c(-lim, lim),
#     name = "Effet spatial"
#   ) +
#   
#   coord_equal() +
#   theme_minimal()
# 
# if (SavePlots) {
#   ggsave(
#     filename = file.path("data/", output_plot, "Effet_spatial_m_0.png"),
#     plot     = Spat_effect,
#     width    = 8,
#     height   = 6,
#     dpi      = 300
#   )
# }
# 
# 
# # Effet Spatio-temporel
# 
# pred_spatiotemp <- predict(m_0, newdata = grid_all_years)
# 
# grid_all_years$effet_spatiotemp <- pred_spatiotemp$epsilon_st
# 
# 
# annees <- sort(unique(grid_all_years$Year))
# 
# liste_cartes_spatiotemp <- list()
# 
# lim <- max(abs(grid_all_years$effet_spatiotemp), na.rm = TRUE)
# 
# for (an in annees) {
#   
#   pred_year <- grid_all_years %>%
#     filter(Year == an)
#   
#   carte <- ggplot() +
#     
#     annotate(
#       "rect",
#       xmin = min(grid_all_years$X),
#       xmax = max(grid_all_years$X),
#       ymin = min(grid_all_years$Y),
#       ymax = max(grid_all_years$Y),
#       fill = "grey90",
#       alpha = 0.4
#     ) +
#     
#     geom_raster(
#       data = pred_year,
#       aes(X, Y, fill = effet_spatiotemp)
#     ) +
#     
#     scale_fill_gradient2(
#       low = "blue",
#       mid = "white",
#       high = "red",
#       midpoint = 0,
#       limits = c(-lim, lim),
#       name = "Effet spatio-temporel"
#     ) +
#     
#     coord_equal() +
#     
#     labs(
#       title = paste("Effet spatio-temporel - année", an),
#       x = "X (km)",
#       y = "Y (km)"
#     ) +
#     
#     theme_minimal()
#   
#   if (SavePlots) {
#     ggsave(
#       filename = file.path("data/", output_plot, paste0("map_Effet_Spatiotemp_m_0_", an, ".png")),
#       plot = carte,
#       width = 8,
#       height = 6,
#       dpi = 300
#     )}
# }
# 
# 
# pred_fix <- predict(
#   m_spatemp,
#   newdata = grid_all_years,
#   re_form = NA
# )
# 
# grid_all_years$effet_fixe <- pred_fix$est
# unique(grid_all_years$effet_fixe)
# 
# 
# #Effet année
# pred_fix <- predict(
#   m_spatemp,
#   newdata = grid_all_years,
#   re_form = NA
# )
# 
# grid_all_years$effet_annee <- pred_fix$est
# 
# year_effect <- grid_all_years %>%
#   group_by(Year) %>%
#   summarise(effet = mean(effet_annee))
# 
# Effet_An <- ggplot(year_effect, aes(Year, effet)) +
#   geom_line() +
#   geom_point() +
#   theme_minimal()
# 
# 
# if (SavePlots) {
#   ggsave(
#     filename = file.path("data/", output_plot, "Effet_annuel_m_0.png"),
#     plot     = Effet_An,
#     width    = 8,
#     height   = 6,
#     dpi      = 300
#   )
# }





#d. m_1
# Vérifications rapides
stopifnot("Year" %in% names(dat_pos))
stopifnot(is.numeric(dat_pos$Year))

#Préparer grid avec un identifiant de cellule
grid_with_id <- grid %>%
  mutate(cell_id = row_number())

grid_with_id <- dplyr::as_tibble(grid_with_id)

#Réplication temporelle : une ligne par cellule × année
years <- sort(unique(dat_pos$Year))
grid_all_years <- grid_with_id %>%
  tidyr::uncount(length(years)) %>%
  mutate(Year = rep(years, times = nrow(grid_with_id))) %>%
  mutate(nb_ham = 1) %>%
  relocate(cell_id)

#Vérifier colonnes requises
required_cols <- c("prof_bathy", "Year", "X", "Y")
missing_cols <- setdiff(required_cols, names(grid_all_years))
if (length(missing_cols) > 0) stop("Colonnes manquantes dans newdata : ", paste(missing_cols, collapse = ", "))

# Prédiction (utilise chunking si nécessaire)
pred_grid_raw <- predict_in_chunks(m_1, grid_all_years, chunk_size = 30000)



# Rattacher cell_id et transformer sur l'échelle réponse
pred_grid <- pred_grid_raw %>%
  mutate(cell_id = grid_all_years$cell_id,
         CPUE_pred = est) %>%
  left_join(grid_with_id %>% dplyr::select(cell_id, nb_ham, tps_h, slope_scaled), by = "cell_id")

# Moyenne temporelle par cellule spatiale
pred_mean_cell <- pred_grid %>%
  group_by(cell_id, X, Y, prof_bathy) %>%
  summarise(CPUE_mean = mean(CPUE_pred, na.rm = TRUE), .groups = "drop")

#  Carte moyenne -> sortie tout effets confondus
mean_map <- ggplot(pred_mean_cell, aes(x = X, y = Y, fill = exp(CPUE_mean))) +
  geom_raster() +
  scale_fill_viridis_c(option = "viridis", name = "Biomasse moyenne") +
  coord_equal() +
  labs(title = "Biomasse moyenne (moyenne sur les années)", x = "X (km)", y = "Y (km)") +
  theme_minimal()

if (SavePlots) {
  ggsave(
    filename = file.path("data/", output_plot, "mean_map_m_1.png"),
    plot     = mean_map,
    width    = 8,
    height   = 6,
    dpi      = 300
  )
}


annees <- sort(unique(dat_pos$Year))

liste_cartes <- list()


for (an in annees) {
  # filtrer l'année
  pred_year <- pred_grid %>%
    filter(Year == an)
  
  # carte année courante
  carte <- ggplot(pred_year, aes(x = X, y = Y, fill = exp(CPUE_pred))) +
    geom_tile() +
    scale_fill_viridis_c(option = "viridis", name = "Biomasse prédite") +
    coord_equal() +
    labs(
      title = paste("Biomasse prédite - année", an),
      x = "X (km)",
      y = "Y (km)"
    ) +
    theme_minimal()
  
  if (SavePlots) {
    ggsave(
      filename = file.path("data/", output_plot, paste0("map_Biomasse_pred_m_1_", an, ".png")),
      plot = carte,
      width = 8,
      height = 6,
      dpi = 300
    )}
}
# 
# mois <- sort(unique(dat_pos$Month))
# 
# liste_cartes <- list()
# 
# 
# for (mo in mois) {
#   # filtrer l'année
#   pred_year <- pred_grid %>%
#     filter(Month == mo)
#   
#   # carte année courante
#   carte <- ggplot(pred_year, aes(x = X, y = Y, fill = exp(CPUE_pred))) +
#     geom_tile() +
#     scale_fill_viridis_c(option = "viridis", name = "Biomasse prédite") +
#     coord_equal() +
#     labs(
#       title = paste("Biomasse prédite - Mois", an),
#       x = "X (km)",
#       y = "Y (km)"
#     ) +
#     theme_minimal()
#   
#   if (SavePlots) {
#     ggsave(
#       filename = file.path("data/", output_plot, paste0("map_Biomasse_pred_mois_m_1_", mo, ".png")),
#       plot = carte,
#       width = 8,
#       height = 6,
#       dpi = 300
#     )}
# }




# Effet spatial
grid_spatial <- grid %>%
  mutate(Year = 2015)
pred_spatial <- predict(m_1, newdata = grid_spatial)

grid_spatial$effet_spatial <- pred_spatial$omega_s

lim <- max(abs(grid_spatial$effet_spatial), na.rm = TRUE)

Spat_effect <- ggplot() +
  
  # fond gris sur toute l'emprise
  annotate(
    "rect",
    xmin = min(grid_spatial$X),
    xmax = max(grid_spatial$X),
    ymin = min(grid_spatial$Y),
    ymax = max(grid_spatial$Y),
    fill = "grey80",
    alpha = 0.4
  ) +
  
  # raster prédit
  geom_raster(
    data = grid_spatial,
    aes(X, Y, fill = effet_spatial)
  ) +
  
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-lim, lim),
    name = "Effet spatial"
  ) +
  
  coord_equal() +
  theme_minimal()

if (SavePlots) {
  ggsave(
    filename = file.path("data/", output_plot, "Effet_spatial_m_1.png"),
    plot     = Spat_effect,
    width    = 8,
    height   = 6,
    dpi      = 300
  )
}


# Effet Spatio-temporel

pred_spatiotemp <- predict(m_1, newdata = grid_all_years)

grid_all_years$effet_spatiotemp <- pred_spatiotemp$epsilon_st


annees <- sort(unique(grid_all_years$Year))

liste_cartes_spatiotemp <- list()

lim <- max(abs(grid_all_years$effet_spatiotemp), na.rm = TRUE)

for (an in annees) {
  
  pred_year <- grid_all_years %>%
    filter(Year == an)
  
  carte <- ggplot() +
    
    annotate(
      "rect",
      xmin = min(grid_all_years$X),
      xmax = max(grid_all_years$X),
      ymin = min(grid_all_years$Y),
      ymax = max(grid_all_years$Y),
      fill = "grey90",
      alpha = 0.4
    ) +
    
    geom_raster(
      data = pred_year,
      aes(X, Y, fill = effet_spatiotemp)
    ) +
    
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      limits = c(-lim, lim),
      name = "Effet spatio-temporel"
    ) +
    
    coord_equal() +
    
    labs(
      title = paste("Effet spatio-temporel - année", an),
      x = "X (km)",
      y = "Y (km)"
    ) +
    
    theme_minimal()
  
  if (SavePlots) {
    ggsave(
      filename = file.path("data/", output_plot, paste0("map_Effet_Spatiotemp_m_1_", an, ".png")),
      plot = carte,
      width = 8,
      height = 6,
      dpi = 300
    )}
}


pred_fix <- predict(
  m_1,
  newdata = grid_all_years,
  re_form = NA
)

grid_all_years$effet_fixe <- pred_fix$est
unique(grid_all_years$effet_fixe)


#Effet année
pred_fix <- predict(
  m_spatemp,
  newdata = grid_all_years,
  re_form = NA
)

grid_all_years$effet_annee <- pred_fix$est

year_effect <- grid_all_years %>%
  group_by(Year) %>%
  summarise(effet = mean(effet_annee))

Effet_An <- ggplot(year_effect, aes(Year, effet)) +
  geom_line() +
  geom_point() +
  theme_minimal()


if (SavePlots) {
  ggsave(
    filename = file.path("data/", output_plot, "Effet_annuel_m_1.png"),
    plot     = Effet_An,
    width    = 8,
    height   = 6,
    dpi      = 300
  )
}





#Indice d'abondance

grid_small <- grid_with_id 
grid_yrs <- tidyr::crossing(
  grid_small,
  Year = unique(dat_pos$Year))
# ) %>%
#   mutate(
#     nb_ham_std = median(dat_pos$nb_ham_std, na.rm = TRUE),
#     trempe_std = median(dat_pos$trempe_std, na.rm = TRUE)
#   )

p_st <- predict(m_1, newdata = grid_yrs, 
                return_tmb_object = TRUE)
index <- get_index(
  p_st,
  area = rep(1, nrow(grid_yrs))
)


Abundance_Index <- ggplot(index, aes(Year, est)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey90") +
  geom_line(linewidth = 1, colour = "grey30") +
  labs(x = "Year", y = "Biomass") +
  theme_bw()

if (SavePlots) {
  ggsave(
    filename = file.path("data/", output_plot, "Indice_abondance_m_1.png"),
    plot     = Abundance_Index,
    width    = 8,
    height   = 6,
    dpi      = 300
  )
}


#Visualisation des effets
# Effet marginal de Year, Month marginalisé
# prédiction
Year_effect <- ggpredict(
  m_1,
  terms = "Year [all]",
  condition = c(nb_ham = 1)
)

# créer un objet ggplot
Year_plot <- plot(Year_effect) + theme_bw()

# afficher
Year_plot

# sauvegarder
if (SavePlots) {
  ggsave(
    filename = file.path("data/", output_plot, "Effet_annuel_m_1.png"),
    plot     = Year_plot,
    width    = 8,
    height   = 6,
    dpi      = 300
  )
}


# # prédiction
# ham_effect <- ggpredict(
#   m_1,
#   terms = c("nb_ham_std [all]", "Year [all]")
# )
# 
# meanham_effect <- ham_effect %>%
#   group_by(x) %>%
#   summarise(
#     predicted = mean(predicted),
#     conf.low = mean(conf.low),
#     conf.high = mean(conf.high)
#   )
# 
# ham_plot <- ggplot(meanham_effect, aes(x = x, y = predicted)) +
#   geom_line(size = 1.2, color = "black") +
#   geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
#               alpha = 0.1) +
#   theme_bw() 
# 
# ham_plot
# 
# # sauvegarder
# if (SavePlots) {
#   ggsave(
#     filename = file.path("data/", output_plot, "Effet_ham_m_1.png"),
#     plot     = ham_plot,
#     width    = 8,
#     height   = 6,
#     dpi      = 300
#   )
# }
# 
# 
# 
# trempe_effect <- ggpredict(m_1, terms = c("trempe_std [all]", "Year [all]"))
# 
# mean_effect <- trempe_effect %>%
#   group_by(x) %>%
#   summarise(
#     predicted = mean(predicted),
#     conf.low = mean(conf.low),
#     conf.high = mean(conf.high)
#   )
# 
# # créer un objet ggplot
# trempe_plot <- ggplot(mean_effect, aes(x = x, y = predicted)) +
#   geom_line(size = 1.2, color = "black") +
#   geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
#               alpha = 0.1) +
#   theme_bw() 
# 
# # afficher
# trempe_plot
# 
# # sauvegarder
# if (SavePlots) {
#   ggsave(
#     filename = file.path("data/", output_plot, "Effet_trempe_m_1.png"),
#     plot     = trempe_plot,
#     width    = 8,
#     height   = 6,
#     dpi      = 300
#   )
# }








