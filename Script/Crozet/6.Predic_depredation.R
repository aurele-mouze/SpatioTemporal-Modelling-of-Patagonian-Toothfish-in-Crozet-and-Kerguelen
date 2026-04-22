
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
      !is.na(NB_MIN_CACHALOT) & !is.na(NB_MAX_CACHALOT) ~ (NB_MIN_CACHALOT + NB_MAX_CACHALOT) / 2,
      is.na(NB_MIN_CACHALOT) & !is.na(NB_MAX_CACHALOT) ~ NB_MAX_CACHALOT,
      !is.na(NB_MIN_CACHALOT) & is.na(NB_MAX_CACHALOT) ~ NB_MIN_CACHALOT,
      TRUE ~ NA_real_
    )
  )

dat_pos <- dat_pos %>%
  mutate(
    NB_MEAN_ORQUE = case_when(
      !is.na(NB_MIN_ORQUE) & !is.na(NB_MAX_ORQUE) ~ (NB_MIN_ORQUE + NB_MAX_ORQUE) / 2,
      is.na(NB_MIN_ORQUE) & !is.na(NB_MAX_ORQUE) ~ NB_MAX_ORQUE,
      !is.na(NB_MIN_ORQUE) & is.na(NB_MAX_ORQUE) ~ NB_MIN_ORQUE,
      TRUE ~ NA_real_
    )
  )

dat_pos <- dat_pos %>%
  mutate(
    # indicateur d'observation
    cachalot_obs = as.integer(PRES_CACHALOT != "non obsé"),
    orque_obs = as.integer(PRES_ORQUE != "non obsé"),
    
    # remplacer NA uniquement quand non observé
    NB_MEAN_CACHALOT = ifelse(PRES_CACHALOT == "non obsé", 0, NB_MEAN_CACHALOT),
    NB_MEAN_ORQUE = ifelse(PRES_ORQUE == "non obsé", 0, NB_MEAN_ORQUE)
  )


dat_pos <- dat_pos %>%
  filter(!Year %in% c(2002, 2025))

dat_pos <- dat_pos %>%
  filter(
    !is.na(PRES_CACHALOT),
    !is.na(PRES_ORQUE)
  )

dat_pos <- dat_pos %>%
  filter(
    !is.na(NB_MEAN_CACHALOT),
    !is.na(NB_MEAN_ORQUE)
  )

# 4. Mesh -----------------------------------------------------------------

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
  crs = st_crs(Study_domain_Hole) 
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




m_1 <- sdmTMB(
  ScanWeightKg ~ s(Year) + s(NB_MEAN_CACHALOT) + s(NB_MEAN_ORQUE) + cachalot_obs + orque_obs,
  offset = log(dat_pos$nb_ham),
  data = dat_pos,
  mesh = bspde_pos,
  family = Gamma(link = "log"),
  spatial = "on",
  time= "Year",
  spatiotemporal = "iid"
)



#5. Contrôle des modèles ---------------------------------

dat_pos$resids <- residuals(m_1)
hist_residual <- ggplot(dat_pos, aes(x = resids)) +
  geom_histogram(
    bins = 30,
    fill = "steelblue",
    color = "black",
    alpha = 0.7
  ) +
  labs(title = "Histogramme des résidus", x = "Résidus", y = "Fréquence") +
  theme_bw()

s <- simulate(m_1, nsim = number_sim, type = "mle-mvn")
dharma_qqplot = dharma_residuals(s, m_1, test_uniformity = TRUE)

QQplot <- ggplot(dharma_qqplot, aes(x = expected, y = observed)) +
  geom_point(alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0) +
  theme_minimal()

if (SavePlots) {
  ggsave(
    filename = file.path("data/", output_plot, "residual-histogramme-m_1.png"),
    plot     = hist_residual,
    width    = 8,
    height   = 6,
    dpi      = 300
  )
}

if (SavePlots) {
  ggsave(
    filename = file.path("data/", output_plot, "qqplot_residual-m_1.png"),
    plot     = QQplot,
    width    = 8,
    height   = 6,
    dpi      = 300
  )
}



#6. Fabrication d'une grille d'extrapolation ------------------------------------------------------
# Import et préparation de la grille d’extrapolation spatiale
Extrapol_grid <- readRDS(here("data/outputs/df/grid_extrapolCRO.rds"))

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
    x = x / 1000 ,
    y = y / 1000 
  ) %>%
  filter(!is.na(slope_scaled)) %>%
  filter(!is.na(depth_scaled)) %>%
  filter(depth > min_depth) %>%
  rename(X = x, Y = y) %>%
  dplyr::select(depth_scaled, slope_scaled, slope, depth, nb_ham, tps_h, X, Y) %>%
  rename(prof_bathy = depth)


#7. Prédictions -----------------------------------

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
  mutate(cachalot_obs = 1,
         orque_obs = 1,
         nb_ham = mean(dat_pos$nb_ham, na.rm = TRUE),
         NB_MEAN_CACHALOT = mean(dat_pos$NB_MEAN_CACHALOT, na.rm = TRUE),
         NB_MEAN_ORQUE = mean(dat_pos$NB_MEAN_ORQUE, na.rm = TRUE)) %>%
  relocate(cell_id)

#Vérifier colonnes requises
required_cols <- c("prof_bathy", "Year", "X", "Y")
missing_cols <- setdiff(required_cols, names(grid_all_years))
if (length(missing_cols) > 0) stop("Colonnes manquantes dans newdata : ", paste(missing_cols, collapse = ", "))

# Prédiction (utilise chunking si nécessaire)
pred_grid_raw <- predict_in_chunks(m_1, grid_all_years, chunk_size = 20000)



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



# prédiction
PresCach_effect <- ggpredict(
  m_1,
  terms = c("NB_MEAN_CACHALOT [all]", "Year [all]")
)

meanPresCach_effect <- PresCach_effect %>%
  group_by(x) %>%
  summarise(
    predicted = mean(predicted),
    conf.low = mean(conf.low),
    conf.high = mean(conf.high)
  )

PresCach_plot <- ggplot(meanPresCach_effect, aes(x = x, y = predicted)) +
  geom_line(size = 1.2, color = "black") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.1) +
  theme_bw()


Cach_plot <- PresCach_plot +
  geom_rug(
    data = dat_pos,
    aes(x = NB_MEAN_CACHALOT),
    inherit.aes = FALSE,
    alpha = 0.2
  )

# sauvegarder
if (SavePlots) {
  ggsave(
    filename = file.path("data/", output_plot, "Effet_Cach_m_1.png"),
    plot     = Cach_plot,
    width    = 8,
    height   = 6,
    dpi      = 300
  )
}


Presorque_effect <- ggpredict(
  m_1,
  terms = c("NB_MEAN_ORQUE [all]", "Year [all]")
)

meanPresorque_effect <- Presorque_effect %>%
  group_by(x) %>%
  summarise(
    predicted = mean(predicted),
    conf.low = mean(conf.low),
    conf.high = mean(conf.high)
  )

Presorque_plot <- ggplot(meanPresorque_effect, aes(x = x, y = predicted)) +
  geom_line(size = 1.2, color = "black") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.1) +
  theme_bw()


Orque_plot <- Presorque_plot +
  geom_rug(
    data = dat_pos,
    aes(x = NB_MEAN_ORQUE),
    inherit.aes = FALSE,
    alpha = 0.2
  )

# sauvegarder
if (SavePlots) {
  ggsave(
    filename = file.path("data/", output_plot, "Effet_Orque_m_1.png"),
    plot     = Orque_plot,
    width    = 8,
    height   = 6,
    dpi      = 300
  )
}



#Indice d'abondance



p_st <- predict(m_1, newdata = grid_all_years, 
                return_tmb_object = TRUE)
index <- get_index(
  p_st,
  area = rep(1, nrow(grid_all_years))
)

index_mean <- index %>%
  mutate(est_mean = est / sum(rep(4, nrow(grid_all_years)) / length(unique(grid_all_years$Year))))

Abundance_Index <- ggplot(index_mean, aes(Year, est)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey90") +
  geom_line(linewidth = 1, colour = "grey30") +
  labs(x = "Year", y = "Biomasse (g/hameçon)") +
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



