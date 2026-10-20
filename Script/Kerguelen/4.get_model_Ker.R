
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

dat_pos <- dat_pos %>%
  filter(HookLength != 140)

mean_value <- mean(dat_pos$HookLength, na.rm = TRUE)


dat_pos$HookLength[is.na(dat_pos$HookLength)] <- mean_value

dat_pos <- dat_pos %>%
  filter(!is.na(bottomT))




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
    NB_MEAN_ORQUE_sc = as.numeric(scale(NB_MEAN_ORQUE)),
    bottomT_sc = as.numeric(scale(bottomT))
  )

dat_pos_std <- dat_pos_std %>% 
  mutate(
    cachalot_obs = factor(cachalot_obs,
                          levels = c(0, 1),
                          labels = c("Absent", "Présent")),
    orque_obs = factor(orque_obs,
                          levels = c(0, 1),
                          labels = c("Absent", "Présent")),
    
  )
# 4. Mesh -----------------------------------------------------------------

Study_domain_Hole <- readRDS("Data/Ker_polygon.rds")%>%
  mutate(geometry = st_geometry(.) / 1000)

# Update CRS to indicate units in kilometers (custom CRS)

st_crs(Study_domain_Hole) <- "+proj=utm +zone=42 +datum=WGS84 +units=km +no_defs"


mesh_pos <- make_mesh(
  dat_pos,
  xy_cols = c("X", "Y"),
  cutoff = 30,
  max_edge = c(20, 100)   # max_edge_in = 20, max_edge_out = 100
)


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
    filename = file.path("data/", output_plot, "mesh_plot_ker.png"),
    plot     = mesh_plot,
    width    = 8,
    height   = 6,
    dpi      = 300
  )
}
cat("→ mesh.png sauvegardé\n\n")


m_1 <- sdmTMB(
  ScanWeightKg ~ s(TAAFseason_sc) + s(Month, bs = "cc") + Vessel + s(prof_bathy_sc) + s(slope_sc) + s(HookLength_sc)
   + s(NB_MEAN_CACHALOT_sc),
  offset = log(
    0.56 * 0.002 * (nb_ham_med - 1) + 0.246
  ),
  data = dat_pos_std,
  mesh = bspde_pos,
  family = Gamma(link = "log"),
  spatial = "on",
  time= "TAAFseason",
  spatiotemporal = "iid"
)


saveRDS(
  m_1,
  file = file.path("data", output_plot, "m_1_ker.rds")
)

cat("→ m_1.rds sauvegardé\n\n")

m_2 <- sdmTMB(
  ScanWeightKg ~ s(TAAFseason_sc) + s(Month, bs = "cc") + Vessel + s(prof_bathy_sc) + s(slope_sc) + s(HookLength_sc),
  offset = log(
    0.56 * 0.002 * (nb_ham_med - 1) + 0.246
  ),
  data = dat_pos_std,
  mesh = bspde_pos,
  family = Gamma(link = "log"),
  spatial = "on",
  time = "TAAFseason",
  spatiotemporal = "iid"
)


saveRDS(
  m_2,
  file = file.path("data", output_plot, "m_2_ker.rds")
)

cat("→ m_2.rds sauvegardé\n\n")



# Prédiction sur la donnée
# 1. Jeux de données pour les scénarios
# -------------------------------------

# Situation réelle (avec déprédation observée)
dat_reel <- dat_pos_std

# Sans aucune déprédation
dat_sans_dep <- dat_pos_std %>%
  mutate(
    NB_MEAN_CACHALOT = 0,
    NB_MEAN_ORQUE = 0
  )

# Sans orques uniquement
dat_sans_orque <- dat_pos_std %>%
  mutate(
    NB_MEAN_ORQUE = 0
  )

# Sans cachalots uniquement
dat_sans_cachalot <- dat_pos_std %>%
  mutate(
    NB_MEAN_CACHALOT = 0
  )

# 2. Offsets
# ----------
offset_reel <- log(0.46*0.002*dat_pos_std$nb_ham + 0.246)
offset_sans_dep <- log(0.46*0.002*dat_pos_std$nb_ham + 0.246)
offset_sans_orque <- log(0.46*0.002*dat_pos_std$nb_ham + 0.246)
offset_sans_cachalot <- log(0.46*0.002*dat_pos_std$nb_ham + 0.246)

# 3. Prédictions du modèle
# ------------------------

pred_reel <- predict(
  m_1,
  newdata = dat_reel,
  offset = offset_reel,  # si on veut une CPUE en sortie
  type = "response"
)

pred_sans_dep <- predict(
  m_1,
  newdata = dat_sans_dep,
  offset = offset_sans_dep,
  type = "response"
)


# 4. Dataframe final
# ------------------

df_station_predictions <- dat_pos_std %>%
  mutate(
    poids_reel = ScanWeightKg,
    
    pred_modele = pred_reel$est,
    
    pred_sans_dep = pred_sans_dep$est,
    
    perte_totale = pred_sans_dep - pred_modele,
    
  ) %>%
  select(
    StationID,
    TAAFseason,
    Month,
    ScanWeightKg,
    nb_ham,
    NB_MEAN_CACHALOT,
    X,
    Y,
    poids_reel,
    pred_modele,
    perte_totale
  )


saveRDS(
  df_station_predictions,
  file = file.path("data", output_plot, "df_station_predictions_ker.rds")
)

cat("→ df_station_predictions.rds sauvegardé\n\n")


get_mean_bottomT_season <- function(lon_sta, lat_sta,
                                    season,
                                    dossier_nc){
  
  jours <- seq(
    as.Date(paste0(season, "-09-01")),
    as.Date(paste0(season + 1, "-09-30")),
    by = "day"
  )
  
  valeurs <- numeric(length(jours))
  
  for(k in seq_along(jours)){
    
    d <- jours[k]
    
    fichier <- file.path(
      dossier_nc,
      sprintf(
        "bottomT_daily_%s_%s.nc",
        format(d, "%Y"),
        format(d, "%m")
      )
    )
    
    if(!file.exists(fichier)){
      valeurs[k] <- NA
      next
    }
    
    dat <- get_nc_data(fichier)
    
    i <- which.min(abs(dat$lon - lon_sta))
    j <- which.min(abs(dat$lat - lat_sta))
    
    jour_mois <- as.integer(format(d, "%d"))
    
    valeurs[k] <- dat$temp[i, j, jour_mois]
  }
  
  mean(valeurs, na.rm = TRUE)
}

#Prédiction sur une grille
#1. Fabrication d'une grille d'extrapolation ------------------------------------------------------
# Import et préparation de la grille d’extrapolation spatiale
Extrapol_grid <- readRDS(here("data/outputs/df/grid_extrapol.rds"))

# Mise à l’échelle de la profondeur et de la pente sur la grille
Mean_SD <- dat_pos %>% ungroup() %>% summarise(Mean = mean(prof_bathy), SD = sd(prof_bathy))
Extrapol_grid = Extrapol_grid %>% rename(depth = profondeur.pred)
Extrapol_grid <- Extrapol_grid %>% mutate(prof_bathy_sc = (depth - Mean_SD$Mean) / Mean_SD$SD)

Mean_SD <- dat_pos %>% ungroup() %>% summarise(Mean = mean(slope), SD = sd(slope))
Extrapol_grid <- Extrapol_grid %>% mutate(slope_sc = (slope - Mean_SD$Mean) / Mean_SD$SD)

# Nettoyage et sélection des colonnes nécessaires pour la prédiction
grid <- Extrapol_grid %>%
  mutate(
    nb_ham = nb_ham,
    tps_h = mean(dat$tps_h),
    x = x / 1000 ,
    y = y / 1000 
  ) %>%
  filter(!is.na(slope_sc)) %>%
  filter(!is.na(prof_bathy_sc)) %>%
  filter(depth > min_depth) %>%
  rename(X = x, Y = y) %>%
  dplyr::select(prof_bathy_sc, slope_sc, slope, depth, nb_ham, tps_h, X, Y) %>%
  rename(prof_bathy = depth)


#2. Prédictions -----------------------------------

# Vérifications rapides
stopifnot("TAAFseason" %in% names(dat_pos_std))
stopifnot(is.numeric(dat_pos_std$TAAFseason))

#Préparer grid avec un identifiant de cellule
grid_with_id <- grid %>%
  mutate(cell_id = row_number())

grid_with_id <- dplyr::as_tibble(grid_with_id)

#Réplication temporelle : une ligne par cellule × année
years <- sort(unique(dat_pos$TAAFseason))
years <- years[years <= 2025]
grid_all_years <- grid_with_id %>%
  tidyr::crossing(TAAFseason = years) %>%
  mutate(
    cachalot_obs = 1,
    Month = 1,
    orque_obs = 1,
    nb_ham = median(dat_pos_std$nb_ham, na.rm = TRUE),
    NB_MEAN_CACHALOT = 0,
    NB_MEAN_ORQUE = 0,
    Vessel = factor("Albius", levels = levels(dat_pos$Vessel)),
    HookLength = mean(dat_pos$HookLength)
  ) %>%
  relocate(cell_id)


grid_all_years <- grid_all_years %>%
  mutate(
    cachalot_obs = factor(
      "1",
      levels = levels(dat_pos_std$cachalot_obs)
    ),
    orque_obs = factor(
      "1",
      levels = levels(dat_pos_std$orque_obs)
    )
  )
Mean_SD_hook <- dat_pos %>%
  summarise(Mean = mean(HookLength), SD = sd(HookLength))

grid_all_years <- grid_all_years %>% mutate(HookLength_sc = (HookLength - Mean_SD_hook$Mean) / Mean_SD_hook$SD)

Mean_SD_cach <- dat_pos %>%
  summarise(Mean = mean(NB_MEAN_CACHALOT), SD = sd(NB_MEAN_CACHALOT))

grid_all_years <- grid_all_years %>% mutate(NB_MEAN_CACHALOT_sc = (NB_MEAN_CACHALOT - Mean_SD_cach$Mean) / Mean_SD_cach$SD)

#Vérifier colonnes requises
required_cols <- c("prof_bathy", "TAAFseason", "X", "Y")
missing_cols <- setdiff(required_cols, names(grid_all_years))
if (length(missing_cols) > 0) stop("Colonnes manquantes dans newdata : ", paste(missing_cols, collapse = ", "))


#3. Prédiction sur m_1
pred_std <- predict(
  m_1,
  newdata = grid_all_years,
  offset = log(0.46*0.002*grid_all_years$nb_ham + 0.246),
  return_tmb_object = TRUE
)

index_std <- get_index(pred_std, bias_correct = TRUE)


#4. Prédiction sur m_2
pred <- predict(
  m_2,
  newdata = grid_all_years,
  offset = log(0.46*0.002*grid_all_years$nb_ham + 0.246),
  return_tmb_object = TRUE
)

index <- get_index(pred, bias_correct = TRUE)

#5. Sauvegarde
saveRDS(
  pred_std,
  file = file.path("data", output_plot, "pred_depred_ker.rds")
)

saveRDS(
  index_std,
  file = file.path("data", output_plot, "index_depred_ker.rds")
)


saveRDS(
  pred,
  file = file.path("data", output_plot, "pred_ker.rds")
)

saveRDS(
  index,
  file = file.path("data", output_plot, "index_ker.rds")
)

cat("→ Prédictions et indices sauvegardés\n\n")
