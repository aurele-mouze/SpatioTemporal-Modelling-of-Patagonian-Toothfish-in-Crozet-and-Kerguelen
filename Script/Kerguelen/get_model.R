
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
    !is.na(pres_cachalot),
    !is.na(pres_orque)
  )

dat_pos <- dat_pos %>%
  filter(
    !is.na(NB_MEAN_CACHALOT),
    !is.na(NB_MEAN_ORQUE),
    !is.na(bottomT)
  )

dat_pos$Vessel <- droplevels(dat_pos$Vessel)


dat_pos$HookLength <- as.numeric(sub("_.*", "", dat_pos$HookSize))


dat_pos <- dat_pos %>%
  filter(HookLength != 7)


mean_value <- mean(dat_pos$HookLength, na.rm = TRUE)


dat_pos$HookLength[is.na(dat_pos$HookLength)] <- mean_value



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
                       labels = c("Absent", "Présent"))
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
  ScanWeightKg ~ s(TAAFseason) + s(Month, bs = "cc") + Vessel + s(prof_bathy_sc) + s(slope_sc) + s(HookLength_sc) + s(bottomT_sc)
   + s(NB_MEAN_CACHALOT_sc),
  offset = log(dat_pos_std$nb_ham*0.246),
  data = dat_pos_std,
  mesh = bspde_pos,
  family = Gamma(link = "log"),
  spatial = "on",
  time= "Year",
  spatiotemporal = "iid"
)

m_2 <- sdmTMB(
  ScanWeightKg ~ s(TAAFseason) + s(Month, bs = "cc") + Vessel + s(prof_bathy_sc) + s(slope_sc) + s(HookLength_sc)
  + s(NB_MEAN_CACHALOT_sc) + cachalot_obs,
  offset = log(dat_pos_std$nb_ham*0.246),
  data = dat_pos_std,
  mesh = bspde_pos,
  family = Gamma(link = "log"),
  spatial = "on",
  time= "Year",
  spatiotemporal = "iid"
)

m_3 <- sdmTMB(
  ScanWeightKg ~ s(TAAFseason) + s(Month, bs = "cc") + Vessel + s(prof_bathy_sc) + s(HookLength_sc)
  + s(NB_MEAN_CACHALOT_sc),
  offset = log(dat_pos_std$nb_ham*0.246),
  data = dat_pos_std,
  mesh = bspde_pos,
  family = Gamma(link = "log"),
  spatial = "on",
  time= "Year",
  spatiotemporal = "iid"
)

m_4 <- sdmTMB(
  ScanWeightKg ~ s(TAAFseason) + s(Month, bs = "cc") + Vessel + s(HookLength_sc)
  + s(NB_MEAN_CACHALOT_sc),
  offset = log(dat_pos_std$nb_ham*0.246),
  data = dat_pos_std,
  mesh = bspde_pos,
  family = Gamma(link = "log"),
  spatial = "on",
  time= "Year",
  spatiotemporal = "iid"
)

m_5 <- sdmTMB(
  ScanWeightKg ~ s(TAAFseason) + s(Month, bs = "cc") + Vessel 
  + s(NB_MEAN_CACHALOT_sc),
  offset = log(dat_pos_std$nb_ham*0.246),
  data = dat_pos_std,
  mesh = bspde_pos,
  family = Gamma(link = "log"),
  spatial = "on",
  time= "Year",
  spatiotemporal = "iid"
)

m_6 <- sdmTMB(
  ScanWeightKg ~ s(TAAFseason) + s(Month, bs = "cc") 
  + s(NB_MEAN_CACHALOT_sc),
  offset = log(dat_pos_std$nb_ham*0.246),
  data = dat_pos_std,
  mesh = bspde_pos,
  family = Gamma(link = "log"),
  spatial = "on",
  time= "Year",
  spatiotemporal = "iid"
)

m_7 <- sdmTMB(
  ScanWeightKg ~ s(TAAFseason) 
  + s(NB_MEAN_CACHALOT_sc),
  offset = log(dat_pos_std$nb_ham*0.246),
  data = dat_pos_std,
  mesh = bspde_pos,
  family = Gamma(link = "log"),
  spatial = "on",
  time= "Year",
  spatiotemporal = "iid"
)

m_8 <- sdmTMB(
  ScanWeightKg ~ s(TAAFseason) ,
  offset = log(dat_pos_std$nb_ham*0.246),
  data = dat_pos_std,
  mesh = bspde_pos,
  family = Gamma(link = "log"),
  spatial = "on",
  time= "Year",
  spatiotemporal = "iid"
)

m_9 <- sdmTMB(
  ScanWeightKg ~ 1,
  offset = log(dat_pos_std$nb_ham*0.246),
  data = dat_pos_std,
  mesh = bspde_pos,
  family = Gamma(link = "log"),
  spatial = "on",
  time= "Year",
  spatiotemporal = "iid"
)

m_10 <- sdmTMB(
  ScanWeightKg ~ s(TAAFseason) + s(Month, bs = "cc") + cachalot_obs 
  + s(NB_MEAN_CACHALOT_sc),
  offset = log(dat_pos_std$nb_ham*0.246),
  data = dat_pos_std,
  mesh = bspde_pos,
  family = Gamma(link = "log"),
  spatial = "on",
  time= "Year",
  spatiotemporal = "iid"
)

AIC(m_1, m_2,m_3, m_4, m_5, m_6, m_7, m_8, m_9)


AIC(m_2)
saveRDS(
  m_1,
  file = file.path("data", output_plot, "m_1_ker.rds")
)

cat("→ m_1.rds sauvegardé\n\n")

# m_2 <- sdmTMB(
#   ScanWeightKg ~ s(TAAFseason) + s(Month, bs = "cc") + Vessel + s(prof_bathy) + s(slope) + s(HookLength) + s(bottomT),
#   offset = log(dat_pos$nb_ham*0.246),
#   data = dat_pos,
#   mesh = bspde_pos,
#   family = Gamma(link = "log"),
#   spatial = "on",
#   time = "Year",
#   spatiotemporal = "iid"
# )
# 
# 
# saveRDS(
#   m_2,
#   file = file.path("data", output_plot, "m_2_ker.rds")
# )
# 
# cat("→ m_2.rds sauvegardé\n\n")
# 



# 2. Jeux de données pour les scénarios
# -------------------------------------
mean_cach <- mean(dat_pos$NB_MEAN_CACHALOT)
sd_cach   <- sd(dat_pos$NB_MEAN_CACHALOT)

# Situation réelle (avec déprédation observée)
dat_reel <- dat_pos_std

# Sans aucune déprédation
dat_sans_dep <- dat_pos_std %>%
  mutate(
    NB_MEAN_CACHALOT = 0,
    NB_MEAN_CACHALOT_sc = (0 - mean_cach) / sd_cach
  )


# 3. Offsets
# ----------
offset_reel <- log(0.246 * dat_reel$nb_ham)
offset_sans_dep <- log(0.246 * dat_sans_dep$nb_ham)

# 4. Prédictions du modèle
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


# 5. Dataframe final
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
    NB_MEAN_CACHALOT_sc,
    X,
    Y,
    poids_reel,
    pred_modele,
    perte_totale,
    pred_sans_dep,
    pred_modele
  )


saveRDS(
  df_station_predictions,
  file = file.path("data", output_plot, "df_station_predictions_ker.rds")
)

cat("→ df_station_predictions.rds sauvegardé\n\n")

df_boxplot_pct <- df_station_predictions %>%
  filter(TAAFseason != 2026) %>%
  mutate(
    perte_totale_pct = 100 * perte_totale / (pred_modele + perte_totale)
  ) %>%
  transmute(
    TAAFseason,
    type_perte = "Déprédation totale",
    pourcentage_perte = perte_totale_pct
  )

df_summary <- df_boxplot_pct %>%
  group_by(TAAFseason, type_perte) %>%
  summarise(
    mean_perte = mean(pourcentage_perte, na.rm = TRUE),
    sd = sd(pourcentage_perte, na.rm = TRUE),
    n = sum(!is.na(pourcentage_perte)),
    se = sd / sqrt(n),
    ic_inf = mean_perte - 1.96 * se,
    ic_sup = mean_perte + 1.96 * se,
    .groups = "drop"
  )

ggplot(df_summary,
       aes(x = TAAFseason,
           y = mean_perte,
           color = type_perte)) +
  
  geom_line(linewidth = 1.1) +
  
  geom_ribbon(
    aes(ymin = ic_inf, ymax = ic_sup, fill = type_perte),
    alpha = 0.15,
    color = NA
  ) +
  
  theme_minimal() +
  labs(
    x = "Year",
    y = "Average loss (%)",
    color = "script",
    fill = "script"
  )



crs_ker_km <- "+proj=utm +zone=42 +south +datum=WGS84 +units=km +no_defs"


stations_sf <- st_as_sf(
  df_station_predictions,
  coords = c("X", "Y"),
  crs = st_crs(crs_ker_km)
)

bbox <- st_bbox(stations_sf)

grid_4km <- st_make_grid(
  st_as_sfc(st_bbox(stations_sf)),
  cellsize = 7,
  square = TRUE
)

grid_4km <- st_sf(
  cell_id = seq_along(grid_4km),
  geometry = grid_4km
)



stations_grid <- st_join(
  stations_sf,
  grid_4km,
  join = st_within
)

grid_summary <- stations_grid %>%
  st_drop_geometry() %>%
  group_by(cell_id) %>%
  summarise(
    poids_reel = sum(pred_modele, na.rm = TRUE),
    poids_sans_dep = sum(pred_sans_dep, na.rm = TRUE),
    n_station = n()
  ) %>%
  mutate(
    perte_totale_pct =
      100 * (poids_sans_dep - poids_reel) / poids_sans_dep
  )

grid_plot <- left_join(
  grid_4km,
  grid_summary,
  by = "cell_id"
)

grid_plot <- grid_plot %>%
  filter(!is.na(perte_totale_pct))

ggplot(grid_plot) +
  geom_sf(aes(fill = perte_totale_pct), color = NA)  +
  scale_fill_viridis_c() +
  theme_bw()

grid_plot_4326 <- st_transform(grid_plot, 4326)
Ker_map_dep <- KER_map +
  geom_sf(
    data = grid_plot_4326,
    aes(fill = perte_totale_pct),
    color = NA,
    alpha = 0.7,
    inherit.aes = FALSE
  ) +
  scale_fill_viridis_c(
    name = "Perte (%)"
  )+
  annotation_scale(
    location = "bl",
    width_hint = 0.3
  ) +
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    style = north_arrow_fancy_orienteering
  ) +
  theme_bw() +
  theme(
    legend.position = c(0.95, 0.10)
  )

Ker_map_dep

stations_grid <- stations_grid %>%
  filter(TAAFseason != 2026) %>%
  mutate(
    periode = case_when(
      TAAFseason >= 2008 & TAAFseason <= 2013 ~ "2008-2013",
      TAAFseason >= 2014 & TAAFseason <= 2019 ~ "2014-2019",
      TAAFseason >= 2020 & TAAFseason <= 2025 ~ "2020-2025"
    )
  )

grid_summary_year <- stations_grid %>%
  st_drop_geometry() %>%
  group_by(periode, TAAFseason, cell_id) %>%
  summarise(
    poids_reel = sum(pred_modele, na.rm = TRUE),
    poids_sans_dep = sum(pred_sans_dep, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    perte_totale_pct =
      100 * (poids_sans_dep - poids_reel) / poids_sans_dep
  )

grid_summary_period <- grid_summary_year %>%
  group_by(periode, cell_id) %>%
  summarise(
    perte_moy = mean(perte_totale_pct, na.rm = TRUE),
    .groups = "drop"
  )

grid_plot_period <- grid_4km %>%
  right_join(
    grid_summary_period,
    by = "cell_id"
  )

grid_plot_period_4326 <- st_transform(
  grid_plot_period,
  4326
)

Ker_map_dep_period <- KER_map +
  geom_sf(
    data = grid_plot_period_4326,
    aes(fill = perte_moy),
    color = NA,
    alpha = 0.7,
    inherit.aes = FALSE
  ) +
  facet_wrap(~periode) +
  scale_fill_viridis_c(
    name = "Perte (%)"
  ) +
  annotation_scale(
    location = "bl",
    width_hint = 0.3
  ) +
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    style = north_arrow_fancy_orienteering
  ) +
  theme_bw(legend.position = c(0.85, 0.15))

Ker_map_dep_period




grid_summary_year <- stations_grid %>%
  st_drop_geometry() %>%
  group_by(cell_id, TAAFseason) %>%
  summarise(
    poids_reel = sum(pred_modele, na.rm = TRUE),
    poids_sans_dep = sum(pred_sans_dep, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    perte_totale_pct =
      100 * (poids_sans_dep - poids_reel) / poids_sans_dep
  )
grid_summary <- grid_summary_year %>%
  group_by(cell_id) %>%
  summarise(
    perte_moy = mean(perte_totale_pct, na.rm = TRUE),
    sd_perte = sd(perte_totale_pct, na.rm = TRUE),
    n_annees = n(),
    .groups = "drop"
  )
Q90 <- quantile(
  grid_summary$perte_moy,
  probs = 0.90,
  na.rm = TRUE
)

Q90
grid_plot <- left_join(
  grid_4km,
  grid_summary,
  by = "cell_id"
)

hotspots <- grid_plot %>%
  mutate(
    hotspot = perte_moy >= Q90
  )
ggplot(hotspots) +
  geom_sf(aes(fill = hotspot), color = NA) +
  scale_fill_manual(values = c("grey90", "red")) +
  theme_bw()

hot <- hotspots %>%
  filter(hotspot)

centres <- st_centroid(hot)
cellsize <- 7

voisins <- st_is_within_distance(
  centres,
  dist = 10
)

g <- graph_from_adj_list(voisins)

hot$hotspot_id <- components(g)$membership

ggplot(hot) +
  geom_sf(aes(fill = factor(hotspot_id)), color = NA) +
  theme_bw()

# Taille de chaque hotspot
taille_hotspots <- hot %>%
  st_drop_geometry() %>%
  count(hotspot_id, name = "n_cellules")

# On rattache la taille à chaque cellule
hot <- hot %>%
  left_join(taille_hotspots, by = "hotspot_id")

# On ne conserve que les hotspots d'au moins 3 cellules
hot <- hot %>%
  filter(n_cellules >= 10)
ggplot(hot) +
  geom_sf(aes(fill = factor(hotspot_id)), color = NA) +
  theme_bw()
hot_4326 <- st_transform(hot, 4326)
Ker_map_hotspots <- KER_map +
  geom_sf(
    data = hot_4326,
    aes(fill = factor(hotspot_id)),
    color = "black",
    linewidth = 0.2,
    alpha = 0.7,
    inherit.aes = FALSE
  ) +
  scale_fill_viridis_d(
    name = "Hotspot"
  ) +
  theme_bw()

Ker_map_hotspots


hotspots_poly <- hot %>%
  group_by(hotspot_id) %>%
  summarise(
    n_cellules = first(n_cellules),
    .groups = "drop"
  )
hotspots_smooth <- hotspots_poly %>%
  st_buffer(dist = 1) %>% 
  st_buffer(dist = -1)
hotspots_smooth_4326 <- st_transform(hotspots_smooth, 4326)

Ker_map_hotspots <- KER_map +
  geom_sf(
    data = hotspots_smooth_4326,
    aes(fill = factor(hotspot_id)),
    color = "black",
    linewidth = 0.5,
    alpha = 0.6,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    name = "Hotspot",
    values = c("red", "orange")
  )

Ker_map_hotspots




stations_grid <- stations_grid %>%
  filter(TAAFseason != 2026)

stations_grid <- stations_grid %>%
  mutate(
    periode = case_when(
      TAAFseason >= 2008 & TAAFseason <= 2013 ~ "2008-2013",
      TAAFseason >= 2014 & TAAFseason <= 2019 ~ "2014-2019",
      TAAFseason >= 2020 & TAAFseason <= 2025 ~ "2020-2025"
    )
  )

grid_summary_period <- stations_grid %>%
  st_drop_geometry() %>%
  group_by(periode, cell_id, TAAFseason) %>%
  summarise(
    poids_reel = sum(pred_modele, na.rm = TRUE),
    poids_sans_dep = sum(pred_sans_dep, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    perte_totale_pct =
      100 * (poids_sans_dep - poids_reel) / poids_sans_dep
  )
grid_mean_period <- grid_summary_period %>%
  group_by(periode, cell_id) %>%
  summarise(
    perte_moy = mean(perte_totale_pct, na.rm = TRUE),
    .groups = "drop"
  )

Q90_period <- grid_mean_period %>%
  group_by(periode) %>%
  summarise(
    Q90 = quantile(perte_moy, 0.90, na.rm = TRUE)
  )

hot_cells_period <- grid_mean_period %>%
  left_join(Q90_period, by = "periode") %>%
  mutate(
    hotspot = perte_moy >= Q90
  ) %>%
  filter(hotspot)

hot_sf_period <- grid_4km %>%
  right_join(
    hot_cells_period,
    by = "cell_id"
  )
find_hotspots <- function(df){
  
  centres <- st_centroid(df)
  
  voisins <- st_is_within_distance(
    centres,
    dist = 10   # une cellule de 7 km
  )
  
  edges <- do.call(
    rbind,
    lapply(seq_along(voisins), function(i){
      if(length(voisins[[i]]) == 0) return(NULL)
      cbind(i, voisins[[i]])
    })
  )
  
  if(is.null(edges)){
    df$hotspot_id <- 1
  } else {
    
    g <- graph_from_edgelist(
      edges,
      directed = FALSE
    )
    
    df$hotspot_id <- components(g)$membership
  }
  
  df
}


hotspots_period <- hot_sf_period %>%
  split(.$periode) %>%
  lapply(find_hotspots)
hotspots_period <- lapply(
  hotspots_period,
  function(x){
    
    x %>%
      add_count(
        periode,
        hotspot_id,
        name = "n_cellules"
      ) %>%
      filter(n_cellules >= 10)
    
  }
)
hotspots_all <- do.call(
  rbind,
  hotspots_period
)
hotspots_poly <- hotspots_all %>%
  group_by(periode, hotspot_id) %>%
  summarise(
    n_cellules = first(n_cellules),
    perte_moy = mean(perte_moy, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  st_buffer(1) %>%
  st_buffer(-1)
hotspots_poly_4326 <- st_transform(
  hotspots_poly,
  4326
)
Ker_map_hotspots <- KER_map +
  geom_sf(
    data = hotspots_poly_4326,
    aes(fill = factor(hotspot_id)),
    color = "black",
    linewidth = 0.4,
    alpha = 0.6,
    inherit.aes = FALSE
  ) +
  facet_wrap(~periode) +
  scale_fill_viridis_d(
    name = "Hotspot"
  ) +
  theme_bw()

Ker_map_hotspots


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
stopifnot("TAAFseason" %in% names(dat_pos))
stopifnot(is.numeric(dat_pos$TAAFseason))

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
    Month = 2,
    orque_obs = 1,
    nb_ham = 1,
    NB_MEAN_CACHALOT = 0,
    NB_MEAN_ORQUE = 0,
    Vessel = factor("Albius", levels = levels(dat_pos$Vessel)),
    HookLength = mean(dat_pos$HookLength)
  ) %>%
  relocate(cell_id)

grid_all_years$bottomT <- mapply(
  get_mean_bottomT_season,
  lon_sta = grid_all_years$X,
  lat_sta = grid_all_years$Y,
  season = grid_all_years$TAAFseason,
  MoreArgs = list(
    dossier_nc = "Data/CROZET_final/bottomT/nc_daily"
  )
)
#Vérifier colonnes requises
required_cols <- c("prof_bathy", "TAAFseason", "X", "Y")
missing_cols <- setdiff(required_cols, names(grid_all_years))
if (length(missing_cols) > 0) stop("Colonnes manquantes dans newdata : ", paste(missing_cols, collapse = ", "))


pred_std <- predict(
  m_1,
  newdata = grid_all_years,
  offset = log(grid_all_years$nb_ham * 0.246),
  return_tmb_object = TRUE
)

index_std <- get_index(pred_std, bias_correct = FALSE)



pred <- predict(
  m_2,
  newdata = grid_all_years,
  offset = log(grid_all_years$nb_ham * 0.246),
  return_tmb_object = TRUE
)

index <- get_index(pred, bias_correct = FALSE)

saveRDS(
  pred_std,
  file = file.path("data", output_plot, "pred_depred_ker.rds")
)

saveRDS(
  index_std,
  file = file.path("data", output_plot, "index_depred_ker.rds")
)

# Modèle m_2
saveRDS(
  pred,
  file = file.path("data", output_plot, "pred_ker.rds")
)

saveRDS(
  index,
  file = file.path("data", output_plot, "index_ker.rds")
)

cat("→ Prédictions et indices sauvegardés\n\n")
