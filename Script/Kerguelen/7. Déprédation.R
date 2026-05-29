
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
  ScanWeightKg ~ s(Year) + s(Month, bs = "cc"),
  offset = log(dat_pos$nb_ham * 0.246),
  data = dat_pos,
  mesh = bspde_pos,
  family = Gamma(link = "log"),
  spatial = "on",
  time= "Year",
  spatiotemporal = "iid"
)


m_2 <- sdmTMB(
  ScanWeightKg ~ s(Year) + s(Month, bs = "cc") + 
    s(NB_MEAN_CACHALOT)  ,
  offset = log(dat_pos$nb_ham * 0.246),
  data = dat_pos,
  mesh = bspde_pos,
  family = Gamma(link = "log"),
  spatial = "on",
  time= "Year",
  spatiotemporal = "iid"
)

visreg(m_2, xvar = "NB_MEAN_CACHALOT", gg = TRUE)
#visreg(m_2, xvar = "cachalot_obs", gg = TRUE)
visreg(m_2, xvar = "Month", gg = TRUE)


m_season_varying <- sdmTMB(
  ScanWeightKg ~ s(Year) + 
    s(NB_MEAN_CACHALOT)   +
    cachalot_obs,
  offset = log(dat_pos$nb_ham * 0.246),
  data = dat_pos,
  mesh = bspde_pos,
  family = Gamma(link = "log"),
  spatial = "on",
  time = "Year",
  spatiotemporal = "iid"
)





visreg(m_season_varying, xvar = "NB_MEAN_CACHALOT", gg = TRUE)
visreg(m_season_varying, xvar = "cachalot_obs", gg = TRUE)
visreg(m_season_varying, xvar = "Month", gg = TRUE)


m_season_varying_1 <- sdmTMB(
  ScanWeightKg ~ s(Year) + 
    s(Month, bs = "cc")  +
    s(NB_MEAN_CACHALOT)  ,
  offset = log(dat_pos$nb_ham * 0.246),
  data = dat_pos,
  mesh = bspde_pos,
  family = Gamma(link = "log"),
  spatial = "on",
  time = "Year",
  spatiotemporal = "iid"
)

visreg(m_season_varying_1, xvar = "NB_MEAN_CACHALOT", gg = TRUE)
visreg(m_season_varying_1, xvar = "Month", gg = TRUE)
visreg2d(m_season_varying_1, xvar = "NB_MEAN_CACHALOT", yvar = "Month", plot.type = "image")


#Comparaison des AIC :
AIC(m_1)
AIC(m_2)
AIC(m_season_varying)
AIC(m_season_varying_1)
 # -> le modèle ne converge pas





#Visualitation des effets

visreg(m_season_varying_1, xvar = "Year", gg = TRUE)
visreg(m_season_varying_1, xvar = "cachalot_obs", gg = TRUE)
visreg(m_season_varying_1, xvar = "NB_MEAN_CACHALOT", gg = TRUE)
visreg(m_season_varying_1, xvar = "Month", gg = TRUE)
visreg2d(m_season_varying_1, xvar = "NB_MEAN_CACHALOT", yvar = "Month", plot.type = "image")








#Calcul des pertes par station
# =========================
# Prédictions par station
# =========================



# 2. Jeux de données pour les scénarios
# -------------------------------------

# Situation réelle (avec déprédation observée)
dat_reel <- dat_pos

# Sans aucune déprédation
dat_sans_dep <- dat_pos %>%
  mutate(
    NB_MEAN_CACHALOT = 0,
    NB_MEAN_ORQUE = 0
  )


# 3. Offsets
# ----------
offset_reel <- log(0.246 * dat_reel$nb_ham)
offset_sans_dep <- log(0.246 * dat_sans_dep$nb_ham)
offset_sans_orque <- log(0.246 * dat_sans_orque$nb_ham)
offset_sans_cachalot <- log(0.246 * dat_sans_cachalot$nb_ham)

# 4. Prédictions du modèle
# ------------------------

pred_reel <- predict(
  m_season_varying_2,
  newdata = dat_reel,
  offset = offset_reel,  # si on veut une CPUE en sortie
  type = "response"
)

pred_sans_dep <- predict(
  m_season_varying_2,
  newdata = dat_sans_dep,
  offset = offset_sans_dep,
  type = "response"
)


# 5. Dataframe final
# ------------------

df_station_predictions <- dat_pos %>%
  mutate(
    poids_reel = ScanWeightKg,
    
    pred_modele = pred_reel$est,
    
    pred_sans_dep = pred_sans_dep$est,
    
    perte_totale = pred_sans_dep - pred_modele
  ) %>%
  select(
    StationID,
    Year,
    Month,
    ScanWeightKg,
    nb_ham,
    NB_MEAN_CACHALOT,
    X,
    Y,
    cachalot_obs,
    poids_reel,
    pred_modele,
    pred_sans_dep,
    perte_totale
  )

# Résultat
head(df_station_predictions)

saveRDS(df_station_predictions,
        "prediction_depredation_kerguelen.rds")

# Création des pourcentages
df_boxplot_pct <- df_station_predictions %>%
  mutate(
    perte_totale_pct =
      100 * (pred_sans_dep - pred_modele) / pred_sans_dep
  ) %>%
  
  select(
    Year,
    perte_totale_pct,
  ) %>%
  
  pivot_longer(
    cols = -Year,
    names_to = "type_perte",
    values_to = "pourcentage_perte"
  ) %>%
  
  mutate(
    type_perte = recode(
      type_perte,
      perte_totale_pct = "Déprédation totale",
      perte_orque_pct = "Perte orques",
      perte_cachalot_pct = "Perte cachalots"
    )
  )

# Boxplots
ggplot(df_boxplot_pct,
       aes(x = factor(Year),
           y = pourcentage_perte,
           fill = type_perte)) +
  
  geom_boxplot(
    position = position_dodge(width = 0.8),
    alpha = 0.7
  ) +
  
  labs(
    title = "Pourcentage de pertes par année",
    x = "Année",
    y = "Perte (%)",
    fill = "Type de perte"
  ) +
  
  theme_minimal()


# -> faire plutot la somme dans le carré et après le pourcentage
# -> interessant sur la résultante déprédaté
# -> peut essayer un kde
# -> série par groupe de 5 ans
# -> peut projeter en retirant effet fixe pour l'intéraction '
# -> sortir tableau et envoyer
# -> 



df_summary <- df_boxplot_pct %>%
  group_by(Year, type_perte) %>%
  summarise(
    mean_perte = mean(pourcentage_perte, na.rm = TRUE),
    q25 = quantile(pourcentage_perte, 0.25, na.rm = TRUE),
    q75 = quantile(pourcentage_perte, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(df_summary,
       aes(x = Year,
           y = mean_perte,
           color = type_perte)) +
  
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  
  geom_ribbon(
    aes(ymin = q25, ymax = q75, fill = type_perte),
    alpha = 0.15,
    color = NA
  ) +
  
  theme_minimal() +
  labs(
    x = "Année",
    y = "Perte moyenne (%)",
    color = "Scénario",
    fill = "Scénario"
  )





map = terra::rast(here(path_rast_bathy, 'ras_pecheker_interpol_tot.tif'))


map_grid_m <- terra::project(map, "EPSG:32743")
print(map_grid_m)

res <- 10000

template <- rast(
  ext = ext(map_grid_m),
  resolution = res,
  crs = crs(map_grid_m)
)

origin(template) <- origin(map_grid_m)
print(template)
cellSize(template)

map_grid_1km <- resample(map_grid_m, template, method = resampling_method)
print(map_grid_1km)


r_id <- rast(map_grid_1km)
values(r_id) <- 1:ncell(r_id)

df_station_predictions <- df_station_predictions %>%
  mutate(
    X = X*1000 ,
    Y = Y *1000
  )
df_station_predictions <- df_station_predictions %>%
  mutate(
    perte_totale_pct =
      100 * (pred_sans_dep - pred_modele) / pred_sans_dep)

pts <- vect(df_station_predictions, geom = c("X", "Y"), crs = crs(map_grid_1km))

df_station_predictions$cell_id <- extract(r_id, pts)[,2]
sum(is.na(df_station_predictions$cell_id))
loss_grid <- df_station_predictions %>%
  group_by(cell_id) %>%
  summarise(
    perte_totale_pct_mean = mean(perte_totale_pct, na.rm = TRUE)
  )
grid_df <- as.data.frame(r_id, xy = TRUE)
names(grid_df)[3] <- "cell_id"
grid_final <- grid_df %>%
  left_join(loss_grid, by = "cell_id")


ggplot(grid_final) +
  geom_raster(aes(x = x, y = y, fill = perte_totale_pct_mean)) +
  scale_fill_viridis_c(na.value = "transparent") +
  coord_equal() +
  coord_cartesian(
    xlim = c(0,1000000),
    #ylim = c(4505000, 4700000)
  )




