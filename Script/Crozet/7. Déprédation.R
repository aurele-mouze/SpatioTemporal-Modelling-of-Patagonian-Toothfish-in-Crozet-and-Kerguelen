
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
  ScanWeightKg ~ s(Year) + s(Month, bs = "cc"),
  offset = log(dat_pos$nb_ham * 0.246),
  data = dat_pos,
  mesh = bspde_pos,
  family = Gamma(link = "log"),
  spatial = "on",
  time= "Year",
  spatiotemporal = "iid"
)




m_season_varying <- sdmTMB(
  ScanWeightKg ~ s(Year) + 
    s(NB_MEAN_CACHALOT) +  s(NB_MEAN_ORQUE)  + s(Month, bs = "cc") +
    cachalot_obs + orque_obs,
  offset = log(dat_pos$nb_ham * 0.246),
  data = dat_pos,
  mesh = bspde_pos,
  family = Gamma(link = "log"),
  spatial = "on",
  time = "Year",
  spatiotemporal = "iid"
)


# m_season_varying_inter <- sdmTMB(
#   ScanWeightKg ~ s(Year) + 
#     s(Month, bs = "cc") +
#     cachalot_obs + orque_obs + s(NB_MEAN_CACHALOT, NB_MEAN_ORQUE),
#   offset = log(dat_pos$nb_ham * 0.246),
#   data = dat_pos,
#   mesh = bspde_pos,
#   family = Gamma(link = "log"),
#   spatial = "on",
#   time = "Year",
#   spatiotemporal = "iid"
# )


m_season_varying_1 <- sdmTMB(
  ScanWeightKg ~ s(Year) + 
    s(Month, bs = "cc") + s(NB_MEAN_ORQUE, Month) +
    s(NB_MEAN_CACHALOT, Month) +
    cachalot_obs + orque_obs,
  offset = log(dat_pos$nb_ham * 0.246),
  data = dat_pos,
  mesh = bspde_pos,
  family = Gamma(link = "log"),
  spatial = "on",
  time = "Year",
  spatiotemporal = "iid"
)


m_season_varying_2 <- sdmTMB(
  ScanWeightKg ~ s(Year) + 
    s(Month, bs = "cc") + s(NB_MEAN_ORQUE, Month) +
    s(NB_MEAN_CACHALOT, Month) +
    cachalot_obs + orque_obs + s(NB_MEAN_CACHALOT, NB_MEAN_ORQUE),
  offset = log(dat_pos$nb_ham * 0.246),
  data = dat_pos,
  mesh = bspde_pos,
  family = Gamma(link = "log"),
  spatial = "on",
  time = "Year",
  spatiotemporal = "iid"
)


#Comparaison des AIC :
AIC(m_2)
AIC(m_season_varying)
AIC(m_season_varying_inter)
AIC(m_season_varying_1)
AIC(m_season_varying_2) # -> le modèle ne converge pas





#Visualitation des effets

visreg(m_season_varying_1, xvar = "NB_MEAN_CACHALOT", gg = TRUE)
visreg(m_season_varying_1, xvar = "NB_MEAN_ORQUE", gg = TRUE)
visreg(m_season_varying_1, xvar = "Month", gg = TRUE)
visreg2d(m_season_varying_1, xvar = "NB_MEAN_CACHALOT", yvar = "Month", plot.type = "image")
visreg2d(m_season_varying_1, xvar = "NB_MEAN_ORQUE", yvar = "Month", plot.type = "image")
visreg2d(m_season_varying_2, xvar = "NB_MEAN_CACHALOT", yvar = "NB_MEAN_ORQUE", plot.type = "image")


visreg(m_season_varying, xvar = "NB_MEAN_CACHALOT", gg = TRUE)
visreg(m_season_varying, xvar = "NB_MEAN_ORQUE", gg = TRUE)





# Contrôle des modèles

dat_pos$resids <- residuals(m_season_varying_1)
hist_residual <- ggplot(dat_pos, aes(x = resids)) +
  geom_histogram(
    bins = 30,
    fill = "steelblue",
    color = "black",
    alpha = 0.7
  ) +
  labs(title = "Histogramme des résidus", x = "Résidus", y = "Fréquence") +
  theme_bw()

s <- simulate(m_season_varying_1, nsim = number_sim, type = "mle-mvn")
dharma_qqplot = dharma_residuals(s, m_season_varying_1, test_uniformity = TRUE)

QQplot <- ggplot(dharma_qqplot, aes(x = expected, y = observed)) +
  geom_point(alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0) +
  theme_minimal()

if (SavePlots) {
  ggsave(
    filename = file.path("data/", output_plot, "residual-histogramme-m_season_varying_1.png"),
    plot     = hist_residual,
    width    = 8,
    height   = 6,
    dpi      = 300
  )
}

if (SavePlots) {
  ggsave(
    filename = file.path("data/", output_plot, "qqplot_residual-m_season_varying_1.png"),
    plot     = QQplot,
    width    = 8,
    height   = 6,
    dpi      = 300
  )
}






# Fabrication d'une grille d'extrapolation ------------------------------------------------------
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




# Prédictions -----------------------------------
#Indice d'abondance


# Vérifications rapides
stopifnot("Year" %in% names(dat_pos))
stopifnot(is.numeric(dat_pos$Year))

#Préparer grid avec un identifiant de cellule
grid_with_id <- grid %>%
  mutate(cell_id = row_number())

grid_with_id <- dplyr::as_tibble(grid_with_id)

#Réplication temporelle : une ligne par cellule × année
years <- sort(unique(dat_pos$Year))
month <- sort(unique(dat_pos$Month))
grid_all_years <- grid_with_id %>%
  tidyr::uncount(length(years)*length(month)) %>%
  mutate(Year = rep(years, times = nrow(grid_with_id)*length(month)),
         Month = rep(month, times = nrow(grid_with_id)*length(years))) %>%
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


grid_index_0 <- grid_with_id %>%
  crossing(Year = years) %>%
  mutate(
    cachalot_obs = 1,
    orque_obs = 1,
    nb_ham = mean(dat_pos$nb_ham, na.rm = TRUE),
    NB_MEAN_CACHALOT = 0,
    NB_MEAN_ORQUE = 0,
    scenario = "Modèle de déprédation"
  )



index_0 <- get_index_split(
  m_season_varying_1,
  newdata = grid_index_0,
  area = 1,
  nsplit = 3
)

grid_index <- grid_with_id %>%
  crossing(Year = years, Month = 1:12) %>%
  mutate(
    nb_ham = mean(dat_pos$nb_ham, na.rm = TRUE),
    scenario = "Modèle sans déprédation"
  )

index_ <- get_index_split(
  m_season_varying_1,
  newdata = grid_index,
  area = 1,
  nsplit = 3
)


index_0_year <- index_0 %>%
  group_by(Year) %>%
  summarise(
    est = weighted.mean(est, w = nb_ham),
    lwr = weighted.mean(lwr, w = nb_ham),
    upr = weighted.mean(upr, w = nb_ham)
  ) %>%
  mutate(scenario = "Modèle avec déprédation")

index_year <- index_ %>%
  group_by(Year) %>%
  summarise(
    est = weighted.mean(est, w = nb_ham),
    lwr = weighted.mean(lwr, w = nb_ham),
    upr = weighted.mean(upr, w = nb_ham)
  ) %>%
  mutate(scenario = "Modèle sans déprédation")


index_comp <- bind_rows(index_0_year, index_year)


Abundance_Index <- ggplot(index_comp, aes(Year, est, color = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2) +
  labs(
    x = "Year",
    y = "Biomasse",
    color = "Scenario",
    fill = "Scenario"
  ) +
  theme_bw()

if (SavePlots) {
  ggsave(
    filename = file.path("data/", output_plot, "Indice_abondance_m_spatially_varying.png"),
    plot     = Abundance_Index,
    width    = 8,
    height   = 6,
    dpi      = 300
  )
}



if (SavePlots) {
  ggsave(
    filename = file.path("data/", output_plot, "Indice_abondance_m_spatially_varying.png"),
    plot     = Abundance_Index,
    width    = 8,
    height   = 6,
    dpi      = 300
  )
}


diff_index <- index_comp %>%
  # garder la référence
  filter(scenario == "Modèle avec déprédation") %>%
  select(Year, est_ref = est) %>%
  
  # joindre avec tous les scénarios
  right_join(index_comp, by = "Year") %>%
  
  # calcul de la différence
  mutate(diff = est - est_ref)




diff_index <- index_comp %>%
  filter(scenario == "Modèle avec déprédation") %>%
  select(Year, est_ref = est) %>%
  right_join(index_comp, by = "Year") %>%
  mutate(
    perc_diff = (est - est_ref) / est_ref * 100
  )


ggplot(diff_index, aes(Year, perc_diff, color = scenario)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = 1.2) +
  labs(
    x = "Year",
    y = "Variation (%) par rapport au scénario sans déprédation",
    color = "Scenario"
  ) +
  theme_bw()





#Calcul du taux de déprédation avec les stations

pred_data <- dat_pos

pred_no_dep <- pred_data %>%
  mutate(
    NB_MEAN_CACHALOT = if_else(cachalot_obs == 1, 0, NB_MEAN_CACHALOT),
    NB_MEAN_ORQUE = if_else(orque_obs == 1, 0, NB_MEAN_ORQUE),
    scenario = "Sans déprédation"
  )

pred_orque <- pred_data %>%
  mutate(
    NB_MEAN_CACHALOT = if_else(cachalot_obs == 1, 0, NB_MEAN_CACHALOT),
    scenario = "Déprédation Orque"
  )

pred_Cachalot <- pred_data %>%
  mutate(
    NB_MEAN_ORQUE = if_else(orque_obs == 1, 0, NB_MEAN_ORQUE),
    scenario = "Déprédation Cachalot"
  )


pred_all <- bind_rows(
  pred_no_dep,
  pred_orque,
  pred_Cachalot
)


pred_all <- pred_all %>%
  mutate(
    offset_val = log(nb_ham * 0.246)
  )

pred_all$poids_attendu <- predict(
  m_season_varying_1,
  newdata = pred_all,
  offset = pred_all$offset_val,
  type = "response"
)$est


comparaison <- pred_all %>%
  select(StationID, scenario, poids_attendu, ScanWeightKg, Year) %>%
  pivot_wider(
    names_from = scenario,
    values_from = poids_attendu
  )

comparaison <- comparaison %>%
  mutate(
    perte_orque =
      `Sans déprédation` - `Déprédation Orque`,
    
    perte_cachalot =
      `Sans déprédation` - `Déprédation Cachalot`,
    
    perte_vs_reel =
      `Sans déprédation` - ScanWeightKg
  )


comparaison <- comparaison %>%
  mutate(
    pct_perte_orque =
      100 * (perte_orque / `Sans déprédation`),
    
    pct_perte_cachalot =
      100 * (perte_cachalot / `Sans déprédation`),
    
    pct_perte_vs_reel =
      100 * ((`Sans déprédation` - ScanWeightKg) / `Sans déprédation`)
  )



plot_data <- comparaison %>%
  select(
    pct_perte_orque,
    pct_perte_cachalot,
    pct_perte_vs_reel
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "scenario",
    values_to = "pct_perte"
  )


ggplot(plot_data, aes(x = scenario, y = pct_perte)) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    x = "Scénario",
    y = "Perte (%)",
    title = "Distribution des pertes de capture selon les scénarios"
  )


plot_data_year <- comparaison %>%
  select(
    Year,
    pct_perte_orque,
    pct_perte_cachalot,
    pct_perte_vs_reel
  ) %>%
  pivot_longer(
    cols = starts_with("pct_perte"),
    names_to = "scenario",
    values_to = "pct_perte"
  ) %>%
  mutate(
    scenario = recode(
      scenario,
      pct_perte_orque = "Orque",
      pct_perte_cachalot = "Cachalot",
      pct_perte_vs_reel = "Modèle vs réel"
    )
  )

ggplot(plot_data_year, aes(x = factor(Year), y = pct_perte)) +
  geom_boxplot() +
  facet_wrap(~ scenario) +
  theme_minimal() +
  labs(
    x = "Année",
    y = "Perte (%)",
    title = "Évolution annuelle des pertes par scénario"
  )




#Prédiction sur la data pour calculer les pertes

# Vérifications
stopifnot("Year" %in% names(dat_pos))
stopifnot(is.numeric(dat_pos$Year))

# =========================
# 1. Données AVEC déprédation
# =========================

dat_pred_dep <- dat_pos %>%
  mutate(
    scenario = "Avec déprédation"
  )

pred_dep <- predict(
  m_season_varying_1,
  newdata = dat_pred_dep,
  return_tmb_object = TRUE
)

index_dep <- get_index(
  pred_dep,
  area = 0.246 * dat_pred_dep$nb_ham
)

# =========================
# 2. Données SANS déprédation
# =========================

dat_pred_nodep <- dat_pos %>%
  mutate(
    NB_MEAN_CACHALOT = 0,
    NB_MEAN_ORQUE = 0,
    scenario = "Sans déprédation"
  )

pred_nodep <- predict(
  m_season_varying_1,
  newdata = dat_pred_nodep,
  return_tmb_object = TRUE
)

index_nodep <- get_index(
  pred_nodep,
  area = 0.246 * dat_pred_nodep$nb_ham
)

# =========================
# 3. Résumé annuel
# =========================

index_dep_year <- index_dep %>%
  group_by(Year) %>%
  summarise(
    est = weighted.mean(est, w = nb_ham, na.rm = TRUE),
    lwr = weighted.mean(lwr, w = nb_ham, na.rm = TRUE),
    upr = weighted.mean(upr, w = nb_ham, na.rm = TRUE)
  ) %>%
  mutate(scenario = "Avec déprédation")

index_nodep_year <- index_nodep %>%
  group_by(Year) %>%
  summarise(
    est = weighted.mean(est, w = nb_ham, na.rm = TRUE),
    lwr = weighted.mean(lwr, w = nb_ham, na.rm = TRUE),
    upr = weighted.mean(upr, w = nb_ham, na.rm = TRUE)
  ) %>%
  mutate(scenario = "Sans déprédation")

# =========================
# 4. Combinaison
# =========================

index_comp <- bind_rows(
  index_dep_year,
  index_nodep_year
)

# =========================
# 5. Figure
# =========================

Abundance_Index <- ggplot(
  index_comp,
  aes(Year, est, color = scenario, fill = scenario)
) +
  geom_ribbon(
    aes(ymin = lwr, ymax = upr),
    alpha = 0.2,
    color = NA
  ) +
  geom_line(linewidth = 1.2) +
  labs(
    x = "Year",
    y = "Biomasse prédite",
    color = "Scénario",
    fill = "Scénario"
  ) +
  theme_bw()

Abundance_Index


index_loss_year <- index_dep_year %>%
  left_join(
    index_nodep_year,
    by = "Year",
    suffix = c("_dep", "_nodep")
  ) %>%
  mutate(
    perte_pct = 100 * (est_nodep - est_dep) / est_nodep
  )

ggplot(index_loss_year, aes(x = Year, y = perte_pct)) +
  geom_line(linewidth = 1.2) +
  geom_point() +
  theme_bw() +
  labs(
    x = "Année",
    y = "% de biomasse déprédatée",
    title = "Impact annuel de la déprédation"
  )












#Calcul des pertes par station

dat_pred_dep <- dat_pos %>%
  mutate(
    NB_MEAN_CACHALOT = NB_MEAN_CACHALOT,
    NB_MEAN_ORQUE = NB_MEAN_ORQUE
  )

dat_pred_nodep <- dat_pos %>%
  mutate(
    NB_MEAN_CACHALOT = 0,
    NB_MEAN_ORQUE = 0
  )

offset_dep <- log(0.246 * dat_pred_dep$nb_ham)
offset_nodep <- log(0.246 * dat_pred_nodep$nb_ham)

pred_dep <- predict(
  m_season_varying_1,
  newdata = dat_pred_dep,
  offset = offset_dep
)

pred_nodep <- predict(
  m_season_varying_1,
  newdata = dat_pred_nodep,
  offset = offset_nodep
)


df_dep <- dat_pred_dep %>%
  select(StationID, Year) %>%
  mutate(est_dep = pred_dep$est)

df_nodep <- dat_pred_nodep %>%
  select(StationID, Year) %>%
  mutate(est_nodep = pred_nodep$est)

df_compare <- dat_pos %>%
  select(StationID, Year) %>%
  left_join(df_dep, by = c("StationID", "Year")) %>%
  left_join(df_nodep, by = c("StationID", "Year")) %>%
  mutate(
    perte_abs = est_nodep - est_dep,
    perte_pct = 100 * (est_nodep - est_dep) / est_nodep
  )


perte_station <- df_compare %>%
  distinct(StationID, Year, .keep_all = TRUE)

ggplot(df_compare, aes(x = Year, y = perte_pct, group = Year)) +
  geom_boxplot() +
  theme_bw() +
  labs(
    x = "Année",
    y = "% de perte",
    title = "Distribution annuelle de l’impact de la déprédation"
  )