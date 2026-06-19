
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
  filter(,
    !is.na(ScanWeightKg)
  )


dat_pos <- dat_pos %>%
  filter(,
         !is.na(NB_MEAN_ORQUE)
  )

dat_pos <- dat_pos %>%
  filter(,
         !is.na(tps_trempe_hm_num)
  )

dat_pos <- dat_pos %>% filter(tps_trempe_hm_num < 100)

dat_pos <- dat_pos %>%
  mutate(
    NB_MEAN_CACHALOT = ifelse(is.na(NB_MEAN_CACHALOT), 0, NB_MEAN_CACHALOT)
  )

dat_pos$Vessel <- droplevels(dat_pos$Vessel)


dat_pos$HookLength <- as.numeric(sub("_.*", "", dat_pos$HookSize))


dat_pos <- dat_pos %>%
  filter(HookLength != 7)


mean_value <- mean(dat_pos$HookLength, na.rm = TRUE)


dat_pos$HookLength[is.na(dat_pos$HookLength)] <- mean_value



# dat_pos$bottomT[is.nan(dat_pos$bottomT)] <- NA
# 
# dat_pos <- dat_pos %>%
#   filter(!is.na(bottomT))

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



# 
# 
# ggplot(dat_pos, aes(x = Vessel, y = ScanWeightKg)) +
#   geom_boxplot() +
#   theme_minimal()
# 
# ggplot(dat_pos, aes(x = Captain, y = ScanWeightKg)) +
#   geom_boxplot() +
#   theme_minimal()
# 
# ggplot(dat_pos, aes(x = HookSize, y = ScanWeightKg)) +
#   geom_boxplot() +
#   theme_minimal()
# 
# dat_count <- dat_pos %>%
#   group_by(Year, HookSize) %>%
#   summarise(n_stations = n(), .groups = "drop")
# 
# ggplot(dat_count,
#        aes(x = Year, y = n_stations, color = factor(HookSize))) +
#   geom_line() +
#   geom_point() +
#   theme_minimal() +
#   labs(
#     x = "Année",
#     y = "Nombre de stations",
#     color = "HookSize"
#   )
# 
# 
# ggplot(dat_pos, aes(x = prof_bathy, y = ScanWeightKg)) +
#   geom_point(alpha = 0.4) +
#   geom_smooth(method = "loess") +
#   theme_minimal()
# 
# ggplot(dat_pos, aes(x = slope, y = ScanWeightKg)) +
#   geom_point(alpha = 0.4) +
#   geom_smooth(method = "loess") +
#   theme_minimal()
# 
# ggplot(dat_pos, aes(x = tps_trempe_hm_num, y = ScanWeightKg)) +
#   geom_point(alpha = 0.4) +
#   geom_smooth(method = "loess") +
#   theme_minimal()
# 



m_1 <- sdmTMB(
  ScanWeightKg ~ s(TAAFseason) + s(Month, bs = "cc") + Vessel + s(prof_bathy) + s(slope) + s(HookLength),
  offset = log(dat_pos$nb_ham * 0.246),
  data = dat_pos,
  mesh = bspde_pos,
  family = Gamma(link = "log"),
  spatial = "on",
  time= "TAAFseason",
  spatiotemporal = "iid"
)


m_2 <- sdmTMB(
  ScanWeightKg ~ s(TAAFseason) + s(Month, bs = "cc") + Vessel + s(prof_bathy) + s(slope) + s(HookLength) +
    s(NB_MEAN_CACHALOT)  ,
  offset = log(dat_pos$nb_ham * 0.246),
  data = dat_pos,
  mesh = bspde_pos,
  family = Gamma(link = "log"),
  spatial = "on",
  time= "TAAFseason",
  spatiotemporal = "iid"
)
# 
# visreg(m_2, xvar = "NB_MEAN_CACHALOT", gg = TRUE)
# #visreg(m_2, xvar = "cachalot_obs", gg = TRUE)
# visreg(m_2, xvar = "Month", gg = TRUE)
# 
# 
# m_season_varying <- sdmTMB(
#   ScanWeightKg ~ s(Year) + 
#     s(NB_MEAN_CACHALOT)   +
#     cachalot_obs,
#   offset = log(dat_pos$nb_ham * 0.246),
#   data = dat_pos,
#   mesh = bspde_pos,
#   family = Gamma(link = "log"),
#   spatial = "on",
#   time = "Year",
#   spatiotemporal = "iid"
# )
# 
# 
# 
# 
# 
# visreg(m_season_varying, xvar = "NB_MEAN_CACHALOT", gg = TRUE)
# visreg(m_season_varying, xvar = "cachalot_obs", gg = TRUE)
# visreg(m_season_varying, xvar = "Month", gg = TRUE)
# 
# 
# m_season_varying_1 <- sdmTMB(
#   ScanWeightKg ~ s(TAAFseason) + 
#     s(Month, bs = "cc")  +
#     s(NB_MEAN_CACHALOT) ,
#   offset = log(dat_pos$nb_ham * 0.246),
#   data = dat_pos,
#   mesh = bspde_pos,
#   family = Gamma(link = "log"),
#   spatial = "on",
#   time = "TAAFseason",
#   spatiotemporal = "iid"
# )
# 
# 
# m_season_varying_2 <- sdmTMB(
#   ScanWeightKg ~ s(TAAFseason) + 
#     s(Month, bs = "cc")  +
#     s(NB_MEAN_CACHALOT)  +
#     + Vessel,
#   offset = log(dat_pos$nb_ham * 0.246),
#   data = dat_pos,
#   mesh = bspde_pos,
#   family = Gamma(link = "log"),
#   spatial = "on",
#   time = "TAAFseason",
#   spatiotemporal = "iid"
# )
# 
# 
# 
# m_season_varying_3 <- sdmTMB(
#   ScanWeightKg ~ s(TAAFseason) + 
#     s(Month, bs = "cc")  +
#     s(NB_MEAN_CACHALOT)  +
#     + Vessel + s(tps_trempe_hm_num),
#   offset = log(dat_pos$nb_ham * 0.246),
#   data = dat_pos,
#   mesh = bspde_pos,
#   family = Gamma(link = "log"),
#   spatial = "on",
#   time = "TAAFseason",
#   spatiotemporal = "iid"
# )
# 
# 
# m_season_varying_4 <- sdmTMB(
#   ScanWeightKg ~ s(TAAFseason) + 
#     s(Month, bs = "cc")  +
#     s(NB_MEAN_CACHALOT)  +
#     + Vessel + s(prof_bathy) + s(tps_trempe_hm_num),
#   offset = log(dat_pos$nb_ham * 0.246),
#   data = dat_pos,
#   mesh = bspde_pos,
#   family = Gamma(link = "log"),
#   spatial = "on",
#   time = "TAAFseason",
#   spatiotemporal = "iid"
# )
# 
# AIC(m_season_varying_1)
# AIC(m_season_varying_2) 
# AIC(m_season_varying_3)
# AIC(m_season_varying_4)
# 
# 
# visreg(m_season_varying_4, xvar = "prof_bathy", gg = TRUE)
# visreg(m_season_varying_4, xvar = "Vessel", gg = TRUE)
# visreg(m_season_varying_4, xvar = "NB_MEAN_CACHALOT", gg = TRUE)
# visreg(m_season_varying_4, xvar = "Month", gg = TRUE)
# visreg(m_season_varying_4, xvar = "tps_trempe_hm_num", gg = TRUE)
# 
# 
# 
# V <- m_season_varying_4$sd_report$cov.fixed
# Z <- m_season_varying_4$smoothers$Zs[[1]]
# Z <- Z[, 1:8]  # important !
# 
# beta <- m_season_varying_4$parlist$b_smooth[,1]
# beta_taaf <- beta[1:8]
# 
# fit <- as.numeric(Z %*% beta_taaf)
# V_beta <- V[1:8, 1:8]
# 
# ZV <- Z %*% V_beta
# 
# var_fit <- rowSums(ZV * Z)
# se_fit <- sqrt(var_fit)
# upper <- fit + 1.96 * se_fit
# lower <- fit - 1.96 * se_fit
# eff <- data.frame(
#   TAAFseason = dat_pos$TAAFseason,
#   fit = fit,
#   lower = lower,
#   upper = upper
# ) %>%
#   group_by(TAAFseason) %>%
#   summarise(
#     fit = mean(fit),
#     lower = mean(lower),
#     upper = mean(upper),
#     .groups = "drop"
#   )
# 
# ggplot(eff, aes(TAAFseason, fit)) +
#   geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
#   geom_line(linewidth = 1) +
#   theme_bw() +
#   labs(
#     x = "TAAFseason",
#     y = "Effet s(TAAFseason)",
#     title = "Smooth temporel avec IC 95%"
#   )
# 
# #Comparaison des AIC :
# AIC(m_1)
# AIC(m_2)
# AIC(m_season_varying)
# AIC(m_season_varying_1)
#  # -> le modèle ne converge pas
# 
# 
# beta <- m_season_varying_1$parlist$b_smooth[,1]
# beta_taaf <- beta[1:8]
# 
# Z_taaf <- m_season_varying_1$smoothers$Zs[[1]]
# attr(Z_taaf, "s.label")
# # "s(TAAFseason)"
# 
# 
# f_taaf <- Z_taaf %*% beta_taaf
# eff_taaf <- data.frame(
#   TAAFseason = dat_pos$TAAFseason,
#   fit = as.numeric(f_taaf)
# )
# eff_taaf <- eff_taaf %>%
#   group_by(TAAFseason) %>%
#   summarise(fit = mean(fit), .groups = "drop")
# library(ggplot2)
# 
# ggplot(eff_taaf, aes(TAAFseason, fit)) +
#   geom_line(linewidth = 1) +
#   geom_point() +
#   theme_bw() +
#   labs(
#     x = "TAAFseason",
#     y = "Effet partiel s(TAAFseason)"
#   )
# 
# 
# V <- m_season_varying_1$sd_report$cov.fixed
# Z <- m_season_varying_1$smoothers$Zs[[1]]
# Z <- Z[, 1:8]  # important !
# beta <- m_season_varying_1$parlist$b_smooth[,1]
# beta_taaf <- beta[1:8]
# 
# fit <- as.numeric(Z %*% beta_taaf)
# var_fit <- diag(Z %*% V[1:8, 1:8] %*% t(Z))
# se_fit <- sqrt(var_fit)
# upper <- fit + 1.96 * se_fit
# lower <- fit - 1.96 * se_fit
# eff <- data.frame(
#   TAAFseason = dat_pos$TAAFseason,
#   fit = fit,
#   lower = lower,
#   upper = upper
# ) %>%
#   group_by(TAAFseason) %>%
#   summarise(
#     fit = mean(fit),
#     lower = mean(lower),
#     upper = mean(upper),
#     .groups = "drop"
#   )
# 
# ggplot(eff, aes(TAAFseason, fit)) +
#   geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
#   geom_line(linewidth = 1) +
#   geom_point() +
#   theme_bw() +
#   labs(
#     x = "TAAFseason",
#     y = "Effet s(TAAFseason)",
#     title = "Smooth temporel avec IC 95%"
#   )
# 
# 
# 
# #Visualitation des effets
# 
# visreg(m_season_varying_1, xvar = "Year", gg = TRUE)
# visreg(m_season_varying_1, xvar = "cachalot_obs", gg = TRUE)
# visreg(m_season_varying_1, xvar = "NB_MEAN_CACHALOT", gg = TRUE)
# visreg(m_season_varying_1, xvar = "Month", gg = TRUE)
# visreg2d(m_season_varying_1, xvar = "NB_MEAN_CACHALOT", yvar = "Month", plot.type = "image")
# 
# 
# 
# 
# 
# 
# 
# 
# #Calcul des pertes par station
# # =========================
# # Prédictions par station
# # =========================
# 
# 
# 
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

# 4. Prédictions du modèle
# ------------------------

pred_reel <- predict(
  m_2,
  newdata = dat_reel,
  offset = offset_reel,  # si on veut une CPUE en sortie
  type = "response"
)

pred_sans_dep <- predict(
  m_2,
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
    TAAFseason,
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
    TAAFseason,
    perte_totale_pct
  )


# Boxplots
ggplot(df_boxplot_pct,
       aes(x = factor(TAAFseason),
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

# 
# # -> faire plutot la somme dans le carré et après le pourcentage
# # -> interessant sur la résultante déprédaté
# # -> peut essayer un kde
# # -> série par groupe de 5 ans
# # -> peut projeter en retirant effet fixe pour l'intéraction '
# # -> sortir tableau et envoyer
# # -> 
# 
# 
# 
df_resume <- df_station_predictions %>%
  mutate(
    perte_totale_pct = 100 * (pred_sans_dep - pred_modele) / pred_sans_dep
  ) %>%
  filter(TAAFseason >= 2010) %>%
  group_by(TAAFseason) %>%
  summarise(
    moyenne = mean(perte_totale_pct, na.rm = TRUE),
    sd = sd(perte_totale_pct, na.rm = TRUE),
    n = n(),
    se = sd / sqrt(n),
    ic95_inf = moyenne - 1.96 * se,
    ic95_sup = moyenne + 1.96 * se,
    .groups = "drop"
  )


ggplot(df_resume,
       aes(x = TAAFseason, group = 1)) +
  geom_line(aes(y = moyenne), linewidth = 1) +
  geom_line(aes(y = ic95_inf), linetype = "dashed") +
  geom_line(aes(y = ic95_sup), linetype = "dashed") +
  labs(
    x = "Saison",
    y = "Perte totale (%)",
    title = "Pourcentage moyen de perte totale avec IC à 95 %"
  ) +
  theme_classic()

ggplot(df_resume,
       aes(x = TAAFseason, group = 1)) +
  geom_ribbon(
    aes(ymin = ic95_inf, ymax = ic95_sup),
    alpha = 0.2
  ) +
  geom_line(aes(y = ic95_inf), linetype = "dashed") +
  geom_line(aes(y = ic95_sup), linetype = "dashed") +
  geom_line(aes(y = moyenne), linewidth = 1) +
  labs(
    x = "Year",
    y = "Total loss (%)",
    title = "Average percentage of total loss with a 95% confidence interval"
  ) +
  theme_classic()
# 
# 
# 
# map = terra::rast(here(path_rast_bathy, 'ras_pecheker_interpol_tot.tif'))
# 
# 
# map_grid_m <- terra::project(map, "EPSG:32743")
# print(map_grid_m)
# 
# res <- 10000
# 
# template <- rast(
#   ext = ext(map_grid_m),
#   resolution = res,
#   crs = crs(map_grid_m)
# )
# 
# origin(template) <- origin(map_grid_m)
# print(template)
# cellSize(template)
# 
# map_grid_1km <- resample(map_grid_m, template, method = resampling_method)
# print(map_grid_1km)
# 
# 
# r_id <- rast(map_grid_1km)
# values(r_id) <- 1:ncell(r_id)
# 
# df_station_predictions <- df_station_predictions %>%
#   mutate(
#     X = X*1000 ,
#     Y = Y *1000
#   )
# df_station_predictions <- df_station_predictions %>%
#   mutate(
#     perte_totale_pct =
#       100 * (pred_sans_dep - pred_modele) / pred_sans_dep)
# 
# pts <- vect(df_station_predictions, geom = c("X", "Y"), crs = crs(map_grid_1km))
# 
# df_station_predictions$cell_id <- extract(r_id, pts)[,2]
# sum(is.na(df_station_predictions$cell_id))
# loss_grid <- df_station_predictions %>%
#   group_by(cell_id) %>%
#   summarise(
#     perte_totale_pct_mean = mean(perte_totale_pct, na.rm = TRUE)
#   )
# grid_df <- as.data.frame(r_id, xy = TRUE)
# names(grid_df)[3] <- "cell_id"
# grid_final <- grid_df %>%
#   left_join(loss_grid, by = "cell_id")
# 
# 
# ggplot(grid_final) +
#   geom_raster(aes(x = x, y = y, fill = perte_totale_pct_mean)) +
#   scale_fill_viridis_c(na.value = "transparent") +
#   coord_equal() +
#   coord_cartesian(
#     xlim = c(0,1000000),
#     #ylim = c(4505000, 4700000)
#   )
# 
# 
# 
# 
# 
# #Calcul des CPUE par périodes, comparaison pour observer la baisse
# 
# nb_ham_std <- 1
# offset_std <- log(nb_ham_std * 0.246)
# 
# pred_std <- predict(
#   m_season_varying_1,
#   newdata = dat_pos,
#   offset = rep(offset_std, nrow(dat_pos)),
#   type = "response"
# )
# 
# dat_pos$CPUE_std <- pred_std$est * 1000
# 
# 
# dat_pos <- dat_pos %>%
#   mutate(
#     periode = case_when(
#       Year %in% c(2021, 2022) ~ "2021-2022",
#       Year %in% c(2022, 2023) ~ "2022-2023",
#       Year %in% c(2023, 2024) ~ "2023-2024",
#       TRUE ~ NA_character_
#     )
#   )
# 
# 
# cpue_period_2 <- dat_pos %>%
#   filter(!is.na(periode)) %>%
#   group_by(periode) %>%
#   summarise(
#     n = n(),
#     mean_cpue = mean(CPUE_std, na.rm = TRUE),
#     sd_cpue = sd(CPUE_std, na.rm = TRUE),
#     se = sd_cpue / sqrt(n),
#     IC_inf = mean_cpue - qt(0.975, df = n - 1) * se,
#     IC_sup = mean_cpue + qt(0.975, df = n - 1) * se,
#     .groups = "drop"
#   )
# 
# 
# ggplot(cpue_period_2, aes(x = periode, y = mean_cpue)) +
#   geom_point(size = 3) +
#   geom_errorbar(aes(ymin = IC_inf, ymax = IC_sup), width = 0.15) +
#   theme_classic() +
#   labs(x = "", y = "CPUE standardisée (moyenne)")
# 
# 
# ref <- cpue_period_2$mean_cpue[cpue_period_2$periode == "2021-2022"]
# 
# cpue_period_2 <- cpue_period_2 %>%
#   mutate(
#     index_rel = 100 * mean_cpue / ref,
#     pct_change = 100 * (mean_cpue / ref - 1)
#   )
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #6. Fabrication d'une grille d'extrapolation ------------------------------------------------------
# # Import et préparation de la grille d’extrapolation spatiale
# Extrapol_grid <- readRDS(here("data/outputs/df/grid_extrapol.rds"))
# 
# # Mise à l’échelle de la profondeur et de la pente sur la grille
# Mean_SD <- dat_pos %>% ungroup() %>% summarise(Mean = mean(prof_bathy), SD = sd(prof_bathy))
# Extrapol_grid = Extrapol_grid %>% rename(depth = profondeur.pred)
# Extrapol_grid <- Extrapol_grid %>% mutate(depth_scaled = (depth - Mean_SD$Mean) / Mean_SD$SD)
# 
# Mean_SD <- dat_pos %>% ungroup() %>% summarise(Mean = mean(slope), SD = sd(slope))
# Extrapol_grid <- Extrapol_grid %>% mutate(slope_scaled = (slope - Mean_SD$Mean) / Mean_SD$SD)
# 
# # Nettoyage et sélection des colonnes nécessaires pour la prédiction
# grid <- Extrapol_grid %>%
#   mutate(
#     nb_ham = nb_ham,
#     tps_h = mean(dat$tps_h),
#     x = x / 1000 ,
#     y = y / 1000 
#   ) %>%
#   filter(!is.na(slope_scaled)) %>%
#   filter(!is.na(depth_scaled)) %>%
#   filter(depth > min_depth) %>%
#   rename(X = x, Y = y) %>%
#   dplyr::select(depth_scaled, slope_scaled, slope, depth, nb_ham, tps_h, X, Y) %>%
#   rename(prof_bathy = depth)
# 
# 
# #7. Prédictions -----------------------------------
# 
# # Vérifications rapides
# stopifnot("TAAFseason" %in% names(dat_pos))
# stopifnot(is.numeric(dat_pos$TAAFseason))
# 
# #Préparer grid avec un identifiant de cellule
# grid_with_id <- grid %>%
#   mutate(cell_id = row_number())
# 
# grid_with_id <- dplyr::as_tibble(grid_with_id)
# 
# #Réplication temporelle : une ligne par cellule × année
# years <- sort(unique(dat_pos$TAAFseason))
# grid_all_years <- grid_with_id %>%
#   tidyr::uncount(length(years)) %>%
#   mutate(TAAFseason = rep(years, times = nrow(grid_with_id))) %>%
#   mutate(cachalot_obs = 1,
#          Month = 10,
#          nb_ham = 1,
#          NB_MEAN_CACHALOT = 0) %>%
#   relocate(cell_id)
# 
# #Vérifier colonnes requises
# required_cols <- c("prof_bathy", "TAAFseason", "X", "Y")
# missing_cols <- setdiff(required_cols, names(grid_all_years))
# if (length(missing_cols) > 0) stop("Colonnes manquantes dans newdata : ", paste(missing_cols, collapse = ", "))
# 
# 
# 
# # Prédiction (utilise chunking si nécessaire)
# 
# pred_std <- predict(
#   m_season_varying_1,
#   newdata = grid_all_years,
#   offset = log(grid_all_years$nb_ham * 0.246),
#   return_tmb_object = TRUE
# )
# 
# index <- get_index(pred_std, bias_correct = FALSE)
# # index <- get_index_split(
# #   m_season_varying_1,
# #   newdata = grid_all_years,
# #   bias_correct = TRUE,
# #   nsplit = 5
# # )
# 
# index %>%
#   filter(TAAFseason >= 2010, TAAFseason <= 2025) %>%
#   ggplot(aes(x = TAAFseason, y = est)) +
#   geom_line(linewidth = 1) +
#   geom_point(size = 2) +
#   geom_errorbar(
#     aes(ymin = lwr, ymax = upr),
#     width = 0.2
#   ) +
#   scale_x_continuous(breaks = 2010:2025) +
#   labs(
#     x = "Saison",
#     y = "Indice standardisé",
#     title = "Évolution de l'indice (2010–2025)"
#   ) +
#   theme_bw() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )
# 
# 
# ggplot(
#   filter(index, TAAFseason >= 2010, TAAFseason <= 2025),
#   aes(TAAFseason, est)
# ) +
#   geom_ribbon(
#     aes(ymin = lwr, ymax = upr),
#     alpha = 0.2
#   ) +
#   geom_line(linewidth = 1) +
#   geom_point(size = 2) +
#   theme_bw() +
#   labs(
#     x = "Saison de pêche",
#     y = "Indice standardisé de la CPUE"
#   )
# 
# 
# 
# 
# 
# 
# 
# 
# cpue_summary <- dat_pos %>%
#   filter(TAAFseason >= 2010, TAAFseason <= 2025) %>%
#   group_by(TAAFseason) %>%
#   summarise(
#     n = sum(!is.na(CPUE)),
#     mean_CPUE = mean(CPUE, na.rm = TRUE),
#     sd_CPUE = sd(CPUE, na.rm = TRUE),
#     se_CPUE = sd_CPUE / sqrt(n),
#     lwr = mean_CPUE - 1.96 * se_CPUE,
#     upr = mean_CPUE + 1.96 * se_CPUE,
#     .groups = "drop"
#   )
# 
# ggplot(cpue_summary, aes(TAAFseason, mean_CPUE)) +
#   geom_ribbon(
#     aes(ymin = lwr, ymax = upr),
#     alpha = 0.2
#   ) +
#   geom_line(linewidth = 1) +
#   geom_point(size = 2) +
#   theme_bw() +
#   labs(
#     x = "Saison de pêche",
#     y = "CPUE moyenne en g/hameçon"
#   )







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
years <- years[years >= 2010 & years <= 2025]
grid_all_years <- grid_with_id %>%
  tidyr::crossing(TAAFseason = years) %>%
  mutate(
    Month = 2,
    nb_ham = 1,
    NB_MEAN_CACHALOT = 0,
    Vessel = factor("Albius", levels = levels(dat_pos$Vessel)),
    HookLength = mean(dat_pos$HookLength)
  ) %>%
  relocate(cell_id)

#Vérifier colonnes requises
required_cols <- c("prof_bathy", "TAAFseason", "X", "Y")
missing_cols <- setdiff(required_cols, names(grid_all_years))
if (length(missing_cols) > 0) stop("Colonnes manquantes dans newdata : ", paste(missing_cols, collapse = ", "))


pred_std <- predict(
  m_2,
  newdata = grid_all_years,
  offset = log(grid_all_years$nb_ham * 0.246),
  return_tmb_object = TRUE
)

index_std <- get_index(pred_std, bias_correct = FALSE)




pred <- predict(
  m_1,
  newdata = grid_all_years,
  offset = log(grid_all_years$nb_ham * 0.246),
  return_tmb_object = TRUE
)

index <- get_index(pred, bias_correct = FALSE)

index_std$model <- "CPUE Depredation Model"
index$model <- "CPUE Model without depredation"

indices_comp <- bind_rows(index_std, index)

ggplot(indices_comp, 
       aes(x = TAAFseason, y = est, colour = model, group = model)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = model),
              alpha = 0.2, colour = NA) +
  theme_bw() +
  labs(
    x = "Year",
    y = "Standardized CPUE",
    colour = "Model",
    fill = "Model"
  )






grid_comp <- pred_std$data %>%
  mutate(
    est_depred = exp(est),
    est_classique = exp(pred$data$est)
  )

grid_comp <- grid_comp %>%
  mutate(
    diff = est_depred - est_classique,
    ratio = est_depred / est_classique,
    diff_pct = 100 * (est_depred - est_classique) / pmax(est_classique, 1e-6),
    log_ratio = log1p(est_depred) - log1p(est_classique)
  )




map_mean <- grid_comp %>%
  group_by(cell_id, X, Y) %>%
  summarise(
    mean_diff_pct = mean(diff_pct, na.rm = TRUE),
    .groups = "drop"
  )


ggplot(map_mean,
       aes(X, Y, fill = mean_diff_pct)) +
  geom_tile() +
  coord_equal() +
  scale_fill_gradient2(
    midpoint = 0,
    low = "blue",
    mid = "white",
    high = "red",
    name = "% différence"
  ) +
  theme_bw()


map_mean <- map_mean %>%
  mutate(
    X = X * 1000,
    Y = Y * 1000
  )
map_mean_sf <- st_as_sf(map_mean, coords = c("X","Y"), crs = 32742)

map_mean_sf <- st_transform(map_mean_sf, 4326)




KER_map + 
  geom_sf(
    data = map_mean_sf, 
    aes(color = -mean_diff_pct), 
    size = 0.1,
    inherit.aes = FALSE
  ) + 
  scale_color_gradient2( 
    midpoint = 0, 
    low = "red", 
    mid = "white", 
    high = "blue",
    name = "Difference of CPUE in %"
  ) +
  theme(
    legend.position = c(0.85, 0.09),   # <- position dans la carte (x, y)
    legend.background = element_rect(fill = "white", colour = "black"),
    legend.key = element_rect(fill = "white")
  ) +
  guides(
    color = guide_colorbar(
      direction = "vertical",
      barheight = unit(1, "cm"),
      barwidth = unit(0.2, "cm")
    )
  ) +
  geom_sf(
    data = Isobath_KER_sf %>% dplyr::filter(level == "-500"),
    colour = "grey80",
    size = 0.5
  ) +
  geom_sf(
    data = Isobath_KER_sf %>% dplyr::filter(level == "-1000"),
    colour = "grey60",
    size = 0.5
  ) +
  geom_sf(
    data = Isobath_KER_sf %>% dplyr::filter(level == "-1500"),
    colour = "grey30",
    size = 0.5
  ) +
  geom_sf(
    data = Isobath_KER_sf %>% dplyr::filter(level == "-2000"),
    colour = "black",
    size = 0.5
  )

