
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




cpue_summary <- dat_pos %>%
  filter(TAAFseason >= 2010, TAAFseason <= 2025) %>%
  group_by(TAAFseason) %>%
  summarise(
    n = sum(!is.na(CPUE)),
    mean_CPUE = mean(CPUE, na.rm = TRUE),
    sd_CPUE = sd(CPUE, na.rm = TRUE),
    se_CPUE = sd_CPUE / sqrt(n),
    lwr = mean_CPUE - 1.96 * se_CPUE,
    upr = mean_CPUE + 1.96 * se_CPUE,
    .groups = "drop"
  )

ggplot(cpue_summary, aes(TAAFseason, mean_CPUE)) +
  geom_ribbon(
    aes(ymin = lwr, ymax = upr),
    alpha = 0.2
  ) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  theme_bw() +
  labs(
    x = "Saison de pêche",
    y = "CPUE moyenne en g/hameçon"
  )

cpue_summary <- cpue_summary %>%
  arrange(TAAFseason)

last_year <- max(cpue_summary$TAAFseason)

cpue_2025 <- cpue_summary %>%
  filter(TAAFseason == last_year) %>%
  pull(mean_CPUE)

cpue_2024 <- cpue_summary %>%
  filter(TAAFseason == last_year - 1) %>%
  pull(mean_CPUE)

pct_24_25 <- 100 * (cpue_2025 - cpue_2024) / cpue_2024
pct_24_25

mean_4y <- cpue_summary %>%
  filter(TAAFseason >= (last_year - 4),
         TAAFseason < last_year) %>%
  summarise(m = mean(mean_CPUE, na.rm = TRUE)) %>%
  pull(m)

pct_4y <- 100 * (cpue_2025 - mean_4y) / mean_4y
pct_4y

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
  ScanWeightKg ~ s(TAAFseason) + s(Month, bs = "cc") + Vessel + s(prof_bathy) + s(slope) + s(HookLength) +
    cachalot_obs + orque_obs + s(NB_MEAN_CACHALOT, NB_MEAN_ORQUE),
  offset = log(dat_pos$nb_ham * 0.246),
  data = dat_pos,
  mesh = bspde_pos,
  family = Gamma(link = "log"),
  spatial = "on",
  time= "TAAFseason",
  spatiotemporal = "iid"
)



m_2 <- sdmTMB(
  ScanWeightKg ~ s(TAAFseason) + s(Month, bs = "cc") + Vessel + s(prof_bathy) + s(slope) + s(HookLength),
  offset = log(dat_pos$nb_ham * 0.246),
  data = dat_pos,
  mesh = bspde_pos,
  family = Gamma(link = "log"),
  spatial = "on",
  time= "TAAFseason",
  spatiotemporal = "iid"
)

visreg(m_1, xvar = "HookLength", gg = TRUE)
visreg(m_1, xvar = "prof_bathy", gg = TRUE)
visreg(m_1, xvar = "slope", gg = TRUE)
visreg(m_1, xvar = "Vessel", gg = TRUE)
visreg2d(m_1, xvar = "prof_bathy", yvar = "slope", plot.type = "image")



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

#enlever mois puis peut etre ensuite trouver biblio pour justifier (periode ou ils ont plus faim etc)
m_season_varying_2 <- sdmTMB(
  ScanWeightKg ~ s(TAAFseason) +
    s(Month, bs = "cc")  +
    cachalot_obs + orque_obs + s(NB_MEAN_CACHALOT, NB_MEAN_ORQUE) ,
  offset = log(dat_pos$nb_ham * 0.246),
  data = dat_pos,
  mesh = bspde_pos,
  family = Gamma(link = "log"),
  spatial = "on",
  time = "TAAFseason",
  spatiotemporal = "iid"
)


m_season_varying_3 <- sdmTMB(
  ScanWeightKg ~ s(TAAFseason) +
    s(Month, bs = "cc")  +
    cachalot_obs + orque_obs + s(NB_MEAN_CACHALOT, NB_MEAN_ORQUE)  + Vessel,
  offset = log(dat_pos$nb_ham * 0.246),
  data = dat_pos,
  mesh = bspde_pos,
  family = Gamma(link = "log"),
  spatial = "on",
  time = "TAAFseason",
  spatiotemporal = "iid"
)


m_season_varying_4 <- sdmTMB(
  ScanWeightKg ~ s(TAAFseason) +
    s(Month, bs = "cc")  +
    cachalot_obs + orque_obs + s(NB_MEAN_CACHALOT, NB_MEAN_ORQUE) + Vessel + s(tps_trempe_hm_num),
  offset = log(dat_pos$nb_ham * 0.246),
  data = dat_pos,
  mesh = bspde_pos,
  family = Gamma(link = "log"),
  spatial = "on",
  time = "TAAFseason",
  spatiotemporal = "iid"
)


m_season_varying_5 <- sdmTMB(
  ScanWeightKg ~ s(TAAFseason) +
    s(Month, bs = "cc")  +
    cachalot_obs + orque_obs + s(NB_MEAN_CACHALOT, NB_MEAN_ORQUE) + Vessel + s(prof_bathy) + s(tps_trempe_hm_num),
  offset = log(dat_pos$nb_ham * 0.246),
  data = dat_pos,
  mesh = bspde_pos,
  family = Gamma(link = "log"),
  spatial = "on",
  time = "TAAFseason",
  spatiotemporal = "iid"
)



m_null <- sdmTMB(
  ScanWeightKg ~ 1,
  offset = log(dat_pos$nb_ham * 0.246),
  data = dat_pos,
  mesh = bspde_pos,
  family = Gamma(link = "log")
)



#Comparaison des AIC :
AIC(m_2)
AIC(m_season_varying)
AIC(m_season_varying_inter)
AIC(m_season_varying_1)
AIC(m_season_varying_2) 
AIC(m_season_varying_3)
AIC(m_season_varying_4)
AIC(m_season_varying_5)

library(visreg)

#Visualitation des effets
visreg(m_season_varying_5, xvar = "cachalot_obs", gg = TRUE)
visreg(m_season_varying_5, xvar = "orque_obs", gg = TRUE)
visreg(m_season_varying_5, xvar = "NB_MEAN_CACHALOT", gg = TRUE)
visreg(m_season_varying_5, xvar = "NB_MEAN_ORQUE", gg = TRUE)
visreg(m_season_varying_5, xvar = "Month", gg = TRUE)
visreg(m_season_varying_5, xvar = "Vessel", gg = TRUE)
visreg(m_season_varying_5, xvar = "prof_bathy", gg = TRUE)
visreg(m_season_varying_5, xvar = "tps_trempe_hm_num", gg = TRUE)
visreg2d(m_season_varying_5, xvar = "NB_MEAN_CACHALOT", yvar = "NB_MEAN_ORQUE", plot.type = "image")

V <- m_season_varying_5$sd_report$cov.fixed
Z <- m_season_varying_5$smoothers$Zs[[1]]
Z <- Z[, 1:8]  # important !
beta <- m_season_varying_5$parlist$b_smooth[,1]
beta_taaf <- beta[1:8]

fit <- as.numeric(Z %*% beta_taaf)
var_fit <- diag(Z %*% V[1:8, 1:8] %*% t(Z))
se_fit <- sqrt(var_fit)
upper <- fit + 1.96 * se_fit
lower <- fit - 1.96 * se_fit
eff <- data.frame(
  TAAFseason = dat_pos$TAAFseason,
  fit = fit,
  lower = lower,
  upper = upper
) %>%
  group_by(TAAFseason) %>%
  summarise(
    fit = mean(fit),
    lower = mean(lower),
    upper = mean(upper),
    .groups = "drop"
  )

ggplot(eff, aes(TAAFseason, fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(linewidth = 1) +
  theme_bw() +
  labs(
    x = "TAAFseason",
    y = "Effet s(TAAFseason)",
    title = "Smooth temporel avec IC 95%"
  )



visreg(m_season_varying, xvar = "NB_MEAN_CACHALOT", gg = TRUE)
visreg(m_season_varying, xvar = "NB_MEAN_ORQUE", gg = TRUE)





# Contrôle des modèles

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
  m_season_varying_2,
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
  m_season_varying_2,
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
  group_by(TAAFseason) %>%
  summarise(
    est = weighted.mean(est, w = nb_ham, na.rm = TRUE),
    lwr = weighted.mean(lwr, w = nb_ham, na.rm = TRUE),
    upr = weighted.mean(upr, w = nb_ham, na.rm = TRUE)
  ) %>%
  mutate(scenario = "Avec déprédation")

index_nodep_year <- index_nodep %>%
  group_by(TAAFseason) %>%
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
  aes(TAAFseason, est, color = scenario, fill = scenario)
) +
  geom_ribbon(
    aes(ymin = lwr, ymax = upr),
    alpha = 0.2,
    color = NA
  ) +
  geom_line(linewidth = 1.2) +
  labs(
    x = "TAAFseason",
    y = "Biomasse prédite",
    color = "Scénario",
    fill = "Scénario"
  ) +
  theme_bw()

Abundance_Index


index_loss_year <- index_dep_year %>%
  left_join(
    index_nodep_year,
    by = "TAAFseason",
    suffix = c("_dep", "_nodep")
  ) %>%
  mutate(
    perte_pct = est_nodep - est_dep
  )

ggplot(index_loss_year, aes(x = TAAFseason, y = perte_pct)) +
  geom_line(linewidth = 1.2) +
  geom_point() +
  theme_bw() +
  labs(
    x = "Année",
    y = "Biomasse déprédatée",
    title = "Impact annuel de la déprédation"
  )













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

# Sans orques uniquement
dat_sans_orque <- dat_pos %>%
  mutate(
    NB_MEAN_ORQUE = 0
  )

# Sans cachalots uniquement
dat_sans_cachalot <- dat_pos %>%
  mutate(
    NB_MEAN_CACHALOT = 0
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

pred_sans_orque <- predict(
  m_1,
  newdata = dat_sans_orque,
  offset = offset_sans_orque,
  type = "response"
)

pred_sans_cachalot <- predict(
  m_1,
  newdata = dat_sans_cachalot,
  offset = offset_sans_cachalot,
  type = "response"
)

# 5. Dataframe final
# ------------------

df_station_predictions <- dat_pos %>%
  mutate(
    poids_reel = ScanWeightKg,
    
    pred_modele = pred_reel$est,
    
    pred_sans_dep = pred_sans_dep$est,
    
    pred_sans_orque = pred_sans_orque$est,
    
    pred_sans_cachalot = pred_sans_cachalot$est,
    
    perte_totale = pred_sans_dep - pred_modele,
    
    perte_orque = pred_sans_orque - pred_modele,
    
    perte_cachalot = pred_sans_cachalot - pred_modele
  ) %>%
  select(
    StationID,
    TAAFseason,
    Month,
    ScanWeightKg,
    nb_ham,
    NB_MEAN_CACHALOT,
    NB_MEAN_ORQUE,
    X,
    Y,
    cachalot_obs,
    orque_obs,
    poids_reel,
    pred_modele,
    pred_sans_dep,
    pred_sans_orque,
    pred_sans_cachalot,
    perte_totale,
    perte_orque,
    perte_cachalot
  )

# Résultat
head(df_station_predictions)

saveRDS(df_station_predictions,
        "prediction_depredation_crozet.rds")

# Création des pourcentages
df_boxplot_pct <- df_station_predictions %>%
  mutate(
    perte_totale_pct =
      100 * (pred_sans_dep - pred_modele) / pred_sans_dep,
    
    perte_orque_pct =
      100 * (pred_sans_dep - pred_sans_cachalot) / pred_sans_dep,
    
    perte_cachalot_pct =
      100 * (pred_sans_dep - pred_sans_orque) / pred_sans_dep
  ) %>%
  
  select(
    TAAFseason,
    perte_totale_pct,
    perte_orque_pct,
    perte_cachalot_pct
  ) %>%
  
  pivot_longer(
    cols = -TAAFseason,
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


# -> faire plutot la somme dans le carré et après le pourcentage
# -> interessant sur la résultante déprédaté
# -> peut essayer un kde
# -> série par groupe de 5 ans
# -> peut projeter en retirant effet fixe pour l'intéraction '
# -> sortir tableau et envoyer
# -> 



df_summary <- df_boxplot_pct %>%
  group_by(TAAFseason, type_perte) %>%
  summarise(
    mean_perte = mean(pourcentage_perte, na.rm = TRUE),
    q25 = quantile(pourcentage_perte, 0.25, na.rm = TRUE),
    q75 = quantile(pourcentage_perte, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(df_summary,
       aes(x = TAAFseason,
           y = mean_perte,
           color = type_perte)) +
  
  geom_line(linewidth = 1.1) +
  
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




offset_pos <- log(0.246 * dat_pos$nb_ham)
pred_reel <- predict(
  m_1,
  newdata = dat_pos,
  offset = offset_pos,  # si on veut une CPUE en sortie
  type = "response",
  return_tmb_object = TRUE
)

pred_sans_dep <- predict(
  m_2,
  newdata = dat_pos,
  offset = offset_pos,
  type = "response",
  return_tmb_object = TRUE
)

idx_reel <- get_index(pred_reel)

idx_sans_dep <- get_index(pred_sans_dep)


idx_plot <- bind_rows(
  idx_reel %>% mutate(scenario = "Réel"),
  idx_sans_dep %>% mutate(scenario = "Sans déprédation")
)





map = terra::rast(here(path_rast_bathy, 'ras_pecheker_interpol_totCRO.tif'))


map_grid_m <- terra::project(map, "EPSG:32743")
print(map_grid_m)

res <- 5000

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
# df_station_predictions <- df_station_predictions %>%
#   mutate(
#     perte_totale_pct =
#       100 * (pred_sans_dep - pred_modele) / pred_sans_dep,
#     
#     perte_orque_pct =
#       100 * perte_orque / pred_sans_dep,
#     
#     perte_cachalot_pct =
#       100 * perte_cachalot / pred_sans_dep
#   )

pts <- vect(df_station_predictions, geom = c("X", "Y"), crs = crs(map_grid_1km))

df_station_predictions$cell_id <- extract(r_id, pts)[,2]
sum(is.na(df_station_predictions$cell_id))
loss_grid <- df_station_predictions %>%
  group_by(cell_id) %>%
  summarise(
    perte_totale_pct_mean = 100*sum(perte_totale, na.rm = TRUE)/sum(pred_sans_dep, na.rm = TRUE),
    perte_orque_pct_mean = 100*sum(perte_orque, na.rm = TRUE)/sum(pred_sans_dep, na.rm = TRUE),
    perte_cachalot_pct_mean = 100*sum(perte_cachalot, na.rm = TRUE)/sum(pred_sans_dep, na.rm = TRUE)
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
    xlim = c(-1550000, -1150000),
    ylim = c(4505000, 4700000)
  )



ggplot(grid_final) +
  geom_raster(aes(x = x, y = y, fill = perte_cachalot_pct_mean)) +
  scale_fill_viridis_c(na.value = "transparent") +
  coord_equal() +
  coord_cartesian(
    xlim = c(-1550000, -1150000),
    ylim = c(4505000, 4700000)
  )



ggplot(grid_final) +
  geom_raster(aes(x = x, y = y, fill = perte_orque_pct_mean)) +
  scale_fill_viridis_c(na.value = "transparent") +
  coord_equal() +
  coord_cartesian(
    xlim = c(-1550000, -1150000),
    ylim = c(4505000, 4700000)
  )







#Calcul des CPUE par périodes, comparaison pour observer la baisse

nb_ham_std <- 1
offset_std <- log(nb_ham_std * 0.246)

pred_std <- predict(
  m_season_varying_2,
  newdata = dat_pos,
  offset = rep(offset_std, nrow(dat_pos)),
  type = "response"
)

dat_pos$CPUE_std <- pred_std$est * 1000


dat_pos <- dat_pos %>%
  mutate(
    periode = case_when(
      Year %in% c(2021, 2022) ~ "2021-2022",
      Year %in% c(2022, 2023) ~ "2022-2023",
      Year %in% c(2023, 2024) ~ "2023-2024",
      TRUE ~ NA_character_
    )
  )


cpue_period <- dat_pos %>%
  filter(!is.na(periode)) %>%
  group_by(periode) %>%
  summarise(
    n = n(),
    mean_cpue = mean(CPUE_std, na.rm = TRUE),
    sd_cpue = sd(CPUE_std, na.rm = TRUE),
    se = sd_cpue / sqrt(n),
    IC_inf = mean_cpue - qt(0.975, df = n - 1) * se,
    IC_sup = mean_cpue + qt(0.975, df = n - 1) * se,
    .groups = "drop"
  )


ggplot(cpue_period, aes(x = periode, y = mean_cpue)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = IC_inf, ymax = IC_sup), width = 0.15) +
  theme_classic() +
  labs(x = "", y = "CPUE standardisée (moyenne)")


ref <- cpue_period$mean_cpue[cpue_period$periode == "2021-2022"]

cpue_period <- cpue_period %>%
  mutate(
    index_rel = 100 * mean_cpue / ref,
    pct_change = 100 * (mean_cpue / ref - 1)
  )




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
  tidyr::uncount(length(years)) %>%
  mutate(TAAFseason = rep(years, times = nrow(grid_with_id))) %>%
  mutate(cachalot_obs = 1,
         Month = 2,
         orque_obs = 1,
         nb_ham = 1,
         NB_MEAN_CACHALOT = 0,
         NB_MEAN_ORQUE = 0,
         Vessel = "Albius",
         HookLength = mean(dat_pos$HookLength)) %>%
  relocate(cell_id)

#Vérifier colonnes requises
required_cols <- c("prof_bathy", "TAAFseason", "X", "Y")
missing_cols <- setdiff(required_cols, names(grid_all_years))
if (length(missing_cols) > 0) stop("Colonnes manquantes dans newdata : ", paste(missing_cols, collapse = ", "))



index <- get_index_split(
  m_1,
  newdata = grid_all_years,
  area = grid_all_years$nb_ham * 0.246,
  bias_correct = FALSE,
  nsplit = 8
)


index_2 <- get_index_split(
  m_2,
  newdata = grid_all_years,
  area = grid_all_years$nb_ham * 0.246,
  bias_correct = FALSE,
  nsplit = 8
)



index %>%
  filter(TAAFseason >= 2010, TAAFseason <= 2025) %>%
  ggplot(aes(x = TAAFseason, y = est)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(ymin = lwr, ymax = upr),
    width = 0.2
  ) +
  scale_x_continuous(breaks = 2010:2025) +
  labs(
    x = "Saison",
    y = "Indice standardisé",
    title = "Évolution de l'indice (2010–2025)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


ggplot(
  filter(index, TAAFseason >= 2010, TAAFseason <= 2025),
  aes(TAAFseason, est)
) +
  geom_ribbon(
    aes(ymin = lwr, ymax = upr),
    alpha = 0.2
  ) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  theme_bw() +
  labs(
    x = "Saison de pêche",
    y = "Indice standardisé de la CPUE"
  )

index_comp <- index %>%
  filter(TAAFseason >= 2010, TAAFseason <= 2025) %>%
  arrange(TAAFseason)

last_year <- max(index_comp$TAAFseason)

est_last <- index_comp %>%
  filter(TAAFseason == last_year) %>%
  pull(est)

est_prev <- index_comp %>%
  filter(TAAFseason == last_year - 1) %>%
  pull(est)

pct_vs_prev_year <- 100 * (est_last - est_prev) / est_prev

pct_vs_prev_year


index_comp <- index %>%
  filter(TAAFseason >= 2010, TAAFseason <= 2025) %>%
  arrange(TAAFseason)

last_year <- max(index_comp$TAAFseason)

est_last <- index_comp %>%
  filter(TAAFseason == last_year) %>%
  pull(est)

mean_4y <- index_comp %>%
  filter(TAAFseason >= (last_year - 4), TAAFseason < last_year) %>%
  summarise(m = mean(est)) %>%
  pull(m)

pct_vs_4y <- 100 * (est_last - mean_4y) / mean_4y

pct_vs_4y





beta <- m_season_varying_2$parlist$b_smooth[,1]

beta_taaf <- beta[1:8]

Z_taaf <- m_season_varying_2$smoothers$Zs[[1]]
attr(Z_taaf, "s.label")
# "s(TAAFseason)"

f_taaf <- Z_taaf %*% beta_taaf
eff_taaf <- data.frame(
  TAAFseason = dat_pos$TAAFseason,
  fit = as.numeric(f_taaf)
)
eff_taaf <- eff_taaf %>%
  group_by(TAAFseason) %>%
  summarise(fit = mean(fit), .groups = "drop")
library(ggplot2)

ggplot(eff_taaf, aes(TAAFseason, fit)) +
  geom_line(linewidth = 1) +
  geom_point() +
  theme_bw() +
  labs(
    x = "TAAFseason",
    y = "Effet partiel s(TAAFseason)"
  )

eff_taaf <- eff_taaf %>%
  arrange(TAAFseason)

last_year <- max(eff_taaf$TAAFseason)

fit_2025 <- eff_taaf %>%
  filter(TAAFseason == last_year) %>%
  pull(fit)

fit_2024 <- eff_taaf %>%
  filter(TAAFseason == last_year - 1) %>%
  pull(fit)

pct_24_25 <- 100 * (fit_2025 - fit_2024) / abs(fit_2024)
pct_24_25

mean_4y <- eff_taaf %>%
  filter(TAAFseason >= (last_year - 4),
         TAAFseason < last_year) %>%
  summarise(m = mean(fit, na.rm = TRUE)) %>%
  pull(m)

pct_4y <- 100 * (fit_2025 - mean_4y) / abs(mean_4y)
pct_4y











smooth_sdmTMB <- function(model, term, se = FALSE) {
  
  # labels des smooths
  labs <- unlist(model$smoothers$labels)
  
  # index du smooth demandé
  i <- which(labs == paste0("s(", term, ")"))
  if (length(i) == 0) {
    stop("Term not found. Available: ",
         paste(gsub("s\\(|\\)", "", labs), collapse = ", "))
  }
  
  # matrice Z du smooth
  Z <- model$smoothers$Zs[[i]]
  
  # indices des coefficients dans b_smooth
  starts <- model$smoothers$b_smooth_start + 1
  dims   <- model$smoothers$sm_dims
  idx <- starts[i]:(starts[i] + dims[i] - 1)
  
  # coefficients
  beta <- model$parlist$b_smooth[,1]
  b <- beta[idx]
  
  # fit
  fit <- as.numeric(Z %*% b)
  
  out <- data.frame(
    fit = fit
  )
  
  # SE (approx uniquement si demandé)
  if (se) {
    
    V <- model$sd_report$cov.fixed
    
    # ⚠️ approximation (car V ≠ espace b_smooth)
    Vb <- diag(length(b)) * mean(diag(V))
    
    se_fit <- sqrt(diag(Z %*% Vb %*% t(Z)))
    
    out$se <- se_fit
    out$lower <- fit - 1.96 * se_fit
    out$upper <- fit + 1.96 * se_fit
  }
  
  return(out)
}

eff <- smooth_sdmTMB(m_season_varying_2, "prof_bathy", se = TRUE)
eff$prof_bathy <- dat_pos$prof_bathy

ggplot(eff, aes(prof_bathy, fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(linewidth = 1)  +
  theme_bw() +
  labs(
    x = "TAAFseason",
    y = "Effet s(TAAFseason)",
    title = "Smooth temporel avec IC 95%"
  )






eff <- smooth_sdmTMB(m_season_varying_2, "tps_trempe_hm_num", se = TRUE)
eff$tps_trempe_hm_num <- dat_pos$tps_trempe_hm_num

ggplot(eff, aes(tps_trempe_hm_num, fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(linewidth = 1)  +
  theme_bw() +
  labs(
    x = "TAAFseason",
    y = "Effet s(TAAFseason)",
    title = "Smooth temporel avec IC 95%"
  )




eff <- smooth_sdmTMB(m_season_varying_2, "Month", se = TRUE)
eff$Month <- dat_pos$Month

ggplot(eff, aes(Month, fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(linewidth = 1)  +
  theme_bw() +
  labs(
    x = "TAAFseason",
    y = "Effet s(TAAFseason)",
    title = "Smooth temporel avec IC 95%"
  )



eff <- smooth_sdmTMB(m_season_varying_2, "TAAFseason", se = TRUE)
eff$TAAFseason <- dat_pos$TAAFseason

ggplot(eff, aes(TAAFseason, fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(linewidth = 1)  +
  theme_bw() +
  labs(
    x = "TAAFseason",
    y = "Effet s(TAAFseason)",
    title = "Smooth temporel avec IC 95%"
  )




eff_cach <- data.frame(
  level = c("absence","non obsé","présence"),
  beta = c(
    0,
    coef(m_season_varying_2)["cachalot_obsnon obsé"],
    coef(m_season_varying_2)["cachalot_obsprésence"]
  )
)
ggplot(eff_cach, aes(level, beta)) +
  geom_col() +
  theme_bw() +
  labs(x = "cachalot_obs", y = "Effet")


eff_orque <- data.frame(
  level = c("absence","non obsé","présence"),
  beta = c(
    0,
    coef(m_season_varying_2)["orque_obsnon obsé"],
    coef(m_season_varying_2)["orque_obsprésence"]
  )
)
ggplot(eff_orque, aes(level, beta)) +
  geom_col() +
  theme_bw() +
  labs(x = "orque_obs", y = "Effet")
