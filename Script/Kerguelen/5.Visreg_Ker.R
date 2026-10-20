
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


vars_std <- c(
  "TAAFseason",
  "prof_bathy",
  "slope",
  "HookLength",
  "NB_MEAN_CACHALOT",
  "NB_MEAN_ORQUE",
  "bottomT"
)

stats_std <- dat_pos %>%
  summarise(
    across(
      all_of(vars_std),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd   = ~sd(.x, na.rm = TRUE)
      )
    )
  )


nb_ham_med <- median(dat_pos$nb_ham, na.rm = TRUE)

offset_med <- log(
  0.56 * 0.002 * (nb_ham_med - 1) + 0.246
)

# Pour chaque variable standardisé on :
#   Fait une prédiction via visreg
#   Remet la variable sur l échelle réelle
#   crée le graphique de l effet

#Hameçon
H <- visreg(
  m_1,
  "HookLength_sc",
  plot = FALSE
)

y_ref <- H$fit$visregFit[1]

H$fit$indice_relatif <-
  H$fit$visregFit / y_ref - 1

H$fit$indice_lwr <-
  H$fit$visregLwr / y_ref - 1

H$fit$indice_upr <-
  H$fit$visregUpr / y_ref - 1

H$fit$indice_pct <-
  100 * (H$fit$visregFit / y_ref - 1)

H$fit$indice_lwr_pct <-
  100 * (H$fit$visregLwr / y_ref - 1)

H$fit$indice_upr_pct <-
  100 * (H$fit$visregUpr / y_ref - 1)

H$fit$HookLength_sc_real <-
  H$fit$HookLength_sc *
  stats_std$HookLength_sd +
  stats_std$HookLength_mean


p_hook <- ggplot(H$fit, aes(x = HookLength_sc_real, y = indice_pct)) +
  geom_ribbon(
    aes(
      ymin = indice_lwr_pct,
      ymax = indice_upr_pct
    ),
    alpha = 0.2
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_rug(
    data = dat_pos_std,
    aes(x = HookLength),
    inherit.aes = FALSE,
    sides = "b",
    alpha = 0.2
  )+
  geom_line(linewidth = 1) +
  labs(
    x = "Longueur d'hameçon",
    y = "Variation relative de la CPUE (%)"
  ) +
  theme_bw()


#bathymétrie
B <- visreg(
  m_1,
  "prof_bathy_sc",
  plot = FALSE,
  predict = list(offset = offset_med)
)

# Retour à l'échelle réelle
B$fit$prof_bathy_real <-
  B$fit$prof_bathy_sc *
  stats_std$prof_bathy_sd +
  stats_std$prof_bathy_mean

# Valeur de référence = première prédiction
y_ref <- B$fit$visregFit[1]

# Indice relatif (%)
B$fit$indice_pct <-
  100 * (B$fit$visregFit / y_ref - 1)

B$fit$indice_lwr_pct <-
  100 * (B$fit$visregLwr / y_ref - 1)

B$fit$indice_upr_pct <-
  100 * (B$fit$visregUpr / y_ref - 1)

# Graphique
p_bathy <- ggplot(
  B$fit,
  aes(x = prof_bathy_real, y = indice_pct)
) +
  geom_ribbon(
    aes(
      ymin = indice_lwr_pct,
      ymax = indice_upr_pct
    ),
    alpha = 0.2
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_line(linewidth = 1) +
  geom_rug(
    data = dat_pos_std,
    aes(x = prof_bathy),
    inherit.aes = FALSE,
    sides = "b",
    alpha = 0.2
  ) +
  labs(
    x = "Bathymétrie (m)",
    y = "Variation relative de la CPUE (%)"
  ) +
  theme_bw()


#slope
S <- visreg(
  m_1,
  "slope_sc",
  plot = FALSE,
  predict = list(offset = offset_med)
)

# Retour à l'échelle réelle
S$fit$slope_real <-
  S$fit$slope_sc *
  stats_std$slope_sd +
  stats_std$slope_mean

# Valeur de référence = première prédiction
y_ref <- S$fit$visregFit[1]

# Indice relatif (%)
S$fit$indice_pct <-
  100 * (S$fit$visregFit / y_ref - 1)

S$fit$indice_lwr_pct <-
  100 * (S$fit$visregLwr / y_ref - 1)

S$fit$indice_upr_pct <-
  100 * (S$fit$visregUpr / y_ref - 1)

# Graphique
p_slope <- ggplot(
  S$fit,
  aes(x = slope_real, y = indice_pct)
) +
  geom_ribbon(
    aes(
      ymin = indice_lwr_pct,
      ymax = indice_upr_pct
    ),
    alpha = 0.2
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_line(linewidth = 1) +
  geom_rug(
    data = dat_pos_std,
    aes(x = slope),
    inherit.aes = FALSE,
    sides = "b",
    alpha = 0.2
  ) +
  labs(
    x = "Pente",
    y = "Variation relative de la CPUE (%)"
  ) +
  theme_bw()



#Cachalot
C <- visreg(
  m_1,
  "NB_MEAN_CACHALOT_sc",
  plot = FALSE,
  predict = list(offset = offset_med)
)

# Retour à l'échelle réelle
C$fit$NB_MEAN_CACHALOT_real <-
  C$fit$NB_MEAN_CACHALOT_sc *
  stats_std$NB_MEAN_CACHALOT_sd +
  stats_std$NB_MEAN_CACHALOT_mean

# Valeur de référence = première prédiction
y_ref <- C$fit$visregFit[1]

# Indice relatif (%)
C$fit$indice_pct <-
  100 * (C$fit$visregFit / y_ref - 1)

C$fit$indice_lwr_pct <-
  100 * (C$fit$visregLwr / y_ref - 1)

C$fit$indice_upr_pct <-
  100 * (C$fit$visregUpr / y_ref - 1)

# Graphique
p_cachalot <- ggplot(
  C$fit,
  aes(x = NB_MEAN_CACHALOT_real, y = indice_pct)
) +
  geom_ribbon(
    aes(
      ymin = indice_lwr_pct,
      ymax = indice_upr_pct
    ),
    alpha = 0.2
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_line(linewidth = 1) +
  geom_rug(
    data = dat_pos_std,
    aes(x = NB_MEAN_CACHALOT),
    inherit.aes = FALSE,
    sides = "b",
    alpha = 0.2
  ) +
  labs(
    x = "Nombre de cachalots",
    y = "Variation relative de la CPUE (%)"
  ) +
  theme_bw()

#Année
A <- visreg(
  m_1,
  "TAAFseason_sc",
  plot = FALSE,
  predict = list(offset = offset_med)
)

# Retour à l'échelle réelle
A$fit$TAAFseason_real <-
  A$fit$TAAFseason_sc *
  stats_std$TAAFseason_sd +
  stats_std$TAAFseason_mean

# Valeur de référence = première prédiction
y_ref <- A$fit$visregFit[1]

# Indice relatif (%)
A$fit$indice_pct <-
  100 * (A$fit$visregFit / y_ref - 1)

A$fit$indice_lwr_pct <-
  100 * (A$fit$visregLwr / y_ref - 1)

A$fit$indice_upr_pct <-
  100 * (A$fit$visregUpr / y_ref - 1)

# Graphique
p_annee <- ggplot(
  A$fit,
  aes(x = TAAFseason_real, y = indice_pct)
) +
  geom_ribbon(
    aes(
      ymin = indice_lwr_pct,
      ymax = indice_upr_pct
    ),
    alpha = 0.2
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_line(linewidth = 1) +
  geom_rug(
    data = dat_pos_std,
    aes(x = TAAFseason),
    inherit.aes = FALSE,
    sides = "b",
    alpha = 0.2
  ) +
  labs(
    x = "Année",
    y = "Variation relative de la CPUE (%)"
  ) +
  theme_bw()


#Mois
M <- visreg(
  m_1,
  "Month",
  plot = FALSE,
  predict = list(offset = offset_med)
)

# Valeur de référence = première prédiction
y_ref <- M$fit$visregFit[1]

# Indice relatif (%)
M$fit$indice_pct <-
  100 * (M$fit$visregFit / y_ref - 1)

M$fit$indice_lwr_pct <-
  100 * (M$fit$visregLwr / y_ref - 1)

M$fit$indice_upr_pct <-
  100 * (M$fit$visregUpr / y_ref - 1)

# Graphique
p_month <- ggplot(
  M$fit,
  aes(x = Month, y = indice_pct)
) +
  geom_ribbon(
    aes(
      ymin = indice_lwr_pct,
      ymax = indice_upr_pct
    ),
    alpha = 0.2
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_line(linewidth = 1) +
  scale_x_continuous(
    breaks = 1:12
  ) +
  labs(
    x = "Mois",
    y = "Variation relative de la CPUE (%)"
  ) +
  theme_bw()

#Vessel
# Vessel
V <- visreg(
  m_1,
  "Vessel",
  plot = FALSE,
  predict = list(offset = offset_med)
)

# Modalité de référence = première modalité
y_ref <- V$fit$visregFit[1]

# Indice relatif (%)
V$fit$indice_pct <-
  100 * (V$fit$visregFit / y_ref - 1)

V$fit$indice_lwr_pct <-
  100 * (V$fit$visregLwr / y_ref - 1)

V$fit$indice_upr_pct <-
  100 * (V$fit$visregUpr / y_ref - 1)

# Anonymisation des navires
V$fit$Vessel_anonyme <- LETTERS[
  seq_along(unique(V$fit$Vessel))
]

p_vessel <- ggplot(
  V$fit,
  aes(x = Vessel_anonyme, y = indice_pct)
) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(
      ymin = indice_lwr_pct,
      ymax = indice_upr_pct
    ),
    width = 0.2
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    x = "Navire",
    y = "Variation relative de la CPUE (%)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


#Création des groupement de figure en 2 catégories (variable de capturabilité ou de densité)

figure_effets_catch <- (
  p_cachalot +
    p_month +
    p_hook +
    p_vessel 
) +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "A") &
  theme(
    axis.title.y = element_text(size = 8),
    axis.title.x = element_text(size = 9)
  )



figure_effets_dens <- (
  p_bathy +
    p_slope +
    p_annee )+
  plot_layout(ncol = 3) +
  plot_annotation(tag_levels = "A")



#Sauvegarde

ggsave(
  here("data", output_plot, "figure_effets_catch_Ker.jpg"),
  figure_effets_catch,
  width = 12,
  height = 8,
  units = "in",
  dpi = 300
)

ggsave(
  here("data", output_plot, "figure_effets_dens_Ker.jpg"),
  figure_effets_dens,
  width = 12,
  height = 4.8,
  units = "in",
  dpi = 300
)
