
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

dat_pres <- dat_clean %>%
  mutate(pres = as.integer(CPUE > 0))

dat_pos <- dat_pres %>% filter(CPUE > 0, prof_bathy >= 500)

# Calcul du 95e percentile de la CPUE
cpue_95 <- quantile(dat_pos$CPUE, probs = 0.95, na.rm = TRUE)

# Filtrer les données pour ne garder que celles <= 95e percentile
dat_pos <- dat_pos %>%
  filter(CPUE <= cpue_95)


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
    xlim = range(dat_pos$X) + c(-30, 30),
    ylim = range(dat_pos$Y) + c(-30, 30)
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


m_0 <- sdmTMB(
  CPUE ~ s(Year),
  data = dat_pos,
  mesh = bspde_pos,
  family = Gamma(link = "log"),
  spatial = "on",
  time= "Year",
  spatiotemporal = "iid"
)







# 4. Mesh -----------------------------------------------------------------

# Création du maillage spatial pour le modèle (choix du nombre de nœuds)
noeuds <- 120


# Supprime les lignes avec CPUE manquante
dat_clean <- dat %>% filter(!is.na(CPUE), !is.na(prof_bathy), !is.na(Year))

dat_sf <- st_as_sf(dat_clean, coords = c("x","y"), crs = 4326)

dat_sf <- st_transform(dat_sf, 32742)

coords <- st_coordinates(dat_sf)

dat_pres <- dat_clean %>%
  mutate(pres = as.integer(CPUE > 0))

dat_pos <- dat_pres %>% filter(CPUE > 0, prof_bathy >= 500)

# Calcul du 95e percentile de la CPUE
cpue_95 <- quantile(dat_pos$CPUE, probs = 0.95, na.rm = TRUE)

# Filtrer les données pour ne garder que celles <= 95e percentile
dat_pos <- dat_pos %>%
  filter(CPUE <= cpue_95)


Study_domain_Hole <- readRDS("Data/Ker_polygon.rds")%>%
  mutate(geometry = st_geometry(.) / 1000)

# Update CRS to indicate units in kilometers (custom CRS)

st_crs(Study_domain_Hole) <- "+proj=utm +zone=42 +datum=WGS84 +units=km +no_defs"


mesh_pos <- make_mesh(dat_pos, xy_cols = c("X", "Y"), n_knots = 5)


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
    xlim = range(dat_pos$X) + c(-30, 30),
    ylim = range(dat_pos$Y) + c(-30, 30)
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
    filename = file.path("data/", output_plot, "mesh_plot2.png"),
    plot     = mesh_plot,
    width    = 8,
    height   = 6,
    dpi      = 300
  )
}


cat("→ mesh.png sauvegardé\n\n")


m_02 <- sdmTMB(
  CPUE ~ s(Year),
  data = dat_pos,
  mesh = bspde_pos,
  family = Gamma(link = "log"),
  spatial = "on",
  time= "Year",
  spatiotemporal = "iid"
)







# 4. Mesh -----------------------------------------------------------------

# Création du maillage spatial pour le modèle (choix du nombre de nœuds)
noeuds <- 120


# Supprime les lignes avec CPUE manquante
dat_clean <- dat %>% filter(!is.na(CPUE), !is.na(prof_bathy), !is.na(Year))

dat_sf <- st_as_sf(dat_clean, coords = c("x","y"), crs = 4326)

dat_sf <- st_transform(dat_sf, 32742)

coords <- st_coordinates(dat_sf)

dat_pres <- dat_clean %>%
  mutate(pres = as.integer(CPUE > 0))

dat_pos <- dat_pres %>% filter(CPUE > 0, prof_bathy >= 500)

# Calcul du 95e percentile de la CPUE
cpue_95 <- quantile(dat_pos$CPUE, probs = 0.95, na.rm = TRUE)

# Filtrer les données pour ne garder que celles <= 95e percentile
dat_pos <- dat_pos %>%
  filter(CPUE <= cpue_95)


Study_domain_Hole <- readRDS("Data/Ker_polygon.rds")%>%
  mutate(geometry = st_geometry(.) / 1000)

# Update CRS to indicate units in kilometers (custom CRS)

st_crs(Study_domain_Hole) <- "+proj=utm +zone=42 +datum=WGS84 +units=km +no_defs"


mesh_pos <- make_mesh(dat_pos, xy_cols = c("X", "Y"), n_knots = 30)


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
    xlim = range(dat_pos$X) + c(-30, 30),
    ylim = range(dat_pos$Y) + c(-30, 30)
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
    filename = file.path("data/", output_plot, "mesh_plot3.png"),
    plot     = mesh_plot,
    width    = 8,
    height   = 6,
    dpi      = 300
  )
}


cat("→ mesh.png sauvegardé\n\n")


m_03 <- sdmTMB(
  CPUE ~ s(Year),
  data = dat_pos,
  mesh = bspde_pos,
  family = Gamma(link = "log"),
  spatial = "on",
  time= "Year",
  spatiotemporal = "iid"
)





# 4. Mesh -----------------------------------------------------------------

# Création du maillage spatial pour le modèle (choix du nombre de nœuds)
noeuds <- 120


# Supprime les lignes avec CPUE manquante
dat_clean <- dat %>% filter(!is.na(CPUE), !is.na(prof_bathy), !is.na(Year))

dat_sf <- st_as_sf(dat_clean, coords = c("x","y"), crs = 4326)

dat_sf <- st_transform(dat_sf, 32742)

coords <- st_coordinates(dat_sf)

dat_pres <- dat_clean %>%
  mutate(pres = as.integer(CPUE > 0))

dat_pos <- dat_pres %>% filter(CPUE > 0, prof_bathy >= 500)

# Calcul du 95e percentile de la CPUE
cpue_95 <- quantile(dat_pos$CPUE, probs = 0.95, na.rm = TRUE)

# Filtrer les données pour ne garder que celles <= 95e percentile
dat_pos <- dat_pos %>%
  filter(CPUE <= cpue_95)


Study_domain_Hole <- readRDS("Data/Ker_polygon.rds")%>%
  mutate(geometry = st_geometry(.) / 1000)

# Update CRS to indicate units in kilometers (custom CRS)

st_crs(Study_domain_Hole) <- "+proj=utm +zone=42 +datum=WGS84 +units=km +no_defs"


mesh_pos <- make_mesh(dat_pos, xy_cols = c("X", "Y"), n_knots = 2000)


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
    xlim = range(dat_pos$X) + c(-30, 30),
    ylim = range(dat_pos$Y) + c(-30, 30)
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
    filename = file.path("data/", output_plot, "mesh_plot4.png"),
    plot     = mesh_plot,
    width    = 8,
    height   = 6,
    dpi      = 300
  )
}


cat("→ mesh.png sauvegardé\n\n")


m_04 <- sdmTMB(
  CPUE ~ s(Year),
  data = dat_pos,
  mesh = bspde_pos,
  family = Gamma(link = "log"),
  spatial = "on",
  time= "Year",
  spatiotemporal = "iid"
)


sdmTMB::tidy(m_0, "ran_pars")
sdmTMB::tidy(m_02, "ran_pars")
sdmTMB::tidy(m_03, "ran_pars")
sdmTMB::tidy(m_04, "ran_pars")

AIC(m_0, m_02, m_03, m_04)
