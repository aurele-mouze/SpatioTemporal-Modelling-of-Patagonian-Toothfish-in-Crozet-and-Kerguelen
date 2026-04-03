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

dat_pos <- dat_pres %>%
  filter(CPUE > 0, prof_bathy >= 500) %>%
  st_as_sf(coords = c("x", "y"), crs = 4326)

dat_pos <- st_transform(dat_pos, 32742)

# Calcul de la CPUE moyenne par station ou unité si besoin
dat_summary <- dat_pos %>%
  group_by(StationID) %>%   # ou juste summariser sur tout le jeu de données
  summarise(CPUE_mean = mean(CPUE, na.rm = TRUE)) 

# Boxplot avec intervalle de confiance 95%
ggplot(dat_summary, aes(x = "", y = CPUE_mean)) +
  geom_boxplot(outlier.color = "red", outlier.shape = 16, width = 0.3) +
  labs(
    x = "",
    y = "CPUE moyenne",
    title = "Boxplot de la CPUE moyenne avec IC 95%"
  ) +
  theme_minimal()

# 4. Faire une carte avec les CPUE moyennes sur des carrées de (définir la résolution) ------------------------------

#D'abord en moyenne sans prendre en compte l'année
#Fais une grille de carré



grille_10km <- st_make_grid(dat_sf,
                            cellsize = 10000,
                            square = TRUE)

grille_10km <- st_sf(
  id = 1:length(grille_10km),
  geometry = grille_10km
)


dat_join <- st_join(dat_pos, grille_10km)

grille_cpue <- dat_join %>%
  group_by(id) %>%
  summarise(
    CPUE_moy = mean(CPUE, na.rm = TRUE),
    n_station = n()
  )

#grille_cpue <- left_join(grille_10km, st_drop_geometry(cpue_grille), by = "id")





mean_obs <- ggplot() + xlim(48.5, 53.5) + ylim(-47.5, -45) + xlab("Longitude") + ylab("Latitude") +
  # geom_sf(
  #   data = RNN_AMP_sf %>% dplyr::filter(protection %in% c("R")),
  #   colour = "#b2182b",
  #   inherit.aes = F,
  #   linewidth = 0.5,
  #   # fill = "white"   
  #   fill = "transparent"
  #   ) +
  geom_sf(
    data = Isobath_KER_sf  %>% dplyr::filter(level %in% c("-500")),
    colour = "grey80",
    inherit.aes = F,
    size = 0.5
  ) +
  geom_sf(
    data = Isobath_KER_sf  %>% dplyr::filter(level %in% c("-1000")),
    colour = "grey60",
    inherit.aes = F,
    size = 0.5
  ) +
  # geom_sf(
  #   data = Isobath_KER_sf  %>% dplyr::filter(level %in% c("-1500")),
  #   colour = "grey30",
  #   inherit.aes = F,
  #   size = 0.5
  # ) +
  # geom_sf(
  #   data = Isobath_KER_sf  %>% dplyr::filter(level %in% c("-2000")),
  #   colour = "black",
  #   inherit.aes = F,
  #   size = 0.5
  # ) +
  geom_sf(
    data = EEZ_sf,
    fill = "transparent",
    inherit.aes = F,
    size = 0.5
  ) +
  geom_sf(data = FR_island_sf ,
          fill = "darkgrey",
          inherit.aes = F
  ) + 
  geom_sf(data = grille_cpue, aes(fill = CPUE_moy), color = NA) +
  # geom_sf(data = dat_sf, size = 1, color = "black") +
  scale_fill_viridis_c(na.value = "transparent") +
  theme_bw()

if (SavePlots) {
  ggsave(
    filename = file.path("data/", output_plot, "mean_obs.png"),
    plot     = mean_obs,
    width    = 8,
    height   = 6,
    dpi      = 300
  )
}


#Idem mais une carte par année
annees <- sort(unique(dat_pos$Year))

liste_cartes <- list()

cpue_moyenne_global <- dat_pos %>%
  st_join(grille_10km) %>%
  group_by(id) %>%
  summarise(CPUE_moy = mean(CPUE, na.rm = TRUE))


limites <- quantile(cpue_moyenne_global$CPUE_moy, probs = c(0, 0.95), na.rm = TRUE)

for (an in annees) {
  
  dat_an <- dat_pos %>%
    filter(Year == an)
  
  dat_join <- st_join(dat_an, grille_10km)
  
  cpue_grille <- dat_join %>%
    group_by(id) %>%
    summarise(
      CPUE_moy = mean(CPUE, na.rm = TRUE),
      n_station = n()
    )
  
  grille_cpue <- left_join(
    grille_10km,
    st_drop_geometry(cpue_grille),
    by = "id"
  )
  
  grille_cpue_ll <- st_transform(grille_cpue, 4326)
  
  carte <- ggplot() +
    xlim(48.5, 53.5) + ylim(-47.5, -45) +
    xlab("Longitude") + ylab("Latitude") +
    
    geom_sf(
      data = Isobath_KER_sf %>% filter(level %in% c("-500")),
      colour = "grey80",
      inherit.aes = FALSE,
      size = 0.5
    ) +
    geom_sf(
      data = Isobath_KER_sf %>% filter(level %in% c("-1000")),
      colour = "grey60",
      inherit.aes = FALSE,
      size = 0.5
    ) +
    geom_sf(
      data = EEZ_sf,
      fill = "transparent",
      inherit.aes = FALSE,
      size = 0.5
    ) +
    geom_sf(
      data = FR_island_sf,
      fill = "darkgrey",
      inherit.aes = FALSE
    ) +
    geom_sf(
      data = grille_cpue_ll,
      aes(fill = CPUE_moy),
      color = NA,
      alpha = 0.7
    ) +
    scale_fill_viridis_c(
      limits = limites,
      na.value = "transparent",
      name = "CPUE"
    ) +
    ggtitle(paste("Année :", an))
  
  liste_cartes[[as.character(an)]] <- carte
}

for (an in annees) {
  ggsave(
    filename = file.path("data/", output_plot, paste0("CPUE_", an, ".png")),
    plot     = liste_cartes[[as.character(an)]],
    width    = 8,
    height   = 6,
    dpi      = 300
  )
}





# Mois

month <- ggplot(dat_pos, aes(x = factor(Month))) +
  geom_bar(fill = "steelblue") +
  labs(
    x = "Mois",
    y = "Nombre d'observations",
    title = "Nombre d'observations par mois"
  ) +
  theme_minimal()

if (SavePlots) {
  ggsave(
    filename = file.path("data/", output_plot, "obs_month.png"),
    plot     = month,
    width    = 8,
    height   = 6,
    dpi      = 300
  )
}




mois <- sort(unique(dat_pos$Month))

liste_cartes <- list()

cpue_moyenne_global <- dat_pos %>%
  st_join(grille_10km) %>%
  group_by(id) %>%
  summarise(CPUE_moy = mean(CPUE, na.rm = TRUE))


limites <- quantile(cpue_moyenne_global$CPUE_moy, probs = c(0, 0.95), na.rm = TRUE)

for (an in mois) {
  
  dat_an <- dat_pos %>%
    filter(Month == an)
  
  dat_join <- st_join(dat_an, grille_10km)
  
  cpue_grille <- dat_join %>%
    group_by(id) %>%
    summarise(
      CPUE_moy = mean(CPUE, na.rm = TRUE),
      n_station = n()
    )
  
  grille_cpue <- left_join(
    grille_10km,
    st_drop_geometry(cpue_grille),
    by = "id"
  )
  
  grille_cpue_ll <- st_transform(grille_cpue, 4326)
  
  carte <- ggplot() +
    xlim(63.6, 73) + ylim(-52, -46) +
    xlab("Longitude") + ylab("Latitude") +
    
    geom_sf(
      data = Isobath_KER_sf %>% filter(level %in% c("-500")),
      colour = "grey80",
      inherit.aes = FALSE,
      size = 0.5
    ) +
    geom_sf(
      data = Isobath_KER_sf %>% filter(level %in% c("-1000")),
      colour = "grey60",
      inherit.aes = FALSE,
      size = 0.5
    ) +
    geom_sf(
      data = EEZ_sf,
      fill = "transparent",
      inherit.aes = FALSE,
      size = 0.5
    ) +
    geom_sf(
      data = FR_island_sf,
      fill = "darkgrey",
      inherit.aes = FALSE
    ) +
    geom_sf(
      data = grille_cpue_ll,
      aes(fill = CPUE_moy),
      color = NA,
      alpha = 0.7
    ) +
    scale_fill_viridis_c(
      limits = limites,
      na.value = "transparent",
      name = "CPUE"
    ) +
    ggtitle(paste("Mois :", an))
  
  liste_cartes[[as.character(an)]] <- carte
}

for (an in mois) {
  ggsave(
    filename = file.path("data/", output_plot, paste0("CPUE_Mois_", an, ".png")),
    plot     = liste_cartes[[as.character(an)]],
    width    = 8,
    height   = 6,
    dpi      = 300
  )
}