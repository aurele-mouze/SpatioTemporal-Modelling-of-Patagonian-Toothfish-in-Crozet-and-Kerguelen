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

dat_map <- dat_pos %>%
  mutate(
    PRES_CACHALOT = case_when(
      PRES_CACHALOT == "présence" ~ 1,
      PRES_CACHALOT == "absence" ~ 0,
      PRES_CACHALOT == "non obsé" ~ -1
    )
  )


dat_Orque <- dat_pos %>%
  mutate(
    PRES_CACHALOT = case_when(
      PRES_ORQUE == "présence" ~ 1,
      PRES_ORQUE == "absence" ~ 0,
      PRES_ORQUE == "non obsé" ~ -1
    )
  )

dat_pos <- st_transform(dat_pos, 32742)

dat_map <- st_transform(dat_map, 32742)

dat_Orque <- st_transform(dat_Orque, 32742)

# 4. Faire une carte avec les CPUE moyennes sur des carrées de (définir la résolution) ------------------------------

#D'abord en moyenne sans prendre en compte l'année
#Fais une grille de carré

grille_10km <- st_make_grid(dat_pos,
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





mean_obs <- ggplot() + xlim(63.6, 73) + ylim(-52, -46) + xlab("Longitude") + ylab("Latitude") +
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
  geom_sf(data = grille_cpue, aes(fill = CPUE_moy), color = NA, alpha = 0.7) +
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


# 5. Visualitsation du nombre d'hameçon

ggplot(dat_pos, aes(x = nb_ham)) +
  geom_bar(fill = "steelblue") +
  labs(
    title = "Fréquence des nombres d'hameçons",
    x = "Nombre d'hameçons",
    y = "Fréquence"
  ) +
  theme_minimal()



histo_ham <- ggplot(dat_pos, aes(x = nb_ham)) +
  geom_histogram(binwidth = 1000, fill = "steelblue", color = "black") +
  labs(
    title = "Histogramme du nombre d'hameçons",
    x = "Nombre d'hameçons",
    y = "Nombre de cellules"
  ) +
  theme_minimal()

if (SavePlots) {
  ggsave(
    filename = file.path("data/", output_plot, "Histo_Ham.png"),
    plot     = histo_ham,
    width    = 8,
    height   = 6,
    dpi      = 300
  )
}



dat_pos_ham <- dat_pos %>%
  mutate(
    nb_ham_cat = cut(
      nb_ham,
      breaks = seq(min(nb_ham), max(nb_ham), by = 1000),  # intervalle de 5 hameçons
      include.lowest = TRUE,
      right = FALSE,
      labels = paste(head(seq(min(nb_ham), max(nb_ham)-1000, by=1000), -1),
                     tail(seq(min(nb_ham)+5, max(nb_ham), by=1000), -1),
                     sep = "-")
    )
  )

CPUE_ham <- ggplot(dat_pos_ham, aes(x = nb_ham_cat, y = CPUE)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "red") + 
  labs(
    title = "Distribution de la CPUE selon le nombre d'hameçons",
    x = "Intervalle du nombre d'hameçons",
    y = "CPUE"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


if (SavePlots) {
  ggsave(
    filename = file.path("data/", output_plot, "CPUE_Ham.png"),
    plot     = CPUE_ham,
    width    = 8,
    height   = 6,
    dpi      = 300
  )
}


ham_categories <- sort(unique(dat_pos_ham$nb_ham_cat))

for(cat in ham_categories){
  
  data_cat <- dat_pos_ham %>% filter(nb_ham_cat == cat)
  
  points_sf <- st_as_sf(data_cat, coords = c("X", "Y"), crs = 32742)
  points_sf <- st_transform(points_sf, 4326)
  
  carte <- ggplot() +
    
    geom_sf(data = Isobath_KER_sf %>% filter(level == "-500"), colour = "grey80") +
    geom_sf(data = Isobath_KER_sf %>% filter(level == "-1000"), colour = "grey60") +
    geom_sf(data = EEZ_sf, fill = "transparent") +
    geom_sf(data = FR_island_sf, fill = "darkgrey") +
    
    geom_sf(data = points_sf, color = "steelblue", size = 2) +
    
    coord_sf(
      xlim = c(63.6, 73),
      ylim = c(-52, -46),
      expand = FALSE
    ) +
    
    theme_minimal()
  
  if (SavePlots) {
    ggsave(
      filename = file.path("data/", output_plot, paste0("map_Hameçon_", cat, ".png")),
      plot = carte,
      width = 8,
      height = 6,
      dpi = 300
    )}
}




# 6. Temps de trempe

ggplot(dat_pos, aes(x = tps_trempe_hm_num)) +
  geom_bar(fill = "steelblue") +
  labs(
    title = "Fréquence des nombres d'hameçons",
    x = "Nombre d'hameçons",
    y = "Fréquence"
  ) +
  theme_minimal()



histo_ham <- ggplot(dat_pos, aes(x = tps_trempe_hm_num)) +
  geom_histogram(binwidth = 2, fill = "steelblue", color = "black") +
  labs(
    title = "Histogramme dTemps de trempe",
    x = "Temps de trempe",
    y = "Nombre de cellules"
  ) +
  theme_minimal()

if (SavePlots) {
  ggsave(
    filename = file.path("data/", output_plot, "Histo_trempe.png"),
    plot     = histo_ham,
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







dat_join <- st_join(dat_pos, grille_10km)

ham_grille <- dat_join %>%
  group_by(id) %>%
  summarise(
    nb_ham_total = sum(nb_ham, na.rm = TRUE),
    n_station = n()
  )

grille_ham <- left_join(
  grille_10km,
  st_drop_geometry(ham_grille),
  by = "id"
)

grille_ham_ll <- st_transform(grille_ham, 4326)



mean_obs <- ggplot() + xlim(63.6, 73) + ylim(-52, -46) + xlab("Longitude") + ylab("Latitude") +
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
  geom_sf(data = grille_ham_ll, aes(fill = nb_ham_total), color = NA) +
  # geom_sf(data = dat_sf, size = 1, color = "black") +
  scale_fill_viridis_c(
    trans = "log10",
    na.value = "transparent",
    name = "Nb hameçons (log)"
  ) +
  theme_bw()

if (SavePlots) {
  ggsave(
    filename = file.path("data/", output_plot, "sum_ham_obs.png"),
    plot     = mean_obs,
    width    = 8,
    height   = 6,
    dpi      = 300
  )
}



#Idem mais une carte par année
annees <- sort(unique(dat_pos$Year))

liste_Ham <- list()

for (an in annees) {
  
  dat_an <- dat_pos %>%
    filter(Year == an)
  
  dat_join <- st_join(dat_an, grille_10km)
  
  ham_grille <- dat_join %>%
    group_by(id) %>%
    summarise(
      nb_ham_total = sum(nb_ham, na.rm = TRUE),
      n_station = n()
    )
  
  grille_ham <- left_join(
    grille_10km,
    st_drop_geometry(ham_grille),
    by = "id"
  )
  
  grille_ham_ll <- st_transform(grille_ham, 4326)
  
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
      data = grille_ham_ll,
      aes(fill = nb_ham_total),
      color = NA,
      alpha = 0.7
    ) +
    scale_fill_viridis_c(
      na.value = "transparent") +
    ggtitle(paste("Année :", an))
  
  liste_Ham[[as.character(an)]] <- carte
}

for (an in annees) {
  ggsave(
    filename = file.path("data/", output_plot, paste0("Ham_", an, ".png")),
    plot     = liste_Ham[[as.character(an)]],
    width    = 8,
    height   = 6,
    dpi      = 300
  )
}






dat_pres <- dat_map %>%
  filter(PRES_CACHALOT == 1)


dat_join_pred <- st_join(dat_pres, grille_10km)

Cach_grille <- dat_join_pred %>%
  group_by(id) %>%
  summarise(
    NB_OBS = n()   # 👉 nombre d'observations par cellule
  )

Cach_grille <- left_join(
  grille_10km,
  st_drop_geometry(Cach_grille),
  by = "id"
)

Cach_grille_ll <- st_transform(Cach_grille, 4326)

mean_obs <- ggplot() + xlim(63.6, 73) + ylim(-52, -46) + xlab("Longitude") + ylab("Latitude") +
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
  geom_sf(data = Cach_grille_ll, aes(fill = NB_OBS), color = NA) +
  # geom_sf(data = dat_sf, size = 1, color = "black") +
  scale_fill_viridis_c(
    trans = "log10",
    na.value = "transparent",
    name = "Nb de Cachalot Observé"
  ) +
  theme_bw()

if (SavePlots) {
  ggsave(
    filename = file.path("data/", output_plot, "sum_Cach_obs.png"),
    plot     = mean_obs,
    width    = 8,
    height   = 6,
    dpi      = 300
  )
}


df_count <- dat_map %>%
  mutate(
    pres = case_when(
      PRES_CACHALOT == 1 ~ "Présence",
      PRES_CACHALOT == 0 ~ "Absence",
      PRES_CACHALOT == -1 ~ "Non observé",
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(Year, pres) %>%
  summarise(n = n(), .groups = "drop")

ggplot(df_count, aes(x = Year, y = n, fill = pres)) +
  geom_col() +
  scale_fill_manual(values = c(
    "Présence" = "red",
    "Absence" = "blue",
    "Non observé" = "grey70"
  )) +
  labs(
    x = "Year",
    y = "Number of data",
    fill = "Status"
  ) +
  theme_bw()



dat_pres <- dat_pres %>%
  mutate(
    NB_MEAN_CACHALOT = (NB_MIN_CACHALOT + NB_MAX_CACHALOT) / 2
  )

dat_sf <- st_as_sf(dat_pres, coords = c("x", "y"), crs = 4326)
dat_sf <- st_transform(dat_sf, 32742)  # même UTM que ton modèle
coords <- st_coordinates(dat_sf)

x <- coords[,1]
y <- coords[,2]

w <- dat_sf$NB_MEAN_CACHALOT
library(MASS)

# arrondir les poids pour éviter explosion mémoire
w_round <- round(w)
w_clean <- ifelse(is.na(w), 0, w)   # remplace NA par 0
w_clean <- pmax(w_clean, 0)         # évite valeurs négatives
w_round <- round(w_clean)
# réplication des points
x_rep <- rep(x, w_round)
y_rep <- rep(y, w_round)

kde <- kde2d(
  x_rep,
  y_rep,
  n = 200
)

kde_df <- expand.grid(
  X = kde$x,
  Y = kde$y
)

kde_df$Z <- as.vector(kde$z)



ggplot(kde_df, aes(x = X/1000, y = Y/1000, fill = Z)) +
  geom_raster() +
  scale_fill_viridis_c(name = "Intensity") +
  coord_equal() +
  theme_bw()


FR_island_utm <- st_transform(FR_island_sf, 32742)
EEZ_utm       <- st_transform(EEZ_sf, 32742)
Isobath_utm   <- st_transform(Isobath_KER_sf, 32742)
kde_df$Z[kde_df$Z < 1e-12] <- NA
mean_obs <- ggplot() +
  xlim(range(kde_df$X))+
  ylim(range(kde_df$Y)) +
  geom_sf(data = Isobath_utm %>% filter(level == "-500"), color = "grey80") +
  geom_sf(data = EEZ_utm, fill = NA) +
  geom_sf(data = FR_island_utm, fill = "darkgrey") +
  geom_raster(data = kde_df, aes(x = X, y = Y, fill = Z)) +
  coord_sf() +
  scale_fill_viridis_c(
    na.value = "transparent"
  )+
  theme_bw()










dat_or <- dat_Orque %>%
  filter(PRES_ORQUE == "présence")


dat_join_pred <- st_join(dat_or, grille_10km)

Or_grille <- dat_join_pred %>%
  group_by(id) %>%
  summarise(
    NB_OBS = n()   # 👉 nombre d'observations par cellule
  )

Or_grille <- left_join(
  grille_10km,
  st_drop_geometry(Or_grille),
  by = "id"
)

Or_grille_ll <- st_transform(Or_grille, 4326)

mean_obs <- ggplot() + xlim(63.6, 73) + ylim(-52, -46) + xlab("Longitude") + ylab("Latitude") +
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
  geom_sf(data = Or_grille_ll, aes(fill = NB_OBS), color = NA) +
  # geom_sf(data = dat_sf, size = 1, color = "black") +
  scale_fill_viridis_c(
    trans = "log10",
    na.value = "transparent",
    name = "Number of observation of Killer whales (in log)"
  ) +
  theme_bw()

if (SavePlots) {
  ggsave(
    filename = file.path("data/", output_plot, "sum_Orq_obs.png"),
    plot     = mean_obs,
    width    = 8,
    height   = 6,
    dpi      = 300
  )
}




ggplot(dat_Orque, aes(x = Year, fill = PRES_ORQUE)) +
  geom_bar() +
  scale_fill_manual(values = c(
    "présence" = "red",
    "absence" = "blue",
    "non obsé" = "grey70"
  )) +
  labs(
    x = "Year",
    y = "Number of data",
    fill = "Status"
  ) +
  theme_bw()

