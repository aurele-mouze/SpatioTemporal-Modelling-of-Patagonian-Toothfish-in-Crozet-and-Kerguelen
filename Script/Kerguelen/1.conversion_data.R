### load les datas ---------------------------------------------------------------

# Chargement des jeux de données contenant les captures et les données biologiques
biometrie <- readRDS("Data/biometrie.rds")
capture <- readRDS("Data/capture_Agg.rds")
Station <- readRDS("Data/Station.rds")

### Fusionner les datas pour avoir 1 seul dataset----------------------------------

capture <- Station %>%
  left_join(
    capture %>%
      ungroup() %>%
      dplyr::select(-c(Zone, GearType, TAAFseason)),
    by = "StationID"
  )

biometrie <- biometrie %>%
  filter(
    Type == "Bio",
  )


biometrie <- biometrie %>%
  left_join(
    capture %>%
      dplyr::select(StationID, HaulEndDate),
    by = "StationID"
  )

### Trier dans un premier temps sur les données de Kerguelent pour une année (exemple : 2015)------

data_all_cap <- capture %>%
  filter(
    Plateau == "kerguelen",
    # between(HaulStartLon, 66, 67.5),
    # between(HaulStartLat,  -50, -48.5)
    
  )

data_all_bio <- biometrie


# Calcul du temps de trempe -----------------------------------------------------

# Extraction des colonnes de temps de début et fin
tps_trempe = data_all_cap %>%
  dplyr::select(StationID, HaulStartDate, HaulEndDate) %>%
  mutate(
    # Conversion des dates en format POSIXct
    datetime_fil = as.POSIXct(HaulStartDate, format = "%d/%m/%y %H:%M:%S"),
    datetime_vir = as.POSIXct(HaulEndDate, format = "%d/%m/%y %H:%M:%S"),
    # Calcul de la durée en minutes
    duration_mins = as.numeric(difftime(datetime_vir, datetime_fil, units = "mins")),
    # Conversion en heures et minutes
    heures = floor(duration_mins / 60),
    minutes = round(duration_mins %% 60),
    # Format final "hh:mm"
    tps_trempe_hm = sprintf("%02d:%02d", heures, minutes)
  )

# Garde uniquement le temps de trempe final
tps_trempe = tps_trempe %>%
  dplyr::select(StationID, tps_trempe_hm)

# Ajout du temps de trempe dans les données de capture
data_all_cap = left_join(data_all_cap, tps_trempe, by = c("StationID"))

# Calcul de la distance réelle des palangres ------------------------------------

# Récupère les coordonnées de début et de fin
coord_fil <- data_all_cap %>%
  dplyr::select(HaulStartLon, HaulStartLat) %>%
  rename(x = HaulStartLon, y = HaulStartLat)

coord_vir <- data_all_cap %>%
  dplyr::select(HaulEndLon, HaulEndLat) %>%
  rename(x = HaulEndLon, y = HaulEndLat)

distance <- as.data.frame(
  terra::distance(
    coord_fil,
    coord_vir,
    pairwise = TRUE,
    lonlat = TRUE
  )
)

colnames(distance) <- "distance"


# Ajout des distances aux données de capture
data_all_cap <- cbind(data_all_cap, distance)

# Calcul de la différence des classes de taille ---------------------------------

# Calcule l'écart interquantile (Q95 - Q05) des longueurs totales par station
diff_taille = data_all_bio %>%
  dplyr::select(StationID, Length) %>%
  group_by(StationID) %>%
  summarise(
    q95 = quantile(Length, probs = 0.95, na.rm = TRUE),
    q05 = quantile(Length, probs = 0.05, na.rm = TRUE),
    diff_taille = q95 - q05,
    .groups = "drop"
  ) %>%
  dplyr::select(StationID, diff_taille)

# Ajout des écarts de taille à data_all_cap
data_all_cap = left_join(data_all_cap, diff_taille, by = c("StationID"))


# Calcul de la variation de profondeur ------------------------------------------

# Chargement du raster bathymétrique
map = rast(here(path_rast_bathy, 'ras_bathy_pecheker_interpol.tif'))
plot(map)

# Initialisation de la liste pour stocker les lignes densifiées
data_all_cap_nonnull <- data_all_cap %>%
  filter(distance > 0)

lines_dens_list <- list()

# Boucle sur chaque ligne de data_all_cap pour créer des lignes densifiées
for (i in seq_len(nrow(data_all_cap_nonnull))) {
  start_coords <- c(data_all_cap_nonnull$HaulStartLon[i], data_all_cap_nonnull$HaulStartLat[i])
  end_coords   <- c(data_all_cap_nonnull$HaulEndLon[i], data_all_cap_nonnull$HaulEndLat[i])
  
  segment <- rbind(start_coords, end_coords)
  line <- vect(segment, type = "lines", crs = crs(map))
  
  spacing <- data_all_cap_nonnull$distance[i] 
  line_densified <- densify(line, spacing)
  line_densified$StationID <- data_all_cap_nonnull$StationID[i]
  
  lines_dens_list[[i]] <- line_densified
}

# Fusion de toutes les lignes densifiées
lines_densified_all <- do.call(rbind, lines_dens_list)

# Extraction des points des lignes densifiées
points_df <- geom(lines_densified_all, df = TRUE)[, c("x", "y", "geom")]

points_df <- points_df %>%
  dplyr::mutate(
    StationID = lines_densified_all$StationID[geom],
    x_fin = x,
    y_fin = y
  ) %>%
  dplyr::select(StationID, x_fin, y_fin) %>%
  dplyr::group_by(StationID) %>%
  dplyr::mutate(ID = row_number()) %>%
  dplyr::arrange(StationID, ID) %>%
  dplyr::ungroup()

# Calcul des points médians entre paires successives
points_mid <- points_df %>%
  group_by(StationID) %>%
  mutate(x_next = lead(x_fin), y_next = lead(y_fin)) %>%
  filter(!is.na(x_next) & !is.na(y_next)) %>%
  mutate(x_mid = (x_fin + x_next) / 2, y_mid = (y_fin + y_next) / 2) %>%
  dplyr::select(StationID, ID, x_mid, y_mid)

# Extraction des profondeurs sur les points médians
coords_mid <- cbind(points_mid$x_mid, points_mid$y_mid)
points_mid$valeur_raster <- terra::extract(map[["profondeur.pred"]], coords_mid)[[1]]

# Calcul des profondeurs min et max pour chaque station
depth_by_station <- points_mid %>%
  filter(!is.na(valeur_raster)) %>%
  group_by(StationID) %>%
  summarise(
    prof_min = min(valeur_raster),
    prof_max = max(valeur_raster),
    diff_profondeur = prof_max - prof_min,
    .groups = "drop"
  )

depth_by_station <- depth_by_station %>%
  mutate(StationID = as.character(StationID))

data_all_cap <- data_all_cap_nonnull %>%
  mutate(StationID = as.character(StationID)) %>%
  left_join(depth_by_station, by = "StationID")



### Extraction des valeurs de pente ---------------------------------------------

# Copier points_mid dans points_mid_slope
points_mid_slope <- points_mid

# Charger le raster de pente
slope_map <- rast(here(path_rast_bathy, 'ras_slope_pecheker_interpol.tif'))

# Extraction des valeurs de pente
coords_mid <- cbind(points_mid_slope$x_mid, points_mid_slope$y_mid)
points_mid_slope$valeur_raster <- terra::extract(slope_map[["slope"]], coords_mid)[[1]]

# Préparer le tableau des pentes par station et par point (avec ID)
slope_by_station <- points_mid_slope %>%
  dplyr::select(StationID, ID, valeur_raster) %>%
  dplyr::rename(slope = valeur_raster)







### Calul des flags -------------------------------------------------------------

# Nettoyage des valeurs infinies
data_all_cap <- data_all_cap %>%
  mutate(across(
    c(tps_trempe_hm, distance, diff_taille, diff_profondeur),
    ~ ifelse(is.infinite(.), NA, .)
  ))

# Conversion du temps de trempe "hh:mm" en nombre d'heures
data_all_cap <- data_all_cap %>%
  mutate(
    tps_trempe_hm_num = as.numeric(lubridate::hm(tps_trempe_hm)) / 3600
  )

# Calcul des quantiles pour les variables de contrôle qualité
q_vals <- data_all_cap %>%
  summarise(
    q95_tps = quantile(tps_trempe_hm_num, probs = 0.95, na.rm = TRUE),
    q05_tps = quantile(tps_trempe_hm_num, probs = 0.05, na.rm = TRUE),
    q95_dist = quantile(distance, probs = 0.95, na.rm = TRUE),
    q05_dist = quantile(distance, probs = 0.05, na.rm = TRUE),
    q95_taille = quantile(diff_taille, probs = 0.95, na.rm = TRUE),
    q05_taille = quantile(diff_taille, probs = 0.05, na.rm = TRUE),
    q95_prof = quantile(diff_profondeur, probs = 0.95, na.rm = TRUE),
    q05_prof = quantile(diff_profondeur, probs = 0.05, na.rm = TRUE)
  )

# Création du flag binaire (1 si la condition est dans l’intervalle, 0 sinon)
data_flag <- data_all_cap %>%
  mutate(
    flag =
      replace_na(as.integer(distance >= q_vals$q05_dist & distance <= q_vals$q95_dist), 0) +
      replace_na(as.integer(diff_taille >= q_vals$q05_taille & diff_taille <= q_vals$q95_taille), 0) +
      replace_na(as.integer(diff_profondeur <= q_vals$q95_prof), 0)
  )

# Création d’un tableau de détail des conditions de flag par variable
df_flag_detail <- data_all_cap %>%
  mutate(
    cond_tps    = replace_na(as.integer(tps_trempe_hm_num >= q_vals$q05_tps & tps_trempe_hm_num <= q_vals$q95_tps), 0),
    cond_dist   = replace_na(as.integer(distance         >= q_vals$q05_dist & distance         <= q_vals$q95_dist), 0),
    cond_taille = replace_na(as.integer(diff_taille      >= q_vals$q05_taille & diff_taille    <= q_vals$q95_taille), 0),
    cond_prof   = replace_na(as.integer(diff_profondeur  <= q_vals$q95_prof), 0),
    flag_total  = cond_tps + cond_dist + cond_taille + cond_prof
  ) %>%
  dplyr::select(StationID, cond_tps, cond_dist, cond_taille, cond_prof, flag_total)

# Choix du jeu de données pour le modèle ----------------------------------------

# Extraction des flags
data_modelo = data_flag %>%
  dplyr::select(StationID, flag)


points_mid <- points_mid %>%
  mutate(StationID = as.character(StationID))


# Ajout des coordonnées et profondeur bathy
data_modelo = left_join(data_modelo, points_mid, by = c("StationID"))

# Ajout des données de captures (CPUE)
data_CPUE <- data_all_cap %>%
  dplyr::select(
    StationID,
    ScanNumber,
    ScanWeightKg,
    HookHaul,
    TargetSpecies,
    HaulEndDate
  )

data_modelo = left_join(data_modelo, data_CPUE, by = c("StationID"))


data_modelo <- data_modelo %>%
  mutate(Year = year(HaulEndDate))

# Renommage et transformation finale
data_modelo = data_modelo %>%
  rename(x = x_mid,
         y = y_mid,
         prof_bathy = valeur_raster) %>%
  mutate()

# Ajout des pentes
slope_by_station <- slope_by_station %>%
  mutate(StationID = as.character(StationID))

data_modelo <- left_join(data_modelo, slope_by_station, by = c("StationID", "ID"))


# Ajout du temps de trempe numérique
data_tps = data_all_cap %>%
  dplyr::select(StationID, tps_trempe_hm_num)

data_modelo <- left_join(data_modelo, data_tps, by = c("StationID"))

data_modelo = left_join(data_modelo, df_flag_detail, by = c("StationID"))

# Sauvegarde de l’objet final

saveRDS(data_modelo, file = here("data", output_df, "data_modelo.rds"))



















### Trier dans un premier temps sur les données de Kerguelent pour une année (exemple : 2015)------

data_all_cap <- capture %>%
  filter(
    Plateau == "crozet"
    
  )

data_all_bio <- biometrie


# Calcul du temps de trempe -----------------------------------------------------

# Extraction des colonnes de temps de début et fin
tps_trempe = data_all_cap %>%
  dplyr::select(StationID, HaulStartDate, HaulEndDate) %>%
  mutate(
    # Conversion des dates en format POSIXct
    datetime_fil = as.POSIXct(HaulStartDate, format = "%d/%m/%y %H:%M:%S"),
    datetime_vir = as.POSIXct(HaulEndDate, format = "%d/%m/%y %H:%M:%S"),
    # Calcul de la durée en minutes
    duration_mins = as.numeric(difftime(datetime_vir, datetime_fil, units = "mins")),
    # Conversion en heures et minutes
    heures = floor(duration_mins / 60),
    minutes = round(duration_mins %% 60),
    # Format final "hh:mm"
    tps_trempe_hm = sprintf("%02d:%02d", heures, minutes)
  )

# Garde uniquement le temps de trempe final
tps_trempe = tps_trempe %>%
  dplyr::select(StationID, tps_trempe_hm)

# Ajout du temps de trempe dans les données de capture
data_all_cap = left_join(data_all_cap, tps_trempe, by = c("StationID"))

# Calcul de la distance réelle des palangres ------------------------------------

# Récupère les coordonnées de début et de fin
coord_fil <- data_all_cap %>%
  dplyr::select(HaulStartLon, HaulStartLat) %>%
  rename(x = HaulStartLon, y = HaulStartLat)

coord_vir <- data_all_cap %>%
  dplyr::select(HaulEndLon, HaulEndLat) %>%
  rename(x = HaulEndLon, y = HaulEndLat)

distance <- as.data.frame(
  terra::distance(
    coord_fil,
    coord_vir,
    pairwise = TRUE,
    lonlat = TRUE
  )
)

colnames(distance) <- "distance"


# Ajout des distances aux données de capture
data_all_cap <- cbind(data_all_cap, distance)

# Calcul de la différence des classes de taille ---------------------------------

# Calcule l'écart interquantile (Q95 - Q05) des longueurs totales par station
diff_taille = data_all_bio %>%
  dplyr::select(StationID, Length) %>%
  group_by(StationID) %>%
  summarise(
    q95 = quantile(Length, probs = 0.95, na.rm = TRUE),
    q05 = quantile(Length, probs = 0.05, na.rm = TRUE),
    diff_taille = q95 - q05,
    .groups = "drop"
  ) %>%
  dplyr::select(StationID, diff_taille)

# Ajout des écarts de taille à data_all_cap
data_all_cap = left_join(data_all_cap, diff_taille, by = c("StationID"))


# Calcul de la variation de profondeur ------------------------------------------

# Chargement du raster bathymétrique
map = rast(here(path_rast_bathy, 'ras_bathy_pecheker_interpolCRO.tif'))
plot(map)

# Initialisation de la liste pour stocker les lignes densifiées
data_all_cap_nonnull <- data_all_cap %>%
  filter(distance > 0)

lines_dens_list <- list()

# Boucle sur chaque ligne de data_all_cap pour créer des lignes densifiées
for (i in seq_len(nrow(data_all_cap_nonnull))) {
  start_coords <- c(data_all_cap_nonnull$HaulStartLon[i], data_all_cap_nonnull$HaulStartLat[i])
  end_coords   <- c(data_all_cap_nonnull$HaulEndLon[i], data_all_cap_nonnull$HaulEndLat[i])
  
  segment <- rbind(start_coords, end_coords)
  line <- vect(segment, type = "lines", crs = crs(map))
  
  spacing <- data_all_cap_nonnull$distance[i]
  line_densified <- densify(line, spacing)
  line_densified$StationID <- data_all_cap_nonnull$StationID[i]
  
  lines_dens_list[[i]] <- line_densified
}

# Fusion de toutes les lignes densifiées
lines_densified_all <- do.call(rbind, lines_dens_list)

# Extraction des points des lignes densifiées
points_df <- geom(lines_densified_all, df = TRUE)[, c("x", "y", "geom")]

points_df <- points_df %>%
  dplyr::mutate(
    StationID = lines_densified_all$StationID[geom],
    x_fin = x,
    y_fin = y
  ) %>%
  dplyr::select(StationID, x_fin, y_fin) %>%
  dplyr::group_by(StationID) %>%
  dplyr::mutate(ID = row_number()) %>%
  dplyr::arrange(StationID, ID) %>%
  dplyr::ungroup()

# Calcul des points médians entre paires successives
points_mid <- points_df %>%
  group_by(StationID) %>%
  mutate(x_next = lead(x_fin), y_next = lead(y_fin)) %>%
  filter(!is.na(x_next) & !is.na(y_next)) %>%
  mutate(x_mid = (x_fin + x_next) / 2, y_mid = (y_fin + y_next) / 2) %>%
  dplyr::select(StationID, ID, x_mid, y_mid)

# Extraction des profondeurs sur les points médians
coords_mid <- cbind(points_mid$x_mid, points_mid$y_mid)
points_mid$valeur_raster <- terra::extract(map[["profondeur.pred"]], coords_mid)[[1]]

# Calcul des profondeurs min et max pour chaque station
depth_by_station <- points_mid %>%
  filter(!is.na(valeur_raster)) %>%
  group_by(StationID) %>%
  summarise(
    prof_min = min(valeur_raster),
    prof_max = max(valeur_raster),
    diff_profondeur = prof_max - prof_min,
    .groups = "drop"
  )

depth_by_station <- depth_by_station %>%
  mutate(StationID = as.character(StationID))

data_all_cap <- data_all_cap_nonnull %>%
  mutate(StationID = as.character(StationID)) %>%
  left_join(depth_by_station, by = "StationID")



### Extraction des valeurs de pente ---------------------------------------------

# Copier points_mid dans points_mid_slope
points_mid_slope <- points_mid

# Charger le raster de pente
slope_map <- rast(here(path_rast_bathy, 'ras_slope_pecheker_interpolCRO.tif'))

# Extraction des valeurs de pente
coords_mid <- cbind(points_mid_slope$x_mid, points_mid_slope$y_mid)
points_mid_slope$valeur_raster <- terra::extract(slope_map[["slope"]], coords_mid)[[1]]

# Préparer le tableau des pentes par station et par point (avec ID)
slope_by_station <- points_mid_slope %>%
  dplyr::select(StationID, ID, valeur_raster) %>%
  dplyr::rename(slope = valeur_raster)







### Calul des flags -------------------------------------------------------------

# Nettoyage des valeurs infinies
data_all_cap <- data_all_cap %>%
  mutate(across(
    c(tps_trempe_hm, distance, diff_taille, diff_profondeur),
    ~ ifelse(is.infinite(.), NA, .)
  ))

# Conversion du temps de trempe "hh:mm" en nombre d'heures
data_all_cap <- data_all_cap %>%
  mutate(
    tps_trempe_hm_num = as.numeric(lubridate::hm(tps_trempe_hm)) / 3600
  )

# Calcul des quantiles pour les variables de contrôle qualité
q_vals <- data_all_cap %>%
  summarise(
    q95_tps = quantile(tps_trempe_hm_num, probs = 0.95, na.rm = TRUE),
    q05_tps = quantile(tps_trempe_hm_num, probs = 0.05, na.rm = TRUE),
    q95_dist = quantile(distance, probs = 0.95, na.rm = TRUE),
    q05_dist = quantile(distance, probs = 0.05, na.rm = TRUE),
    q95_taille = quantile(diff_taille, probs = 0.95, na.rm = TRUE),
    q05_taille = quantile(diff_taille, probs = 0.05, na.rm = TRUE),
    q95_prof = quantile(diff_profondeur, probs = 0.95, na.rm = TRUE),
    q05_prof = quantile(diff_profondeur, probs = 0.05, na.rm = TRUE)
  )

# Création du flag binaire (1 si la condition est dans l’intervalle, 0 sinon)
data_flag <- data_all_cap %>%
  mutate(
    flag =
      replace_na(as.integer(distance >= q_vals$q05_dist & distance <= q_vals$q95_dist), 0) +
      replace_na(as.integer(diff_taille >= q_vals$q05_taille & diff_taille <= q_vals$q95_taille), 0) +
      replace_na(as.integer(diff_profondeur <= q_vals$q95_prof), 0)
  )

# Création d’un tableau de détail des conditions de flag par variable
df_flag_detail <- data_all_cap %>%
  mutate(
    cond_tps    = replace_na(as.integer(tps_trempe_hm_num >= q_vals$q05_tps & tps_trempe_hm_num <= q_vals$q95_tps), 0),
    cond_dist   = replace_na(as.integer(distance         >= q_vals$q05_dist & distance         <= q_vals$q95_dist), 0),
    cond_taille = replace_na(as.integer(diff_taille      >= q_vals$q05_taille & diff_taille    <= q_vals$q95_taille), 0),
    cond_prof   = replace_na(as.integer(diff_profondeur  <= q_vals$q95_prof), 0),
    flag_total  = cond_tps + cond_dist + cond_taille + cond_prof
  ) %>%
  dplyr::select(StationID, cond_tps, cond_dist, cond_taille, cond_prof, flag_total)

# Choix du jeu de données pour le modèle ----------------------------------------

# Extraction des flags
data_modelo = data_flag %>%
  dplyr::select(StationID, flag)


points_mid <- points_mid %>%
  mutate(StationID = as.character(StationID))


# Ajout des coordonnées et profondeur bathy
data_modelo = left_join(data_modelo, points_mid, by = c("StationID"))

# Ajout des données de captures (CPUE)
data_CPUE <- data_all_cap %>%
  dplyr::select(
    StationID,
    ScanNumber,
    ScanWeightKg,
    HookHaul,
    TargetSpecies,
    HaulEndDate
  )

data_modelo = left_join(data_modelo, data_CPUE, by = c("StationID"))


data_modelo <- data_modelo %>%
  mutate(Year = year(HaulEndDate))

# Renommage et transformation finale
data_modelo = data_modelo %>%
  rename(x = x_mid,
         y = y_mid,
         prof_bathy = valeur_raster) %>%
  mutate()

# Ajout des pentes
slope_by_station <- slope_by_station %>%
  mutate(StationID = as.character(StationID))

data_modelo <- left_join(data_modelo, slope_by_station, by = c("StationID", "ID"))


# Ajout du temps de trempe numérique
data_tps = data_all_cap %>%
  dplyr::select(StationID, tps_trempe_hm_num)

data_modelo <- left_join(data_modelo, data_tps, by = c("StationID"))

data_modelo = left_join(data_modelo, df_flag_detail, by = c("StationID"))

# Sauvegarde de l’objet final

saveRDS(data_modelo, file = here("data", output_df, "data_modelo-crozet.rds"))