#Figure de l'évolution de l'indice d'abondance relatif

index_mean <- index_depred_ker %>%
  mutate(
    mean_index = mean(est),
    est_rel = est / mean_index,
    lwr_rel = lwr / mean_index,
    upr_rel = upr / mean_index
  )

p_index <- ggplot(index_mean, aes(x = TAAFseason, y = est_rel)) +
  geom_hline(yintercept = 1, 
             linetype = "dashed", 
             colour = "grey40") +
  geom_ribbon(aes(ymin = lwr_rel, ymax = upr_rel),
              fill = "firebrick3",
              alpha = 0.25) +
  geom_line(colour = "firebrick3",
            linewidth = 1.1) +
  
  scale_x_continuous(breaks = index_mean$TAAFseason) +
  labs(
    x = "Année",
    y = "Indice d'abondance relatif"
  ) +
  theme_classic()



ggsave(
  here("data", output_plot, "indice_abondance_standardise_Kerguelen.png"),
  p_index,
  width = 10,
  height = 5,
  units = "in",
  dpi = 400
)

#Grille spatiale de l'indice d'abondance

grid_mean_ker <- pred_depred_ker %>%
  filter(TAAFseason >= 2008, TAAFseason <= 2025) %>%
  group_by(cell_id, X, Y) %>%
  summarise(
    est_mean = mean(est, na.rm = TRUE),
    est_sd = sd(est, na.rm = TRUE),
    n_years = n_distinct(TAAFseason),
    .groups = "drop"
  ) %>%
  mutate(
    indice_relatif = est_mean / mean(est_mean, na.rm = TRUE)
  )


# Passage en spatial
grid_mean_ker_sf <- grid_mean_ker %>%
  mutate(
    X_utm = X * 1000,
    Y_utm = Y * 1000
  ) %>%
  st_as_sf(
    coords = c("X_utm", "Y_utm"),
    crs = 32742
  )


# Passage en longitude / latitude
grid_mean_ker_ll <- st_transform(
  grid_mean_ker_sf,
  4326
)

grid_mean_ker_ll <- cbind(
  st_coordinates(grid_mean_ker_ll),
  st_drop_geometry(grid_mean_ker_ll)
)

names(grid_mean_ker_ll)[1:2] <- c("lon", "lat")


# Carte
cell_width <- 0.0265
cell_height <- 0.0265


q05 <- quantile(
  grid_mean_ker_ll$indice_relatif,
  0.01,
  na.rm = TRUE
)

q95 <- quantile(
  grid_mean_ker_ll$indice_relatif,
  0.99,
  na.rm = TRUE
)

Spatial_Ker <- KER_map +
  geom_tile(
    data = grid_mean_ker_ll,
    aes(
      x = lon,
      y = lat,
      fill = indice_relatif
    ),
    width = cell_width,
    height = cell_height,
    alpha = 0.7
  ) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 1,
    limits = c(q05, q95),
    oob = scales::squish,
    name = "Indice relatif\n(moyenne = 1)"
  ) +
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    style = north_arrow_fancy_orienteering
  ) +
  annotation_scale(
    location = "bl",
    width_hint = 0.25,
    text_cex = 0.8
  ) +
  theme(
    legend.position = c(0.90, 0.15),
    legend.background = element_rect(
      fill = "white",
      color = "black"
    )
  )


ggsave(
  here("data", output_plot, "CPUE_pred_Ker.png"),
  Spatial_Ker,
  width = 10,
  height = 10,
  units = "in",
  dpi = 400
)




# Hotspots spatiaux

# 1. Hotspots initiaux


grid_hotspot_ker <- pred_depred_ker %>%
  filter(TAAFseason >= 2008, TAAFseason <= 2025) %>%
  group_by(cell_id, X, Y) %>%
  summarise(
    est_mean = mean(est, na.rm = TRUE),
    n_years = n_distinct(TAAFseason),
    .groups = "drop"
  )

Q90_ker <- quantile(
  grid_hotspot_ker$est_mean,
  probs = 0.90,
  na.rm = TRUE
)

grid_hotspot_ker <- grid_hotspot_ker %>%
  mutate(
    hotspot = est_mean >= Q90_ker
  )


# 2. Filtre pour garder uniquement les cellules hotspot


hotspots <- grid_hotspot_ker %>%
  filter(hotspot)



# 3. Recherche des cellules voisines


# Distance entre deux centres de cellules
dx <- median(diff(sort(unique(hotspots$X))))
dy <- median(diff(sort(unique(hotspots$Y))))

# Tolérance pour éviter les problèmes numériques
tol <- 1e-6

# On cherche les voisins dans les 8 directions
neighbors <- expand.grid(
  dx = c(-dx, 0, dx),
  dy = c(-dy, 0, dy)
) %>%
  filter(!(dx == 0 & dy == 0))


# Création des liens entre cellules voisines
edges <- lapply(seq_len(nrow(neighbors)), function(i) {
  
  hotspots %>%
    select(cell_id, X, Y) %>%
    inner_join(
      hotspots %>%
        select(cell_id, X, Y) %>%
        rename(
          cell_id_2 = cell_id,
          X2 = X,
          Y2 = Y
        ),
      by = character()
    ) %>%
    filter(
      abs((X2 - X) - neighbors$dx[i]) < tol,
      abs((Y2 - Y) - neighbors$dy[i]) < tol
    ) %>%
    transmute(
      from = cell_id,
      to = cell_id_2
    )
}) %>%
  bind_rows() %>%
  distinct()



# 4. Identification des groupes de cellules connectées


g <- graph_from_data_frame(
  edges,
  directed = FALSE,
  vertices = hotspots$cell_id
)

components_g <- components(g)


hotspots <- hotspots %>%
  mutate(
    hotspot_id = components_g$membership[
      match(cell_id, names(components_g$membership))
    ]
  )


# 5. Taille de chaque hotspot


hotspot_sizes <- hotspots %>%
  count(hotspot_id, name = "n_cells")


# 6. Conservation des hotspots > 30 cellules 


hotspots_final <- hotspots %>%
  left_join(
    hotspot_sizes,
    by = "hotspot_id"
  ) %>%
  mutate(
    hotspot_final = n_cells > 30 #peut modifier le nombre de cellules
  )


# Vérification
table(hotspots_final$hotspot_final)

hotspots_final %>%
  filter(hotspot_final) %>%
  count(hotspot_id, n_cells)

grid_hotspot_ker_ll <- hotspots_final %>%
  filter(hotspot_final) %>%
  mutate(
    X_utm = X * 1000,
    Y_utm = Y * 1000
  ) %>%
  st_as_sf(
    coords = c("X_utm", "Y_utm"),
    crs = 32742,
    remove = FALSE
  ) %>%
  st_transform(4326) %>%
  mutate(
    lon = st_coordinates(.)[, 1],
    lat = st_coordinates(.)[, 2]
  ) %>%
  st_drop_geometry()

#Création de la figure
Hotspot_Ker <- KER_map +
  geom_tile(
    data = grid_hotspot_ker_ll,
    aes(
      x = lon,
      y = lat
    ),
    width = 0.0265,
    height = 0.0265,
    fill = "red",
    alpha = 0.8
  ) +
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    style = north_arrow_fancy_orienteering
  ) +
  annotation_scale(
    location = "bl",
    width_hint = 0.25,
    text_cex = 0.8
  ) +
  theme(
    legend.position = "none"
  )


ggsave(
  here("data", output_plot, "Hotspot_Ker.png"),
  Hotspot_Ker,
  width = 10,
  height = 10,
  units = "in",
  dpi = 400
)



#Grille spatio-temporelle de l'indice d'abondance

#Peut modifier les périodes
grid_period_ker <- pred_depred_ker %>%
  filter(TAAFseason >= 2008, TAAFseason <= 2025) %>%
  mutate(
    periode = case_when(
      TAAFseason >= 2008 & TAAFseason <= 2013 ~ "2008-2013",
      TAAFseason >= 2014 & TAAFseason <= 2019 ~ "2014-2019",
      TAAFseason >= 2020 & TAAFseason <= 2025 ~ "2020-2025"
    )
  ) %>%
  group_by(periode, cell_id, X, Y) %>%
  summarise(
    est_mean = mean(est, na.rm = TRUE),
    est_sd = sd(est, na.rm = TRUE),
    n_years = n_distinct(TAAFseason),
    .groups = "drop"
  ) %>%
  group_by(periode) %>%
  mutate(
    indice_relatif = est_mean / mean(est_mean, na.rm = TRUE)
  ) %>%
  ungroup()

# Passage en spatial
grid_period_ker_sf <- grid_period_ker %>%
  mutate(
    X_utm = X * 1000,
    Y_utm = Y * 1000
  ) %>%
  st_as_sf(
    coords = c("X_utm", "Y_utm"),
    crs = 32742
  )


# Passage en longitude / latitude
grid_period_ker_ll <- st_transform(
  grid_period_ker_sf,
  4326
)

grid_period_ker_ll <- cbind(
  st_coordinates(grid_period_ker_ll),
  st_drop_geometry(grid_period_ker_ll)
)

names(grid_period_ker_ll)[1:2] <- c("lon", "lat")


# Carte
cell_width <- 0.0264
cell_height <- 0.0264


q01 <- quantile(
  grid_period_ker_ll$indice_relatif,
  0.01,
  na.rm = TRUE
)

q99 <- quantile(
  grid_period_ker_ll$indice_relatif,
  0.99,
  na.rm = TRUE
)

Spatial_Ker_periods <- KER_map +
  geom_tile(
    data = grid_period_ker_ll,
    aes(
      x = lon,
      y = lat,
      fill = indice_relatif
    ),
    width = cell_width,
    height = cell_height,
    alpha = 0.7
  ) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 1,
    limits = c(q01, q99),
    oob = scales::squish,
    name = "Indice relatif\n(moyenne = 1)"
  ) +
  facet_wrap(
    ~periode,
    nrow = 1
  ) +
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    style = north_arrow_fancy_orienteering
  ) +
  annotation_scale(
    location = "bl",
    width_hint = 0.25,
    text_cex = 0.8
  ) +
  theme(
    legend.position = c(0.93, 0.15),
    legend.background = element_rect(
      fill = "white",
      color = "black"
    ),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(0.3, "cm"),
    strip.text = element_text(
      size = 11,
      face = "bold"
    )
  ) +
  guides(
    fill = guide_colorbar(
      barwidth = 0.35,
      barheight = 2.2
    )
  )


ggsave(
  here("data", output_plot, "CPUE_spatiotemp_Kerguelen.png"),
  Spatial_Ker_periods,
  width = 12,
  height = 4.8,
  units = "in",
  dpi = 500
)


# Hotspots par période


# 1. Moyenne des prédictions par période et cellule


grid_hotspot_period_ker <- pred_depred_ker %>%
  filter(TAAFseason >= 2008, TAAFseason <= 2025) %>%
  mutate(
    periode = case_when(
      TAAFseason >= 2008 & TAAFseason <= 2013 ~ "2008-2013",
      TAAFseason >= 2014 & TAAFseason <= 2019 ~ "2014-2019",
      TAAFseason >= 2020 & TAAFseason <= 2025 ~ "2020-2025"
    )
  ) %>%
  group_by(periode, cell_id, X, Y) %>%
  summarise(
    est_mean = mean(est, na.rm = TRUE),
    n_years = n_distinct(TAAFseason),
    .groups = "drop"
  ) %>%
  mutate(
    cpue_mean = 1000 * est_mean / 4620
  )


# 2. Q90 calculé séparément pour chaque période


Q90_period_ker <- grid_hotspot_period_ker %>%
  group_by(periode) %>%
  summarise(
    Q90 = quantile(
      est_mean,
      probs = 0.90,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


grid_hotspot_period_ker <- grid_hotspot_period_ker %>%
  left_join(
    Q90_period_ker,
    by = "periode"
  ) %>%
  mutate(
    hotspot = est_mean >= Q90
  )


# 3. Fonction pour identifier les groupes connectés


find_groups <- function(nodes, edges) {
  
  nodes <- as.character(nodes)
  edges$cell_id <- as.character(edges$cell_id)
  edges$cell_id_2 <- as.character(edges$cell_id_2)
  
  groups <- setNames(
    integer(length(nodes)),
    nodes
  )
  
  group_id <- 0
  
  for (node in nodes) {
    
    # Si déjà attribué à un groupe
    if (groups[node] != 0) {
      next
    }
    
    group_id <- group_id + 1
    
    current <- node
    groups[node] <- group_id
    
    while (length(current) > 0) {
      
      neighbours <- unique(c(
        edges$cell_id_2[
          edges$cell_id %in% current
        ],
        edges$cell_id[
          edges$cell_id_2 %in% current
        ]
      ))
      
      # On ne garde que les cellules pas encore attribuées
      neighbours <- neighbours[
        neighbours %in% nodes &
          groups[neighbours] == 0
      ]
      
      if (length(neighbours) == 0) {
        current <- character(0)
      } else {
        groups[neighbours] <- group_id
        current <- neighbours
      }
    }
  }
  
  groups
}

# 4. Identification des groupes par période


hotspots_period <- grid_hotspot_period_ker %>%
  filter(hotspot)


hotspots_period_list <- lapply(
  split(hotspots_period, hotspots_period$periode),
  function(dat) {
    
    # Distance entre les cellules
    dx <- median(diff(sort(unique(dat$X))))
    dy <- median(diff(sort(unique(dat$Y))))
    
    # Index de grille
    dat <- dat %>%
      mutate(
        ix = round(X / dx),
        iy = round(Y / dy)
      )
    
    # Recherche des voisins
    edges <- dat %>%
      select(cell_id, ix, iy) %>%
      cross_join(
        dat %>%
          select(
            cell_id_2 = cell_id,
            ix2 = ix,
            iy2 = iy
          )
      ) %>%
      filter(
        abs(ix2 - ix) <= 1,
        abs(iy2 - iy) <= 1,
        !(ix2 == ix & iy2 == iy)
      ) %>%
      select(
        cell_id,
        cell_id_2
      )
    
    # Groupes connectés
    group_ids <- find_groups(
      nodes = dat$cell_id,
      edges = edges
    )
    
    dat %>%
      mutate(
        hotspot_id = group_ids[cell_id]
      )
  }
) %>%
  bind_rows()


# 5. Taille des hotspots


hotspot_sizes_period <- hotspots_period_list %>%
  count(
    periode,
    hotspot_id,
    name = "n_cells"
  )


# 6. Hotspots finaux > 10 cellules

hotspots_period_final <- hotspots_period_list %>%
  left_join(
    hotspot_sizes_period,
    by = c("periode", "hotspot_id")
  ) %>%
  mutate(
    hotspot_final = n_cells > 30
  )

grid_hotspot_period_ker_ll <- hotspots_period_final %>%
  filter(hotspot_final) %>%
  mutate(
    X_utm = X * 1000,
    Y_utm = Y * 1000
  ) %>%
  st_as_sf(
    coords = c("X_utm", "Y_utm"),
    crs = 32742,
    remove = FALSE
  ) %>%
  st_transform(4326) %>%
  mutate(
    lon = st_coordinates(.)[, 1],
    lat = st_coordinates(.)[, 2]
  ) %>%
  st_drop_geometry()

#Création de la figure
Hotspot_period_Ker <- KER_map +
  geom_tile(
    data = grid_hotspot_period_ker_ll,
    aes(
      x = lon,
      y = lat
    ),
    width = 0.0265,
    height = 0.0265,
    fill = "red",
    alpha = 0.5
  ) +
  facet_wrap(~periode) +
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    style = north_arrow_fancy_orienteering
  ) +
  annotation_scale(
    location = "bl",
    width_hint = 0.25,
    text_cex = 0.8
  ) +
  theme(
    legend.position = "none"
  )


ggsave(
  here("data", output_plot, "Hotspot_spatiotemp_Kerguelen.png"),
  Hotspot_period_Ker,
  width = 12,
  height = 4.8,
  units = "in",
  dpi = 500
)
