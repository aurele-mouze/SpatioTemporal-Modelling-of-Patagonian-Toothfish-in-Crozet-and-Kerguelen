index_mean <- index_depred_cro %>%
  filter(TAAFseason >= 2008, TAAFseason <= 2025) %>%
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
              fill = "navyblue",
              alpha = 0.25) +
  geom_line(colour = "navyblue",
            linewidth = 1.1) +
 
  scale_x_continuous(breaks = index_mean$TAAFseason) +
  labs(
    x = "Année",
    y = "Indice d'abondance relatif"
  ) +
  theme_classic()


p_index

ggsave(
  here("data", output_plot, "indice_abondance_standardise_Crozet.png"),
  p_index,
  width = 12,
  height = 4.8,
  units = "in",
  dpi = 300
)

grid_mean_cro <- pred_depred_cro %>%
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


grid_mean_cro_sf <- grid_mean_cro %>%
  mutate(
    X_utm = X * 1000,
    Y_utm = Y * 1000
  ) %>%
  st_as_sf(
    coords = c("X_utm", "Y_utm"),
    crs = 32743
  )

grid_mean_cro_ll <- st_transform(grid_mean_cro_sf, 4326)

grid_mean_cro_ll <- cbind(
  st_coordinates(grid_mean_cro_ll),
  st_drop_geometry(grid_mean_cro_ll)
)

names(grid_mean_cro_ll)[1:2] <- c("lon", "lat")

# Carte

cell_width <-0.0262
cell_height <- 0.0262
Spatial_Cro <- CRO_map +
  geom_tile(
    data = grid_mean_cro_ll,
    aes(
      x = lon,
      y = lat,
      fill = indice_relatif
    ),
    width = cell_width,
    height = cell_height,
    alpha = 0.8
  ) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 1,
    name = "Indice relatif\n(moyenne = 1)",
    guide = guide_colorbar(
      barwidth = 0.8,
      barheight = 4
    )
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
    legend.position = c(0.88, 0.20),
    legend.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.5
    ),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11)
  )

ggsave(
  here("data", output_plot, "CPUE_pred_Cro.png"),
  Spatial_Cro,
  width = 10,
  height = 10,
  units = "in",
  dpi = 400
)


grid_hotspot_cro <- pred_depred_cro %>%
  filter(TAAFseason >= 2008, TAAFseason <= 2025) %>%
  group_by(cell_id, X, Y) %>%
  summarise(
    est_mean = mean(est, na.rm = TRUE),
    n_years = n_distinct(TAAFseason),
    .groups = "drop"
  )

Q90 <- quantile(
  grid_hotspot_cro$est_mean,
  probs = 0.90,
  na.rm = TRUE
)

grid_hotspot_cro <- grid_hotspot_cro %>%
  mutate(
    hotspot = est_mean >= Q90
  )
names(grid_hotspot_cro)
table(grid_hotspot_cro$hotspot)

grid_hotspot_cro_ll <- grid_hotspot_cro %>%
  mutate(
    X_utm = X * 1000,
    Y_utm = Y * 1000
  ) %>%
  st_as_sf(
    coords = c("X_utm", "Y_utm"),
    crs = 32743,
    remove = FALSE
  ) %>%
  st_transform(4326) %>%
  mutate(
    lon = st_coordinates(.)[, 1],
    lat = st_coordinates(.)[, 2]
  ) %>%
  st_drop_geometry()

Hotspot_Cro <- CRO_map +
  geom_tile(
    data = filter(grid_hotspot_cro_ll, hotspot),
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
  here("data", output_plot, "Hotspot_Cro.png"),
  Hotspot_Cro,
  width = 10,
  height = 10,
  units = "in",
  dpi = 400
)



library(dplyr)
library(igraph)

# Uniquement les cellules Q90
hotspots <- grid_hotspot_cro %>%
  filter(hotspot) %>%
  select(cell_id, X, Y)

dx <- median(diff(sort(unique(grid_hotspot_cro$X))))
dy <- median(diff(sort(unique(grid_hotspot_cro$Y))))

dx
dy


hotspots_neighbors <- hotspots %>%
  mutate(
    X_left  = X - dx,
    X_right = X + dx,
    Y_down  = Y - dy,
    Y_up    = Y + dy
  )

n_dist <- 3

neighbor_offsets <- expand.grid(
  ix = -n_dist:n_dist,
  iy = -n_dist:n_dist
) %>%
  filter(!(ix == 0 & iy == 0)) %>%
  mutate(
    dx_offset = ix * dx,
    dy_offset = iy * dy
  )
pairs <- hotspots %>%
  tidyr::crossing(
    neighbor_offsets
  ) %>%
  mutate(
    X_neighbor = X + dx_offset,
    Y_neighbor = Y + dy_offset
  ) %>%
  select(
    cell_id_1 = cell_id,
    X_neighbor,
    Y_neighbor
  ) %>%
  left_join(
    hotspots %>%
      select(
        cell_id_2 = cell_id,
        X,
        Y
      ),
    by = c(
      "X_neighbor" = "X",
      "Y_neighbor" = "Y"
    )
  ) %>%
  filter(!is.na(cell_id_2))
g <- igraph::graph_from_data_frame(
  pairs %>%
    select(cell_id_1, cell_id_2),
  directed = FALSE,
  vertices = hotspots %>% select(cell_id)
)
components <- igraph::components(g)

hotspots <- hotspots %>%
  mutate(
    hotspot_id = components$membership[
      match(cell_id, names(components$membership))
    ]
  )
hotspot_sizes <- hotspots %>%
  count(
    hotspot_id,
    name = "n_cells"
  ) %>%
  arrange(desc(n_cells))

hotspots_final <- hotspots %>%
  left_join(
    hotspot_sizes,
    by = "hotspot_id"
  ) %>%
  filter(n_cells >= 10)
hotspots_final %>%
  distinct(hotspot_id, n_cells) %>%
  arrange(desc(n_cells))
grid_hotspot_cro <- grid_hotspot_cro %>%
  left_join(
    hotspots_final %>%
      select(
        cell_id,
        hotspot_id,
        n_cells
      ),
    by = "cell_id"
  ) %>%
  mutate(
    hotspot_final = !is.na(hotspot_id)
  )
grid_hotspot_cro_ll <- grid_hotspot_cro %>%
  mutate(
    X_utm = X * 1000,
    Y_utm = Y * 1000
  ) %>%
  st_as_sf(
    coords = c("X_utm", "Y_utm"),
    crs = 32743,
    remove = FALSE
  ) %>%
  st_transform(4326) %>%
  mutate(
    lon = st_coordinates(.)[, 1],
    lat = st_coordinates(.)[, 2]
  ) %>%
  st_drop_geometry()
CRO_map +
  geom_tile(
    data = filter(grid_hotspot_cro_ll, hotspot_final),
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




grid_period_cro <- pred_depred_cro %>%
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

grid_period_cro_sf <- grid_period_cro %>%
  mutate(
    X_utm = X * 1000,
    Y_utm = Y * 1000
  ) %>%
  st_as_sf(
    coords = c("X_utm", "Y_utm"),
    crs = 32743
  )

grid_period_cro_ll <- st_transform(grid_period_cro_sf, 4326)

grid_period_cro_ll <- cbind(
  st_coordinates(grid_period_cro_ll),
  st_drop_geometry(grid_period_cro_ll)
)



names(grid_period_cro_ll)[1:2] <- c("lon", "lat")

cell_width <- 0.0265
cell_height <- 0.0265

spatiotemp_Cro <- CRO_map +
  geom_tile(
    data = grid_period_cro_ll,
    aes(
      x = lon,
      y = lat,
      fill = indice_relatif
    ),
    width = cell_width,
    height = cell_height,
    alpha = 0.8
  ) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 1,
    name = "Indice relatif\n(moyenne = 1)",
    guide = guide_colorbar(
      barwidth = 0.4,
      barheight = 3
    )
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
    legend.position = c(0.96, 0.17),
    legend.background = element_rect(
      fill = "white",
      color = "black"
    ),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  )

ggsave(
  here("data", output_plot, "CPUE_spatiotemp_Crozet.png"),
  spatiotemp_Cro,
  width = 12,
  height = 4.8,
  units = "in",
  dpi = 500
)


grid_hotspot_period_cro <- pred_depred_cro %>%
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
  group_by(periode) %>%
  mutate(
    Q90 = quantile(
      est_mean,
      probs = 0.90,
      na.rm = TRUE
    ),
    hotspot = est_mean >= Q90
  ) %>%
  ungroup()

grid_hotspot_period_cro_ll <- grid_hotspot_period_cro %>%
  mutate(
    X_utm = X * 1000,
    Y_utm = Y * 1000
  ) %>%
  st_as_sf(
    coords = c("X_utm", "Y_utm"),
    crs = 32743,
    remove = FALSE
  ) %>%
  st_transform(4326) %>%
  mutate(
    lon = st_coordinates(.)[, 1],
    lat = st_coordinates(.)[, 2]
  ) %>%
  st_drop_geometry()

Hotspot_period_Cro <- CRO_map +
  geom_tile(
    data = filter(grid_hotspot_period_cro_ll, hotspot),
    aes(
      x = lon,
      y = lat
    ),
    width = 0.0265,
    height = 0.0265,
    fill = "red",
    alpha = 0.8
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
  here("data", output_plot, "Hotspot_spatiotemp_Crozet.png"),
  Hotspot_period_Cro,
  width = 12,
  height = 4.8,
  units = "in",
  dpi = 500
)
