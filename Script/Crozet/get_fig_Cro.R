index_mean <- index_depred_cro %>%
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
  filename = "indice_abondance_standardise_Crozet.png",
  plot = p_index,
  width = 8,
  height = 5,
  dpi = 300
)

grid_mean_cro <- pred_depred_cro$data %>%
  group_by(cell_id, X, Y) %>%
  summarise(
    est_mean = mean(est, na.rm = TRUE),
    est_sd = sd(est, na.rm = TRUE),
    n_years = n(),
    .groups = "drop"
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

cell_width <-0.0265
cell_height <- 0.0265
CRO_map +
  geom_tile(
    data = grid_mean_cro_ll,
    aes(
      x = lon,
      y = lat,
      fill = est_mean
    ),
    width = cell_width,
    height = cell_height,
    alpha = 0.9
  ) +
  scale_fill_viridis_c(
    option = "viridis",
    name = "CPUE prédite\n(kg/hameçon)"
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
  )+0
  theme(
    legend.position = c(0.90, 0.15),
    legend.background = element_rect(
      fill = "white",
      color = "black"
    )
  )





grid_period_cro <- pred_depred_cro$data %>%
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
    n_years = n(),
    .groups = "drop"
  )
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

CRO_map +
  geom_tile(
    data = grid_period_cro_ll,
    aes(
      x = lon,
      y = lat,
      fill = est_mean
    ),
    width = cell_width,
    height = cell_height,
    alpha = 0.9
  ) +
  scale_fill_viridis_c(
    option = "viridis",
    name = "CPUE prédite\n(kg/hameçon)"
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
  )+
  theme(
    legend.position = c(0.96, 0.17),
    legend.background = element_rect(
      fill = "white",
      color = "black"
    )
  )
