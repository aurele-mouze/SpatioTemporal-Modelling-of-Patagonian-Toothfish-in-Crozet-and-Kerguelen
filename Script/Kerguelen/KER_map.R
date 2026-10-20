
KER_map <- ggplot() + xlim(63.6,75.9) + ylim(-52.4,-45.3) + xlab("Longitude") + ylab("Latitude") +
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
  geom_sf(
    data = Isobath_KER_sf  %>% dplyr::filter(level %in% c("-1500")),
    colour = "grey30",
    inherit.aes = F,
    size = 0.5
  ) +
  geom_sf(
    data = Isobath_KER_sf  %>% dplyr::filter(level %in% c("-2000")),
    colour = "black",
    inherit.aes = F,
    size = 0.5
  ) +
  geom_sf(
    data = Isobath_KER_sf  %>% dplyr::filter(level %in% c("-2500")),
    colour = "grey80",
    inherit.aes = F,
    size = 0.5
  ) +
  geom_sf(
    data = Isobath_KER_sf  %>% dplyr::filter(level %in% c("-3000")),
    colour = "grey60",
    inherit.aes = F,
    size = 0.5
  ) +
  geom_sf(
    data = Isobath_KER_sf  %>% dplyr::filter(level %in% c("-3500")),
    colour = "grey30",
    inherit.aes = F,
    size = 0.5
  ) +
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
  theme_bw()



KER_SHELF_map <- ggplot() + xlim(63.6, 73) + ylim(-52, -46) + xlab("Longitude") + ylab("Latitude") +
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
  theme_bw()




