CRO_map <- ggplot() +
  xlim(48, 53) +
  ylim(-48, -45) +
  xlab("Longitude") +
  ylab("Latitude") +
  
  geom_sf(
    data = Isobath_CRO_sf %>% dplyr::filter(level == "-500"),
    colour = "grey80",
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_sf(
    data = Isobath_CRO_sf %>% dplyr::filter(level == "-1000"),
    colour = "grey60",
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_sf(
    data = Isobath_CRO_sf %>% dplyr::filter(level == "-1500"),
    colour = "grey30",
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_sf(
    data = Isobath_CRO_sf %>% dplyr::filter(level == "-2000"),
    colour = "black",
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  
  geom_sf(
    data = Isobath_CRO_sf %>% dplyr::filter(level == "-2500"),
    colour = "grey80",
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_sf(
    data = Isobath_CRO_sf %>% dplyr::filter(level == "-3000"),
    colour = "grey60",
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_sf(
    data = Isobath_CRO_sf %>% dplyr::filter(level == "-3500"),
    colour = "grey30",
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  
  geom_sf(
    data = EEZ_sf,
    fill = "transparent",
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  
  geom_sf(
    data = FR_island_sf,
    fill = "darkgrey",
    inherit.aes = FALSE
  ) +
  
  theme_bw()

CRO_map





bathy <- rast("C:/Users/Aurele/Documents/GitHub/stage_AM/Data/gebco_2025_n-27.1582_s-75.8145_w25.1719_e116.1563.nc")
bathy
plot(bathy)
bathy_band <- ifel(bathy <= -500 & bathy >= -2500, 1, NA)
bathy_poly <- as.polygons(bathy_band, dissolve = TRUE)
bathy_poly <- st_as_sf(bathy_poly)

CRO_map +
  geom_sf(
    data = bathy_poly,
    fill = "lightblue",
    colour = NA,
    alpha = 0.5
  )
