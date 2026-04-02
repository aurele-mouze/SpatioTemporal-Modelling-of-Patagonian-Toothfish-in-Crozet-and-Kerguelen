###Script make grid -----------------------------------------------------------

map = terra::rast(here(path_rast_bathy, 'ras_pecheker_interpol_totCRO.tif'))


map_grid_m <- terra::project(map, "EPSG:32743")
print(map_grid_m)
template <- rast(
  ext = ext(map_grid_m),
  resolution = res,
  crs = crs(map_grid_m)
)

origin(template) <- origin(map_grid_m)
print(template)
cellSize(template)

map_grid_1km <- resample(map_grid_m, template, method = resampling_method)
print(map_grid_1km)



plot(map_grid_1km)

map_df <- terra::as.data.frame(map_grid_m, xy = TRUE)
map_df$slope[is.na(map_df$slope)] <- 0

head(map_df)

saveRDS(map_df, file = here("data", output_df, "grid_extrapolCRO.rds"))

