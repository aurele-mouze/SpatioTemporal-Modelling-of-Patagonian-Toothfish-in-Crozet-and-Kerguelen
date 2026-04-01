###Script make grid -----------------------------------------------------------

map = terra::rast(here(path_rast_bathy, 'ras_pecheker_interpol_tot.tif'))


map_grid_m <- terra::project(map, "EPSG:32742")
print(map_grid_m)
template <- rast(ext(map_grid_m), resolution = res, crs = crs(map_grid_m))
print(template)
cellSize(template)

map_grid_1km <- resample(map_grid_m, template, method = resampling_method)
print(map_grid_1km)



plot(map_grid_1km)

map_df = terra::as.data.frame(map_grid_1km, xy = TRUE)
map_df <- map_df %>% filter(!is.na(profondeur.pred))
map_df$slope[is.na(map_df$slope)] <- 0

head(map_df)

saveRDS(map_df, file = here("data", output_df, "grid_extrapol.rds"))

