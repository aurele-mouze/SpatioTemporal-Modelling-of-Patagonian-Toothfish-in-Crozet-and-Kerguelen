
add_boundary_to_mesh <- function(EC_shp, simplify_keep) {
  
  EC_sf <- sf::st_as_sf(EC_shp)
  EC_sf_utm <- sf::st_transform(EC_sf, utm_crs_used)
  
  # Merge polygons into a single valid geometry
  EC_boundary_m <- EC_sf_utm %>%
    sf::st_make_valid() %>%
    dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop")
  
  EC_boundary_simple_m <- rmapshaper::ms_simplify(
    EC_boundary_m,
    keep = simplify_keep,
    keep_shapes = TRUE
  )
  
  # Convert meters -> kilometers because make_mesh uses X/Y units
  EC_boundary_km <- EC_boundary_simple_m
  sf::st_geometry(EC_boundary_km) <- sf::st_geometry(EC_boundary_km) / 1000
  
  EC_sp <- methods::as(methods::as(EC_boundary_km, "Spatial"), "SpatialPolygons")
  bnd_segm <- fmesher::fm_as_segm(EC_sp)
  
  return(bnd_segm)
}