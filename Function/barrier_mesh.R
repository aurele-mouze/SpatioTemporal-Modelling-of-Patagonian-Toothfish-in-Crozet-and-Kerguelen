
mesh_barriere = function(data, nb_noeud){
mesh <- sdmTMB::make_mesh(data, c("X", "Y"), n_knots = nb_noeud)  # Création du maillage
# ---- Barrière ----
bnb = fm_nonconvex_hull(cbind(data$X, data$Y), convex = -0.05)  # Enveloppe concave
coords <- bnb[[1]][[1]]
coords_unique <- coords[!duplicated(coords), ]
coords_closed <- rbind(coords_unique, coords_unique[1, ])  # Fermeture du polygone

bnb_sf <- st_sf(
  geometry = st_sfc(st_polygon(list(coords_closed))),
  crs = "+proj=utm +zone=42 +datum=WGS84 +units=km +no_defs"
)

mesh_coords <- as.data.frame(mesh$loc_xy)

# Création d'une zone tampon pour la barrière
buffer_zone <- st_buffer(
  st_convex_hull(st_union(st_as_sf(mesh_coords, coords = c("X", "Y"), crs = st_crs(bnb_sf)))),
  dist = 5
)

barrier_geom <- st_difference(buffer_zone, st_union(bnb_sf))
barrier_sf <- st_sf(geometry = st_sfc(barrier_geom), crs = st_crs(bnb_sf))

# ---- Ajouter la barrière au maillage ----
bspde <- sdmTMBextra::add_barrier_mesh(
  mesh,
  barrier_sf,
  range_fraction = 0.1,
  plot = TRUE
)
return(bspde)
}

