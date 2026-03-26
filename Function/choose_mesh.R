#Fonction test des différent mesh possible pour le modèle

choose_mesh <- function(data, 
                        formula_sdm   = get("formula_sdm", envir = .GlobalEnv),
                        distrib_family = get("distrib_family", envir = .GlobalEnv),
                        spatial_var   = get("spatial_var", envir = .GlobalEnv),
                        knots_range) {
  
  results <- data.frame(
    n_knots = integer(),
    all_ok  = logical(),
    AIC     = numeric()
  )
  
  for (i in knots_range) {
    cat("Essai avec", i, "nœuds...\n")
    
    # Encapsuler toute la tentative dans tryCatch
    attempt <- tryCatch({
      # Création du maillage
      mesh <- make_mesh(data, c("X", "Y"), n_knots = i)
      
      # Polygone du domaine
      bnb <- fm_nonconvex_hull(cbind(data$X, data$Y), convex = -0.05)
      coords <- bnb[[1]][[1]]
      coords_unique <- coords[!duplicated(coords), ]
      coords_closed <- rbind(coords_unique, coords_unique[1, ])
      
      bnb_sf <- st_sf(
        geometry = st_sfc(st_polygon(list(coords_closed))),
        crs = "+proj=utm +zone=42 +datum=WGS84 +units=km +no_defs"
      )
      
      # Coordonnées du mesh
      mesh_coords <- as.data.frame(mesh$loc_xy)
      colnames(mesh_coords) <- c("X", "Y")
      
      buffer_zone <- st_buffer(
        st_convex_hull(st_union(st_as_sf(mesh_coords, coords = c("X", "Y"), crs = st_crs(bnb_sf)))),
        dist = 5
      )
      
      barrier_geom <- st_difference(buffer_zone, st_union(bnb_sf))
      barrier_sf <- st_sf(geometry = st_sfc(barrier_geom), crs = st_crs(bnb_sf))
      
      bspde <- sdmTMBextra::add_barrier_mesh(
        mesh,
        barrier_sf,
        range_fraction = 0.1,
        plot = FALSE
      )
      
      # Modèle
      m_gamma <- sdmTMB(
        formula = formula_sdm,
        data    = data,
        mesh    = bspde,
        family  = distrib_family,
        spatial = spatial_var
      )
      
      test <- sanity(m_gamma)
      
      data.frame(
        n_knots = i,
        all_ok  = test$all_ok,
        AIC     = AIC(m_gamma)
      )
    },
    error = function(e) {
      cat("⚠️ Erreur avec", i, "nœuds :", conditionMessage(e), "\n")
      # Retourner une ligne "erreur"
      data.frame(
        n_knots = i,
        all_ok  = FALSE,
        AIC     = NA
      )
    })
    
    results <- rbind(results, attempt)
    
    if (attempt$all_ok) {
      cat("Modèle OK avec", i, "nœuds ✅\n")
    } else {
      cat("Modèle NON valide ou erreur avec", i, "nœuds ❌\n")
    }
  }
  
  # Filtrer les modèles valides
  valid_models <- results[results$all_ok == TRUE, ]
  
  return(as.data.frame(valid_models))
}
