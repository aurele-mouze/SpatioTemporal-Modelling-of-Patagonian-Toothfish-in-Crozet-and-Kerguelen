
# ============================================================================== #
####  CREATE THE SPATIAL MESH FOR THE SPDE APPROACH ####
# ============================================================================== #


# Function that builds a mesh from CGFS data
create_meshes_configurations <- function(data_m,
                                         cutoff_values,
                                         max_edge_in,
                                         max_edge_out,
                                         utm_crs_used,
                                         land_domain) {
  
  configs <- vector("list", length(cutoff_values))
  names(configs) <- paste0("cutoff_", cutoff_values)
  
  for (i in seq_along(cutoff_values)) {
    
    cutoff <- cutoff_values[i]
    max_edge <- c(max_edge_in[i], max_edge_out[i])
    
    message(
      "Testing cutoff = ", cutoff,
      " ; max.edge = (", max_edge[1], ", ", max_edge[2], ")"
    )
    
    mesh <- make_mesh(
      data_m,
      c("X", "Y"),
      fmesher_func = fmesher::fm_mesh_2d_inla,
      max.edge = max_edge,
      cutoff = cutoff
    )
    
    # ------------------------------------------------------------------------------#
    #### 4. Plot the mesh with land areas ####
    # ------------------------------------------------------------------------------#
    
    mesh_m <- mesh$mesh
    mesh_m$loc <- mesh_m$loc * 1000
    
    mesh_m_bbox <- st_bbox(
      c(
        xmin = min(mesh_m$loc[, 1]),
        xmax = max(mesh_m$loc[, 1]),
        ymin = min(mesh_m$loc[, 2]),
        ymax = max(mesh_m$loc[, 2])
      ),
      crs = utm_crs_used
    )
    
    land_mesh <- land_domain %>%
      st_transform(utm_crs_used) %>%
      st_crop(mesh_m_bbox)
    
    plot_mesh <- ggplot() +
      geom_sf(data = land_mesh, fill = "grey80", color = NA) +
      gg(mesh_m, edge.color = "grey50", edge.linewidth = 0.2) +
      geom_spatial_point(
        data = data_m,
        aes(lon, lat),
        color = "black",
        size = 1.5,
        alpha = 0.5
      ) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
      coord_sf(crs = st_crs(land_mesh)) +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.ticks = element_line(colour = "black")
      ) +
      labs(x = "", y = "", title = paste0("Mesh with cutoff = ", round(cutoff, 1), " km"))
    
    configs[[i]] <- list(
      cutoff = cutoff,
      max_edge = max_edge,
      mesh = mesh,
      plot_mesh = plot_mesh
    )
  }
  
  configs
}

# ------------------------------------------------------------------------------#
#### 5. Build meshes by region depending on user-defined flags ####
# ------------------------------------------------------------------------------#
# 
# mesh_by_region <- list()
# 
# If East_English_Channel is TRUE, build the east mesh
# if (isTRUE(East_English_Channel)) {
#   mesh_by_region$east <- create_mesh_by_region(region_name = "east")
# }
# 
# # If West_English_Channel is TRUE, build the west mesh
# if (isTRUE(West_English_Channel)) {
#   mesh_by_region$west <- create_mesh_by_region(region_name = "west")
# }
# 





