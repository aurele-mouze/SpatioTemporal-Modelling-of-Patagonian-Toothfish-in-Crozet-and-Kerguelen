# ------------------------------------------------------------------------------#
# Function used to extract only converged models from the diagnostics object
# ------------------------------------------------------------------------------#

get_converged_models <- function(all_fitted_models) {
  
  converged_models <- list()
  
  # Loop over mesh configurations
  for (mesh_name in names(all_fitted_models)) {
    
    mesh_list <- all_fitted_models[[mesh_name]]
    if (is.null(mesh_list)) next
    
    # Loop over responses
    for (response in names(mesh_list)) {
      
      response_list <- mesh_list[[response]]
      
      if (is.null(response_list$models_valid)) next
      
      models_valid <- response_list$models_valid
      
      # Keep only non-empty lists
      if (length(models_valid) > 0) {
        converged_models[[mesh_name]][[response]] <- models_valid
      }
    }
  }
  
  converged_models
}

# converged_models <- get_converged_models(sanity_by_region)
