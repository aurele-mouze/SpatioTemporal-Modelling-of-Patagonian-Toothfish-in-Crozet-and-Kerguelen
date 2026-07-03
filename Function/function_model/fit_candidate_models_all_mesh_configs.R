
    # ============================================================================== #
#### FIT CANDIDATE MODELS ####
# ============================================================================== 

# Function that fits candidate models by region
fit_candidate_models <- function(data_CGFS, meshes_configurations) {
  
  
  families_by_response <- list(
    ScanWeightKg = list(
      tweedie = sdmTMB::tweedie(link = "log"),
      gamma = Gamma(link = "log"),
      lognormal = lognormal(link = "log")
    )
  )
  
  
  fitted_candidate_models <- list()
  
  
  # Log-transformed swept area used as offset for total biomass
  # No offset for density (already standardized per km²)
  offset_by_response <- list(ScanWeightKg = log(data_CGFS$nb_ham * 0.246))
  
  for (mesh_name in names(meshes_configurations)) {
    
    current_mesh <- meshes_configurations[[mesh_name]]$mesh
    
    fitted_candidate_models[[mesh_name]] <- list()
    
  # ------------------------------------------------------------------------------#
  # Loop over response variables
  # ------------------------------------------------------------------------------#
  for (response in responses) {
    
    # Construct the model formula for the current response
    model_formula <- as.formula(
      paste(response, "~", fixed_effect)
    )
    
    # Select the candidate families corresponding to the current response
    families_chosen <- families_by_response[[response]]
    
    # Select the offset argument for the current response (if any)
    offset_arg <- offset_by_response[[response]]
    
    models_fitted <- list()              # Models that fitted successfully 
    models_fit_failed <- list()          # Models that failed during model fitting
    
    
    # ----------------------------------------------------------------------------#
    # Loop over candidate families
    # ----------------------------------------------------------------------------#
    
    for (family_name in names(families_chosen)) {
      
      cat("Fitting:", mesh_name, "-", response, "-", family_name, "\n")
      
      current_candidate_model <- tryCatch({
        
        # List of arguments passed to sdmTMB()
        args <- list(
          data = data_CGFS,
          formula = model_formula,
          mesh = current_mesh,
          family = families_chosen[[family_name]],
          spatial = "on",
          time = "TAAFseason",
          spatiotemporal = "IID"
        )
        
        # Add offset to the model only if it is defined
        if (!is.null(offset_arg)) args$offset <- offset_arg
        
        # Fit the model using arguments
        do.call(sdmTMB, args)
        
      }, error = function(e) e)   # Return the error object instead of stopping execution
      
      if (inherits(current_candidate_model, "error")) {
        models_fit_failed[[family_name]] <- current_candidate_model$message    # Store the error message for this family
        next
      }
      
      models_fitted[[family_name]] <- current_candidate_model
    }
    
    fitted_candidate_models[[mesh_name]][[response]] <- list(
      models_fitted = models_fitted,                # Models that fitted successfully for this response
      models_fit_failed = models_fit_failed         # Models that failed during model fitting for this response
    )
  }
  
  }
  return(fitted_candidate_models)
  
}




