# ============================================================================== #
#### SANITY FILTER FITTED MODELS ####
# ============================================================================== #

# This function takes the output of fit_candidate_models() and applies 
# sdmTMB::sanity() to all successfully fitted models.

sanity_filter_models <- function(fit_outputs) {
  
  sanity_results <- list()
  
  # --------------------------------------------------------------------------- #
  # Loop over mesh configurations
  # --------------------------------------------------------------------------- #
  
  for (mesh_name in names(fit_outputs)) {
    
    sanity_results[[mesh_name]] <- list()
    
    # ------------------------------------------------------------------------- #
    # Loop over response variables
    # ------------------------------------------------------------------------- #
    
    for (response in names(fit_outputs[[mesh_name]])) {
      
      models_fitted <- fit_outputs[[mesh_name]][[response]]$models_fitted
      models_fit_failed <- fit_outputs[[mesh_name]][[response]]$models_fit_failed
      
      models_valid <- list()
      models_sanity_failed <- list()
      

    # ------------------------------------------------------------------------- #
    # Loop over candidate families for this response
    # ------------------------------------------------------------------------- #
    for (family_name in names(models_fitted)) {
      
      # Print progress message
      cat("Sanity:", mesh_name, "-", response, "-", family_name, "\n")
      
      # Extract the fitted model
      current_checked_model <- models_fitted[[family_name]]
      
      # Run sdmTMB sanity checks on the fitted model
      sanity_res <- tryCatch(sanity(current_checked_model), error = function(e) e)
      
      # If sanity() fails or returns all_ok = FALSE, store in sanity_failed
      if (inherits(sanity_res, "error") || !isTRUE(sanity_res$all_ok)) {
        models_sanity_failed[[family_name]] <- sanity_res
      } else {
        # Otherwise, keep the model as valid
        models_valid[[family_name]] <- current_checked_model
      }
    }
    
    # Store results for the current response variable
    sanity_results[[mesh_name]][[response]] <- list(
      models_valid  = models_valid,                    # Models that fitted successfully and passed sanity checks for this response
      models_fit_failed  = models_fit_failed,          # Models that failed during model fitting for this response
      models_sanity_failed  = models_sanity_failed     # Models that failed sanity diagnostics for this responsemodels_sanity_failed  = models_sanity_failed     # Models that failed sanity diagnostics for this response
    )
  }
  
  }
  # Return sanity-filtered models for all responses
  return(sanity_results)
}


# ----------------------------------------------------------------------------- #
# sanity filtering by region
# ----------------------------------------------------------------------------- #

# sanity_by_region <- list()
# 
# # Apply sanity checks to East English Channel models
# if (isTRUE(East_English_Channel)) {
#   sanity_by_region$east <- sanity_filter_models(fitted_models_by_region$east)
# }
# 
# # Apply sanity checks to West English Channel models
# if (isTRUE(West_English_Channel)) {
#   sanity_by_region$west <- sanity_filter_models(fitted_models_by_region$west)
# }








