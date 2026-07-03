# ============================================================================== #
#### COMPUTE AIC WEIGHTS FROM VALID MODELS ####
# ============================================================================== #


compute_cAIC_weights <- function(converged_models) {
  
  models_by_mesh <- converged_models
  
  AIC_list <- list()
  
  for (mesh_name in names(models_by_mesh)) {
    
    models_by_response <- models_by_mesh[[mesh_name]]
    
    for (response_name in names(models_by_response)) {
      
      models_by_family <- models_by_response[[response_name]]
    
      message(
        "Computing cAIC: ",
        mesh_name, " | ",
        response_name
      )
      
    # Compute conditional AIC for each family
    AIC_values <- vapply(models_by_family,
      function(current_model) {
        val <- sdmTMB::cAIC(current_model, what = "cAIC")
        as.numeric(val)
      },
      numeric(1))
    
    # Delta AIC
    deltaAIC <- AIC_values - min(AIC_values)
    
    # AIC weights : difference between the AIC for model i and the model with the minimum AIC
    relative_likelihood <- exp(-0.5 * deltaAIC)
    weights <- relative_likelihood / sum(relative_likelihood)
    
    AIC_df <- data.frame(
      mesh = mesh_name,
      response = response_name,
      family = names(AIC_values),
      AIC = as.numeric(AIC_values),
      deltaAIC = as.numeric(deltaAIC),
      weight = as.numeric(weights)
    )
    
    AIC_df <- AIC_df[order(AIC_df$AIC), ]
    
    AIC_list[[paste(mesh_name, response_name, sep = "_")]] <- AIC_df
    
    }
  }
  
  AIC_weights_df <- do.call(rbind, AIC_list)
  rownames(AIC_weights_df) <- NULL
  
  return(AIC_weights_df)
}



# AIC_by_region <- list()
# 
# if (isTRUE(West_English_Channel)) {
#   AIC_by_region$west <- compute_cAIC_weights(converged_models, "west")
# }
# 
# if (isTRUE(East_English_Channel)) {
#   AIC_by_region$east <- compute_cAIC_weights(converged_models, "east")
# }



