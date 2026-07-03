# ==============================================================================#
# Repeated random cross-validation on convergent sdmTMB models
# Global R² computed from all concatenated CV predictions
# ==============================================================================#

library(dplyr)
library(purrr)
library(tibble)

repeated_random_CV <- function(converged_models,
                               meshes,
                               data_CGFS,
                               n_repeats = 10,
                               k_folds = 10,
                               seed = 123) {
  
  options(future.globals.maxSize = Inf)
  
  data_CGFS <- data_CGFS %>%
    mutate(logSweptAreaKm2 = log(data_CGFS$nb_ham * 0.246))
  
  cv_outputs <- list()
  cv_scores <- list()
  cv_predictions <- list()
  
  models_by_mesh <- converged_models
  
  for (mesh_name in names(models_by_mesh)) {
    
    models_by_response <- models_by_mesh[[mesh_name]]
    
    cv_outputs[[mesh_name]] <- list()
    cv_scores[[mesh_name]] <- list()
    cv_predictions[[mesh_name]] <- list()
    
    for (response_name in names(models_by_response)) {
      
      models_by_family <- models_by_response[[response_name]]
      
      cv_outputs[[mesh_name]][[response_name]] <- list()
      cv_scores[[mesh_name]][[response_name]] <- list()
      cv_predictions[[mesh_name]][[response_name]] <- list()
      
      for (family_name in names(models_by_family)) {
        
        fitted_model <- models_by_family[[family_name]]
        
        model_formula <- formula(fitted_model)
        distr_family <- fitted_model$family
        mesh <- meshes[[mesh_name]]$mesh
        
        cv_outputs[[mesh_name]][[response_name]][[family_name]] <- list()
        cv_scores[[mesh_name]][[response_name]][[family_name]] <- list()
        cv_predictions[[mesh_name]][[response_name]][[family_name]] <- list()
        
        for (rep_id in seq_len(n_repeats)) {
          
          set.seed(seed + rep_id)
          
          fold_ids <- sample(
            rep(seq_len(k_folds), length.out = nrow(data_CGFS))
          )
          
          sdmTMB_cv_args <- list(
            formula = model_formula,
            data = data_CGFS,
            mesh = mesh,
            family = distr_family,
            spatial = "on",
            time = "TAAFseason",
            spatiotemporal = "IID",
            fold_ids = fold_ids,
            k_folds = k_folds
          )
          
          if (response_name == "totalWeightKg") {
            sdmTMB_cv_args$offset <- "logSweptAreaKm2"
          }
          
          cv_fit <- tryCatch(
            do.call(sdmTMB::sdmTMB_cv, sdmTMB_cv_args),
            error = function(e) {
              message("\nCV ERROR:")
              message("mesh: ", mesh_name)
              message("response: ", response_name)
              message("family: ", family_name)
              message("error: ", e$message)
              return(NULL)
            }
          )
          
          cv_outputs[[mesh_name]][[response_name]][[family_name]][[rep_id]] <- cv_fit
          
          if (inherits(cv_fit, "error")) {
            
            cv_scores[[mesh_name]][[response_name]][[family_name]][[rep_id]] <- tibble(
              mesh = mesh_name,
              response = response_name,
              family = family_name,
              repeat_id = rep_id,
              sum_loglik = NA_real_,
              rmse = NA_real_,
              mae = NA_real_,
              n = NA_integer_,
              n_folds_tot = NA_integer_,
              n_converged = NA_integer_,
              percent_converged = NA_real_,
              all_converged = FALSE,
              error = cv_fit$message
            )
            
            next
          }
          
          obs <- cv_fit$data[[response_name]]
          pred <- cv_fit$data$cv_predicted
          
          
          if (is.null(obs) || is.null(pred)) next
          
          
          ok <- is.finite(obs) & is.finite(pred)
          
          cv_predictions[[mesh_name]][[response_name]][[family_name]][[rep_id]] <- tibble(
            mesh = mesh_name,
            response = response_name,
            family = family_name,
            repeat_id = rep_id,
            row_id = seq_along(obs),
            observed = obs,
            predicted = pred
          ) %>%
            filter(is.finite(observed), is.finite(predicted))
          
          cv_scores[[mesh_name]][[response_name]][[family_name]][[rep_id]] <- tibble(
            mesh = mesh_name,
            response = response_name,
            family = family_name,
            repeat_id = rep_id,
            sum_loglik = cv_fit$sum_loglik,
            rmse = sqrt(mean((obs[ok] - pred[ok])^2, na.rm = TRUE)),
            mae = mean(abs(obs[ok] - pred[ok]), na.rm = TRUE),
            n = sum(ok),
            n_folds_tot = length(cv_fit$pdHess),
            n_converged = sum(cv_fit$pdHess),
            percent_converged = 100 * mean(cv_fit$pdHess),
            all_converged = isTRUE(cv_fit$converged),
            error = NA_character_
          )
        }
      }
    }
  }
  
  cv_scores_df <- map_dfr(cv_scores, function(mesh_list) {
    map_dfr(mesh_list, function(response_list) {
      map_dfr(response_list, function(family_list) {
        bind_rows(family_list)
      })
    })
  })
  
  cv_predictions_df <- map_dfr(cv_predictions, function(mesh_list) {
    map_dfr(mesh_list, function(response_list) {
      map_dfr(response_list, function(family_list) {
        if (length(family_list) == 0) return(NULL)
        bind_rows(family_list)
      })
    })
  })
  if (nrow(cv_predictions_df) == 0) {
    stop("cv_predictions_df est vide : aucun CV réussi ou aucune prédiction générée.")
  }
  cv_scores_summary <- cv_predictions_df %>%
    group_by(mesh, response, family) %>%
    summarise(
      n_predictions = n(),
      n_repeats = n_distinct(repeat_id),
      
      sum_loglik_mean = mean(
        cv_scores_df$sum_loglik[
            cv_scores_df$mesh == first(mesh) &
            cv_scores_df$response == first(response) &
            cv_scores_df$family == first(family)
        ],
        na.rm = TRUE
      ),
      
      sum_loglik_sd = sd(
        cv_scores_df$sum_loglik[
            cv_scores_df$mesh == first(mesh) &
            cv_scores_df$response == first(response) &
            cv_scores_df$family == first(family)
        ],
        na.rm = TRUE
      ),
      
      RMSE_global = sqrt(mean((observed - predicted)^2, na.rm = TRUE)),
      MAE_global = mean(abs(observed - predicted), na.rm = TRUE),
      R2_global = cor(observed, predicted, use = "complete.obs")^2,
      
      
      mean_percent_converged = mean(
        cv_scores_df$percent_converged[
            cv_scores_df$mesh == first(mesh) &
            cv_scores_df$response == first(response) &
            cv_scores_df$family == first(family)
        ],
        na.rm = TRUE
      ),
      
      .groups = "drop"
    )
  
  list(
    cv_outputs = cv_outputs,
    cv_scores = cv_scores_df,
    cv_predictions = cv_predictions_df,
    cv_scores_summary = cv_scores_summary
  )
}
