forward_selection_fixed_effects <- function(
    data,
    response,
    mesh,
    family,
    candidate_terms,
    spatial = "on",
    time = "Year",
    spatiotemporal = "IID",
    offset = NULL,
    delta_aic = 2
){
  
  selected_terms <- character(0)
  remaining_terms <- candidate_terms
  
  history <- tibble::tibble(
    step = integer(),
    variable = character(),
    AIC = numeric()
  )
  
  ## --------------------------
  ## Modèle nul
  ## --------------------------
  
  formula_best <- as.formula(
    paste(response, "~ 1")
  )
  
  args <- list(
    formula = formula_best,
    data = data,
    mesh = mesh,
    family = family,
    spatial = spatial,
    time = time,
    spatiotemporal = spatiotemporal
  )
  
  if(!is.null(offset))
    args$offset <- offset
  
  best_model <- do.call(sdmTMB::sdmTMB, args)
  
  best_aic <- AIC(best_model)
  
  step <- 1
  
  repeat{
    
    cat("\n---------------------------\n")
    cat("Step", step, "\n")
    cat("---------------------------\n")
    
    tested_models <- list()
    tested_aic <- c()
    
    ## tester toutes les variables restantes
    
    for(term in remaining_terms){
      
      current_terms <- c(selected_terms, term)
      
      formula_try <- as.formula(
        paste(response,
              "~",
              paste(current_terms, collapse = " + "))
      )
      
      args$formula <- formula_try
      
      fit <- tryCatch(
        do.call(sdmTMB::sdmTMB, args),
        error = function(e) NULL
      )
      
      if(is.null(fit))
        next
      
      tested_models[[term]] <- fit
      tested_aic[term] <- AIC(fit)
      
      cat(term, ": AIC =", round(tested_aic[term],2), "\n")
    }
    
    ## aucun modèle n'a convergé
    
    if(length(tested_aic)==0)
      break
    
    best_term <- names(which.min(tested_aic))
    best_candidate_aic <- min(tested_aic)
    
    ## amélioration suffisante ?
    
    if(best_candidate_aic < best_aic - delta_aic){
      
      cat("==> Keep:", best_term,"\n")
      
      selected_terms <- c(selected_terms, best_term)
      remaining_terms <- setdiff(remaining_terms, best_term)
      
      best_model <- tested_models[[best_term]]
      best_aic <- best_candidate_aic
      
      history <- dplyr::bind_rows(
        history,
        tibble::tibble(
          step = step,
          variable = best_term,
          AIC = best_aic
        )
      )
      
      step <- step + 1
      
    }else{
      
      cat("No variable improves the model enough.\n")
      break
      
    }
    
  }
  
  list(
    best_model = best_model,
    selected_terms = selected_terms,
    final_formula = formula(best_model),
    history = history
  )
  
}