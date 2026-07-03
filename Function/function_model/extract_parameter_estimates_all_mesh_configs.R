
# ------------------------------------------------------------------------------#
#### EXTRACT PARAMETER ESTIMATES FROM CONVERGED MODELS ####
# ------------------------------------------------------------------------------#

extract_params_converged_models <- function(converged_models) {
  
  library(dplyr)
  library(broom.mixed)
  library(purrr)
  
  params <- list()
  
  for (mesh_name in names(converged_models)) {
    
    models_by_response <- converged_models[[mesh_name]]
    
    for (response in names(models_by_response)) {
      
      models_by_family <- models_by_response[[response]]
      
      for (family in names(models_by_family)) {
        
        model <- models_by_family[[family]]
        
        is_delta <- grepl("^delta", family, ignore.case = TRUE)
        
        ## Fixed effects
        fixed_params <- if (is_delta) {
          
          bind_rows(
            tidy(model, conf.int = TRUE, model = 1) %>%
              mutate(component = "delta_1"),
            tidy(model, conf.int = TRUE, model = 2) %>%
              mutate(component = "delta_2")
          )
          
        } else {
          
          tidy(model, conf.int = TRUE) %>%
            mutate(component = "single")
          
        }
        
        fixed_params <- fixed_params %>%
          filter(term == "(Intercept)") %>%
          mutate(param_type = "fixed")
        
        ## Random effects
        random_params <- if (is_delta) {
          
          bind_rows(
            tidy(model,
                 effects = "ran_pars",
                 conf.int = TRUE,
                 model = 1) %>%
              mutate(component = "delta_1"),
            
            tidy(model,
                 effects = "ran_pars",
                 conf.int = TRUE,
                 model = 2) %>%
              mutate(component = "delta_2")
          )
          
        } else {
          
          tidy(model,
               effects = "ran_pars",
               conf.int = TRUE) %>%
            mutate(component = "single")
          
        }
        
        random_params <- random_params %>%
          filter(term %in% c("range", "sigma_O", "sigma_E")) %>%
          mutate(param_type = "random")
        
        params[[length(params) + 1]] <-
          
          bind_rows(fixed_params, random_params) %>%
          mutate(
            mesh = mesh_name,
            response = response,
            family = family,
            parameter = if_else(term == "(Intercept)",
                                "intercept",
                                term)
          ) %>%
          select(
            response,
            mesh,
            family,
            component,
            param_type,
            parameter,
            estimate,
            std.error,
            conf.low,
            conf.high
          )
      }
    }
  }
  
  bind_rows(params)
}



# params_west <- extract_params_converged_models(converged_models, region = "west")
# params_east <- extract_params_converged_models(converged_models, region = "east")
# 
# params_all <- bind_rows(params_east, params_west) %>%
#   mutate(region = factor(region, levels = c("west", "east")))


# params_all %>%
#   ggplot(aes(x = family, y = estimate, shape = component, color = family)) +
#   geom_point(size = 2.8, position = position_dodge(width = 0.6)) +
#   geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
#                 width = 0.2, position = position_dodge(width = 0.6)) +
#   facet_grid(parameter ~ region + response, scales = "free_y") +
#   labs(x = "Family",  y = "Estimate", shape = "Component",
#     title = "sdmTMB parameters estimates") +
#   guides(color = "none") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         legend.position = "top")

