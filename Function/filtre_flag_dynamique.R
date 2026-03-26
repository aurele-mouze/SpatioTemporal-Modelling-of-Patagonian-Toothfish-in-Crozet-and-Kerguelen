#fonction d'application des filtres dynamiques

filter_flag <- function(df, 
                         id_filter   = FALSE, id_value = NULL,
                         var_depth    = FALSE, 
                         var_dist    = FALSE,
                         var_length = FALSE) {
  
  filter_df <- df
  
  # Filtre ID si demandé
  if (id_filter) {
    if (is.null(id_value)) stop("Si id_filter = TRUE, il faut spécifier id_value")
    filter_df <- filter_df %>% filter(ID == id_value)
  }
  
  # Filtre cond_prof si demandé
  if (var_depth) {
    filter_df <- filter_df %>% filter(cond_prof == 1)
  }
  
  # Filtre cond_dist si demandé
  if (var_dist) {
    filter_df <- filter_df %>% filter(cond_dist == 1)
  }
  
  if(var_length) {
    filter_df <- filter_df %>% filter(cond_taille == 1)
  }
  
  return(filter_df)
}
