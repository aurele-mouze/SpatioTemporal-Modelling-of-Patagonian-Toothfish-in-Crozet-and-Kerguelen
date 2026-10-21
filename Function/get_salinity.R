#extrait les salinités journalières du fond marin à la station la plus proche des coordonnées
#fournies, à partir de fichiers NetCDF. 
#Elle calcule ensuite la salinité moyenne sur la période comprise 
#entre le dépôt et la remontée de la ligne de pêche.

get_mean_salinity <- function(lon_sta, lat_sta,
                              depth_sta,
                              date_depot, date_remontee,
                              dossier_nc){
  
  # Ignorer les stations sans coordonnées ou sans profondeur
  if(is.na(lon_sta) || is.na(lat_sta) || is.na(depth_sta)){
    return(NA_real_)
  }
  
  jours <- seq(
    as.Date(date_depot),
    as.Date(date_remontee),
    by = "day"
  )
  
  valeurs <- numeric(length(jours))
  
  for(k in seq_along(jours)){
    
    d <- jours[k]
    
    fichier <- file.path(
      dossier_nc,
      sprintf(
        "so_daily_%s_%s.nc",
        format(d, "%Y"),
        format(d, "%m")
      )
    )
    
    if(!file.exists(fichier)){
      valeurs[k] <- NA
      next
    }
    
    dat <- get_nc_data_salinity(fichier)
    
    i <- which.min(abs(dat$lon - lon_sta))
    j <- which.min(abs(dat$lat - lat_sta))
    z <- which.min(abs(dat$depth - depth_sta))
    
    jour_mois <- as.integer(format(d, "%d"))
    
    valeurs[k] <- dat$salinity[i, j, z, jour_mois]
  }
  
  mean(valeurs, na.rm = TRUE)
}
