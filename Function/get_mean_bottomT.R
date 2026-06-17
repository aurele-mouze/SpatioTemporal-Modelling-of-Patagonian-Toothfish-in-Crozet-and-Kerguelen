
get_mean_bottomT <- function(lon_sta, lat_sta,
                             date_depot, date_remontee,
                             dossier_nc){
  
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
        "bottomT_daily_%s_%s.nc",
        format(d, "%Y"),
        format(d, "%m")
      )
    )
    
    if(!file.exists(fichier)){
      valeurs[k] <- NA
      next
    }
    
    dat <- get_nc_data(fichier)
    
    i <- which.min(abs(dat$lon - lon_sta))
    j <- which.min(abs(dat$lat - lat_sta))
    
    jour_mois <- as.integer(format(d, "%d"))
    
    valeurs[k] <- dat$temp[i, j, jour_mois]
  }
  
  mean(valeurs, na.rm = TRUE)
}