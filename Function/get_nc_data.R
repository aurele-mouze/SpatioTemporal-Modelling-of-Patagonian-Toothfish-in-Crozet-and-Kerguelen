#Lit les températures du fond marin (bottomT) et les coordonnées géographiques d’un fichier NetCDF. 
#Elle met les données en cache dans nc_cache pour éviter de relire 
#plusieurs fois les mêmes fichiers

get_nc_data <- function(fichier){
  
  if(!exists(fichier, envir = nc_cache)){
    
    nc <- nc_open(fichier)
    
    assign(
      fichier,
      list(
        temp = ncvar_get(nc, "bottomT"),
        lon  = ncvar_get(nc, "longitude"),
        lat  = ncvar_get(nc, "latitude")
      ),
      envir = nc_cache
    )
    
    nc_close(nc)
  }
  
  get(fichier, envir = nc_cache)
}