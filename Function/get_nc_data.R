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