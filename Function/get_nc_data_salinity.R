#lit un fichier NetCDF contenant les données de salinité 
#et extrait les coordonnées géographiques ainsi que les profondeurs associées.

get_nc_data_salinity <- function(fichier){
  
  nc <- nc_open(fichier)
  
  lon <- ncvar_get(nc, "longitude")
  lat <- ncvar_get(nc, "latitude")
  depth <- ncvar_get(nc, "depth")
  salinity <- ncvar_get(nc, "so")
  
  nc_close(nc)
  
  list(
    lon = lon,
    lat = lat,
    depth = depth,
    salinity = salinity
  )
}
