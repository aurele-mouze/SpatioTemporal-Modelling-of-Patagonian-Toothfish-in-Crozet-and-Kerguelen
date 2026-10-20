
source("Function/install_or_update_packages.R")


packages_list <- c(
  "ape",
  "DHARMa",
  "dplyr",
  "ggOceanMaps",
  "ggpubr",
  "ggspatial",
  "here",
  "inlabru",
  "fmesher",
  "ncdf4",
  "pak",
  "patchwork",
  "raster",
  "readr",
  "readxl",
  "remotes",
  "rmdformats",
  "rnaturalearth",
  "rnaturalearthhires",
  "scales",
  "sdmTMB",
  "sdmTMBextra",     # ! Load sdmTMBextra after sdmTMB
  "sf",
  "terra",
  "tidyverse",
  "tweedie",
  "viridis"
)

install_or_update_packages(packages_list)



libr <- c(
"tidyr",
"dplyr",
"terra",
"fmesher",
"sdmTMB",
"sdmTMBextra",
"here",
"ggplot2",
"raster",
"sf",
"ggpubr",
"POEPAmap",
"ggeffects",
"tidyterra",
"cowplot",
"lubridate",
"readr",
"ncdf4",
"visreg",
"igraph"
)

install_or_update_packages(libr)

#Création des répertoires

dir.create(here("data/", output_plot),
           showWarnings = F,
           recursive = T)

dir.create(here("data/", output_df),
           showWarnings = F,
           recursive = T)

dir.create(here("data/", output_carte),
           showWarnings = F,
           recursive = T)


source("Function/choose_mesh.R")
source("Function/filtre_flag_dynamique.R")
source("Function/barrier_mesh.R")
source("Function/predict_in_chunks.R")
source("Function/pseudo_R2.R")
source("Function/get_mean_bottomT.R")
source("Function/get_nc_data.R")
source("Function/get_salinity.R")

nc_cache <- new.env(parent = emptyenv())
