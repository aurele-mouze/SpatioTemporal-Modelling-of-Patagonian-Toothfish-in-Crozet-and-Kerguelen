
# Automatic installation, update and loading of required packages
# -----------------------------------------------------------------------------#
# -----------------------------------------------------------------------------#


install_or_update_packages <- function(packages_list){
  
  # Get names of currently installed packages
  installed <- rownames(installed.packages())
  
  # Exclude GitHub-only packages (sdmTMBextra) from CRAN installation
  cran_packages <- setdiff(packages_list, c("sdmTMBextra", "rnaturalearthhires"))
  
  # Identify CRAN packages that are not yet installed
  missing <- cran_packages[!cran_packages %in% installed]
  
  # Install missing CRAN packages 
  if(length(missing) > 0){
    install.packages(missing, dependencies = TRUE)
  }
  
  # Refresh the list of installed packages after installation
  installed <- rownames(installed.packages())
  
  # Ensure 'remotes' is installed (required for GitHub installations)
  if(!"remotes" %in% installed){
    install.packages("remotes")
  }
  
  # Refresh again after installing 'remotes'
  installed <- rownames(installed.packages())
  
  # Install sdmTMBextra from GitHub if not already installed
  if(!"sdmTMBextra" %in% installed){
    remotes::install_github("pbs-assess/sdmTMBextra", dependencies = TRUE)
  }
  if(!"rnaturalearthhires" %in% installed){
    remotes::install_github("ropensci/rnaturalearthhires")
  }
  # Update only the packages listed (not the entire library)
  update.packages(oldPkgs = packages_list, ask = FALSE)
  
  # Load all requested packages silently
  invisible(lapply(packages_list, function(package){
    suppressPackageStartupMessages(library(package, character.only = TRUE))
  }))
}


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
  "visreg"
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


nc_cache <- new.env(parent = emptyenv())