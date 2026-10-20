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
