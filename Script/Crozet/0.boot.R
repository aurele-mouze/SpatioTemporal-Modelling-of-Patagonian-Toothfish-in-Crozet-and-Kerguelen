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
"lubridate"
)

lapply(libr, require, character.only = TRUE)

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