rm(list=ls())

config_name <- "config_template.R"
main <- "main.R"
library(here)

zone <- readline("Choisir une zone (Kerguelen / Crozet) : ")

source(here('config', config_name))

if(zone == "Kerguelen"){
  source(here("main_Ker.R"))
  
} else if(zone == "Crozet"){
  source(here("main_Cro.R"))
  
} else {
  stop("Choix invalide : taper Kerguelen ou Crozet")
}