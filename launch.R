rm(list=ls())

config_name <- "config_template.R"
main <- "main.R"
library(here)

source(here('config', config_name))
source(here('main.R'))
