# SpatioTemporal-Modelling-of-Patagonian-Toothfish-in-Crozet-and-Kerguelen

This project focuses on the spatio-temporal modelling of Patagonian toothfish (Dissostichus eleginoides) populations in the Kerguelen and Crozet EEZs (French Southern and Antarctic Lands).

The work was conducted under the supervision of Félix Massiot-Granier, Clara Péron and Fanny Ouzoulias, from the UMR BOREA laboratory at the French National Museum of Natural History (Muséum National d'Histoire Naturelle, MNHN), Paris.

Using commercial fishing data, the objective is to investigate and model the spatial and temporal distribution of Patagonian toothfish in both study areas. The analyses are carried out using sdmTMB, a statistical framework for spatial and spatio-temporal modelling based on stochastic partial differential equations (SPDEs).

Several environmental and operational covariates are incorporated into the models to investigate the factors influencing toothfish distribution and to account for spatial and temporal variability. The overall aim is to improve our understanding of population distribution patterns and their dynamics in these sub-Antarctic area.

The model structure and the selection of explanatory variables are based on previous work carried out by Mathilde Gudin (https://github.com/ldbk/EBESCO/tree/main/SDM). These preliminary analyses provided the basis for the model specifications used in this project, which were subsequently integrated into the spatial modelling workflow.

# Data Formatting, Spatial Modelling and Analysis of Fishery Data
 
This repository contains all scripts required to **format**, **explore**, and **analyse** data from PECHEKER in the Kerguelen and Crozet EEZs. The workflow includes data preparation, spatial grid creation, mesh and model testing, model fitting, and prediction generation using `sdmTMB`. It also produces **spatial distribution outputs**, **model effect plots**, and figures for *Dissostichus eleginoides*.
 
---
 
## Repository structure and workflow
 
The workflow is managed by a `launch` script, which loads the configuration files and allows the user to select the study area (Kerguelen or Crozet). Each area has its own `main` script, which orchestrates the full workflow.
 
The repository is organised into the following directories:
 
- `Function/`  
  Contains the functions used throughout the workflow.
 
- `Data/`  
  Contains the input data and a dedicated directory for storing results.
 
- `Config/`  
  Contains the configuration files and parameters required to run the analyses.
 
- `Script/`  
  Contains the scripts for each study area, organised into separate `Kerguelen/` and `Crozet/` directories. Each area has its own `main` script and follows the same general workflow:
 
- `0.Boot/`  
  Loads the required R packages and functions.
 
- `1.Conversion_data/`  
  Formats the input data for analysis.
 
- `2.Make_grid/`  
  Creates the spatial grid.
 
- `3.Test_mesh_and_model/`  
  Tests different mesh configurations and model specifications. The initial model structure and explanatory variable selection are based on previous work by Mathilde Gudin (https://github.com/ldbk/EBESCO/tree/main/SDM).
 
- `4.Get_model/`  
  Fits the selected models and generates predictions.
 
- `5.Visreg/`  
  Produces `visreg` plots to examine the effects of model covariates.
 
- `6.Get_fig/`  
  Generates the final model figures.
 
Each study area also contains a script to generate its corresponding map.
 
---
 
## Quick start
 
1. **Prepare the input data**  
  Ensure that all required input files are available in the `Data/` directory.
 
2. **Configure the analysis**  
  Check the configuration files in `Config/` and adjust the parameters if necessary.
 
3. **Select the study area and run the workflow**  
  Use the `launch` script to select Kerguelen or Crozet and run the corresponding analysis. The associated `main` script orchestrates the processing steps.
 
---
 
## Main outputs
 
The workflow produces:
 
- Formatted datasets for model fitting.
- Spatial grids and tested mesh configurations.
- Fitted `sdmTMB` models and spatial predictions.
- `visreg` plots showing the effects of model covariates.
- Maps and figures summarising the model results.</escape>