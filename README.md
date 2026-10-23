# stage_AM


stage Aurèle Mouze

# Data Formatting, Spatial Modelling and Analysis of POKER Data
 
This repository contains all scripts required to **format**, **explore**, and **analyse** data from POKER surveys in the Kerguelen and Crozet EEZs. The workflow includes data preparation, spatial grid creation, mesh and model testing, model fitting, and prediction generation using `sdmTMB`. It also produces **spatial distribution outputs**, **model effect plots**, and figures for *Dissostichus eleginoides*.
 
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
  Tests different mesh configurations and model specifications.
 
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