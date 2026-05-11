# stage_AM

stage Aurèle Mouze

# Data Formatting and Exploration of POKER Data
 
This repository contains all scripts required to **download**, **format**, and **explore** data from **POKER surveys**. These steps are a prerequisite for generating inputs for the **POKER SDM**, as well as producing **spatial distribution outputs** and a **biomass index** of *Dissostichus eleginoides* in the Kerguelen EEZ.
 
---
 
## Repository structure and workflow
 
The workflow is driven by a `launch` script, which:
 
1. calls a `config` file containing all parameters required to run the analysis, and  
2. calls `main.R`, which orchestrates the full workflow.
 
Templates for both the `launch` and `config` files are provided.
 
All processing steps are organised into numbered subfolders within the `Projects` directory:
 
- `0.Boot`  
  Loads required packages and utility functions, and creates the output folder structure.
 
- `1.Extract_data`  
  Scripts to extract the required data frames from PECHEKER. Outputs are saved in `Data`.  
  
 
- `2.Prepare_Grid`  
  Scripts to create the grid. Formatted outputs are saved in `Data`.  
 
 
- `4. to ....`  
  Scripts to generate sdmTMB model and analysis.
 
---
 
## Quick start
 
1. **Prepare required data**  
   Ensure all required input files (see below) are available in the `Data` folder.
 
2. **Configure the analysis**  
   Copy the configuration template and edit parameters as needed:
   ```r
   config_template.R → config.R
