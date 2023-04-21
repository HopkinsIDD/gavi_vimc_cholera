# gavi_vimc_cholera
**VIMC cholera model - estimating impact of oral cholera vaccine (OCV)**

## General information
This repo has all the code necessary to simulate two different modeling projects: the cholera burdens modeling project that estimates vaccination campaign impacts (the GAVI-VIMC core project) on cholera controls and disease burdens in multiple African countries and some countries on other continents, and the enhanced surveillance benefits modeling project to evaluate the influences of higher cholera bacteriological confirmation capacity on vaccination campaign efficiency and cost-effectiveness (the surveillance value project/paper) that's solely focused on African countries. The preprint of this study is accessible at [medRxiv](https://www.medrxiv.org/content/10.1101/2022.11.25.22282776v1). The access to the formal paper will be posted here once it's been published.  

Currently, only the surveillance project has been enabled to run on [ARCH](https://www.arch.jhu.edu/) (will be referred to as "Rockfish" hereafter), but both projects can be run on local machines. We recommend running one project at a time because these two individual projects share folders that have the exactly same names and running both of them at the same time may cause some output files to be overwritten. However, running multiple countries simultaneously within the same project is enabled and encouraged. 

## How to use this repository 
All the information regarding how to git clone this repo (run `git lfs install` first because there are large files stored at this repo), how to set up the model, how to install relevant R packages, how to run simulations, and how to diagnose model output
is described in details on the [Wiki](https://github.com/HopkinsIDD/gavi_vimc_cholera/wiki) page of this repo. Feel free to refer them to get what you need.  

If you have any questions about how to access the restricted data or how to run the model, please contact us at the [IDDynamics](http://www.iddynamics.jhsph.edu/) group.  
  <br />


## Project 1 -- the GAVI-VIMC Core Project
After the whole model is set up following the Wiki page instructions, change the "targeting_strategy" to either "affected_pop" or "incidence" in `scripts/set_all_parameters.R` and then other parameter values that should be user-specified to generate the formal configuration files. Then, submit the Shell script in `sbatch_scripts/gavi_vimc_project/gen_incid/` to generate country-level incidence rate raster files. Next, submit a Shell script that's applicable in the `sbatch_scripts/gavi_vimc_project/` to run the formal simulation. Finally, you'll be able to find all the necessary model outputs inside both `output_raw` folder and `output_final` folder. 

## Project 2 -- the Surveillance Value Study
The process to run the surveillance project is very similar to the steps described above but has some key differences: all the configuration files have already been generated and stored in `configs/202302_survms/` to ensure the replicability of the study results; if you want to make new configurations files, the "targeting_strategy" parameter can only be "threshold_unconstrained" for the model to run a set of different procedures for the study; the Shell scripts are stored within `/sbatch_scripts/surveillance_project/`; upon the completion of the simulations, all the intermediate files (if you choose to save them for future references, which is governed by the "save_final_output_raster" parameter in the configuration) and final model outputs will be saved under the `output_raw` folder, the `output_final` folder however, will still be empty. 
  <br />
  <br />


## Folder Index 
The following is a brief introduction of different folders and files on this repo to guide you through.  
  * configs: auto-generated configuration files containing parameters for both modeling projects 
  * deliverables: methodology description documentation and model comparison documentation for the GAVI-VIMC core project
  * diagnostics: the .Rmd and .html files that help visualize model output for both modeling projects 
  * input_data: public data for both models' simulations
  * montagu: restricted data from Montagu (empty until downloading data from Montagu API onto your local machine) for both projects
  * output_final: for final summarized model output (empty until simulation processes)
  * output_raw: directly generated output files from model during simulation (empty until simulation processes)
  * packages: the project-specific R package containing all functions that work within both models
  * sbatch_scripts: the Shell scripts used to submit jobs to set up the surveillance model on Rockfish, generate incidence rate raster files as model input, simulate both models, and diagnose the surveillance model  
  * scripts: R scripts used to initialize and run both models, as well as to process raw model output to prepare tables and figures for the surveillance project
  * summarize_outputs: code and data that can be used to summarize model output to generate tables and figures for the publication of the surveillance study and the actual tables and figures included in the manuscript 
  * write_up: documents generated and data used before submitting the formal model products to VIMC




