# gavi_vimc_cholera
**VIMC cholera model - estimating impact of oral cholera vaccine (OCV)**

## General information
This repo has all the code necessary to simulate two different modeling projects: the cholera burdens modeling project that estimates vaccination campaign impacts (the GAVI-VIMC core project) on cholera controls and disease burdens in multiple African countries and some countries on other continents, and the enhanced surveillance benefits modeling project to evaluate the influences of higher cholera bacteriological confirmation capacity on vaccination campaign efficiency and cost-effectiveness (the surveillance value project/paper) that's solely focused on African countries. The preprint of this study is accessible at [medRxiv](https://www.medrxiv.org/content/10.1101/2022.11.25.22282776v1). The access to the formal paper will be posted here once it's been published.  

We recommend running one project at a time because these two individual projects share folders that have the exactly same names and running both of them at the same time may cause some output files to be overwritten. However, running multiple countries simultaneously within the same project is enabled and encouraged. 

## How to use this repository 
All the information regarding how to git clone this repo (run `git lfs install` first because there are large files stored at this repo), how to set up the model, how to install relevant R packages, how to run simulations, and how to diagnose model output
is described in details on the [Wiki](https://github.com/HopkinsIDD/gavi_vimc_cholera/wiki) page of this repo. Feel free to refer them to get what you need.  

If you have any questions about how to access the restricted data or how to run the model, please contact us at the [IDDynamics](http://www.iddynamics.jhsph.edu/) group.  
  <br />


## Project 1 -- the GAVI-VIMC Core Project

## Project 2 -- the Surveillance Value Study
  <br />


## Folder/File Index 
The following is a brief introduction of different folders and files on this repo to guide you through.  
  * configs: auto-generated configuration files containing parameters for both modeling projects 
  * deliverables: methodology description documentation and model comparison documentation for the GAVI-VIMC core project
  * diagnostics: the .Rmd and .html files that help visualize model output for both modeling projects 
  * input_data: public data for both models' simulations
  * montagu: restricted data from Montagu (empty until downloading data from Montagu API onto your local machine) for both projects
  * output_final: for final summarized model output (empty until simulation processes)
  * output_raw: directly generated output files from model during simulation (empty until simulation processes)
  * packages: the project-specific R package containing all functions that work within both models
  * sbatch_scripts: the Shell scripts used to submit jobs to set up the surveillance model on Rockfish, generate incidence rate raster files as model input, and simulate both models 
  * scripts: R scripts used to initialize and run both models, as well as to process raw model output to prepare tables and figures for the surveillance project
  * summarize_outputs: scripts to summarize model output to generate visuals for publication and visuals included in the manuscript 
  * write_up: documents and tables used when submitting model products to VIMC
  * diagnostic_report.Rmd: the diagnostic report RMarkdown file for the surveillance project




