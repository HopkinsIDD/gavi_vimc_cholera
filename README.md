# gavi_vimc_cholera
**VIMC cholera model - estimating impact of oral cholera vaccine (OCV)**

This repo has all the code necessary to simulate the spatial cholera model to estimate vaccination campaign impacts (the VIMC project) or 
to evaluate the influences of higher cholera testing capacity on cholera control (the surveillance project).  

Please note to run the model, some data files need to be downloaded using the [VIMC Montagu](https://montagu.vaccineimpact.org/) API and
the access is constrained.  

All the information regarding how to set up the model, how to install relevant R packages, how to run simulations, and how to diagnose model output
is described in details on the [Wiki](https://github.com/HopkinsIDD/gavi_vimc_cholera/wiki) page of this repo. Feel free to refer them to get what you need.  

If you have any questions about how to access the restricted data or how to run the model, please contact us at the [IDDynamics](http://www.iddynamics.jhsph.edu/) group.  
  <br />
  <br />
The following is a brief introduction of different folders and files on this repo to guide you through.  
  * configs: auto-generated configuration files containing parameters for the model
  * deliverables: study methodology description documents or publications
  * diagnostics: `.html` files that visualize model output
  * input_data: public data for model simulation
  * montagu: restricted data from Montagu (empty until downloading data from Montagu API on your local machine)
  * output_final: for final summarized model output
  * output_raw: directly generated output files from model during simulation
  * packages: the R package containing functions that work within the model
  * sbatch_scripts: Shell scripts to submit jobs
  * scripts: R scripts to initialize and run the model
  * summarize_outputs: scripts to summarize model output to generate visuals for publication and visuals included in the manuscript 
  * write_up: documents and tables used when submitting model products to VIMC
  * diagnostic_report.Rmd: the diagnostic report RMarkdown file for the surveillance project



