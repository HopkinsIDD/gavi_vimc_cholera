#!/bin/bash
#SBATCH --job-name=cod_glb
#SBATCH --time=3-00:00:00
#SBATCH --mem=30G
#SBATCH --partition=defq
#SBATCH --account=aazman1

export GCC_VERSION=9.3.0
export R_VERSION=4.0.2

module purge

ml gcc/$GCC_VERSION
ml openmpi
ml gdal
ml r/$R_VERSION
ml udunits
ml proj
ml libjpeg
ml sqlite
ml geos
ml libpng
ml curl

ml r-magrittr
ml r-optparse
ml r-yaml
ml r-rprojroot
ml r-purrr
ml r-jsonlite
ml r-dplyr
ml r-tidyverse
ml r-stringi

echo "Beginning of script"
date

export CHOLERA_DIRECTORY=/data/aazman1/$USER/gavi-modeling/gavi_vimc_cholera/
export CHOLERA_CONFIG_DIRECTORY=/data/aazman1/$USER/gavi-modeling/gavi_vimc_cholera/configs/202110gavi-3/campaign-default/global-estimate/2e-04
export R_LIBRARY_DIRECTORY=$HOME/rlibs/gcm/$R_VERSION/gcc/$GCC_VERSION/

cd $CHOLERA_DIRECTORY
Rscript -e "Sys.setenv(R_LIBRARY_DIRECTORY='$R_LIBRARY_DIRECTORY'); 
            Sys.setenv(RUN_ON_MARCC=TRUE); 
            Sys.setenv(CHOLERA_CONFIG='configs/202110gavi-3/campaign-default/global-estimate/1e-04/COD_campaign-default_global-estimate_50.yml'); 
            source('scripts/run_model.R')" || exit 1

echo "End of script"
date
