#!/bin/bash
#SBATCH --job-name=gen_incid
#SBATCH --time=2-23:59:00
#SBATCH --mem=10G
#SBATCH --array=0-43%10
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
ml gdal/3.5.1
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
export CHOLERA_CONFIG_DIRECTORY=/data/aazman1/$USER/gavi-modeling/gavi_vimc_cholera/configs/202110gavi-3/campaign-default/district-estimate/0.001
export R_LIBRARY_DIRECTORY=$HOME/rlibs/gcm/$R_VERSION/gcc/$GCC_VERSION/
export CONFIGNAMES=($(ls $CHOLERA_CONFIG_DIRECTORY | tr ' ' '\n'))

cd $CHOLERA_DIRECTORY
Rscript -e "Sys.setenv(R_LIBRARY_DIRECTORY='$R_LIBRARY_DIRECTORY'); 
            Sys.setenv(RUN_ON_MARCC=TRUE); 
            Sys.setenv(CHOLERA_CONFIG='$CHOLERA_CONFIG_DIRECTORY/${CONFIGNAMES[$SLURM_ARRAY_TASK_ID]}'); 
            source('scripts/run_country_incid_crop.R')" || exit 1

echo "End of script"
date