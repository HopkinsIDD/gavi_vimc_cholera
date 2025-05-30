#!/bin/bash
#SBATCH --job-name montagu
#SBATCH --time=2-23:59:00
#SBATCH --mem=8G
#SBATCH --ntasks=1
#SBATCH --partition=defq
#SBATCH --account=aazman1

# Specify toolchain
export GCC_VERSION=9.3.0
export R_VERSION=4.0.2

# Set filepaths for:
## Local R libraries and reset
export CHOLERA_DIRECTORY=/data/aazman1/$USER/gavi-modeling/gavi_vimc_cholera
export R_LIBRARY_DIRECTORY=$HOME/rlibs/gcm/$R_VERSION/gcc/$GCC_VERSION/

## Set up modules
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


## Actually install all the 
cd $CHOLERA_DIRECTORY \
 && Rscript -e "require(remotes, lib='$R_LIBRARY_DIRECTORY'); 
                install.packages('devtools', lib='$R_LIBRARY_DIRECTORY'); 
                library(drat, lib='$R_LIBRARY_DIRECTORY'); 
                drat:::add('vimc'); 
                source('$CHOLERA_DIRECTORY/scripts/montagu_handle.R'); 
                install.packages('montagu', lib='$R_LIBRARY_DIRECTORY'); 
                library(montagu, character.only = T, lib='$R_LIBRARY_DIRECTORY'); 
                
                library(desc, lib='$R_LIBRARY_DIRECTORY'); 
                library(pkgload, lib='$R_LIBRARY_DIRECTORY'); 
                library(roxygen2, lib='$R_LIBRARY_DIRECTORY'); 
                roxygen2::roxygenise('$CHOLERA_DIRECTORY/packages/ocvImpact');
                install.packages('$CHOLERA_DIRECTORY/packages/ocvImpact', type='source', repos = NULL, lib='$R_LIBRARY_DIRECTORY'); 
                library('ocvImpact', character.only = T, lib='$R_LIBRARY_DIRECTORY')" 

echo "DONE"

