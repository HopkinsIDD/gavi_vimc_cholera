#!/bin/bash
#SBATCH --job-name setup
#SBATCH --time=00:60:00
#SBATCH --mem=8G
#SBATCH --ntasks=1
#SBATCH --partition=defq
#SBATCH --account=aazman1

# Specify toolchain
export GCC_VERSION=9.3.0
export R_VERSION=4.0.2

# Set filepaths for:
## Local R libraries
export CHOLERA_DIRECTORY=/data/aazman1/$USER/gavi-modeling/gavi_vimc_cholera
export R_LIBRARY_DIRECTORY=$HOME/rlibs/gcm/$R_VERSION/gcc/$GCC_VERSION/

## Set up modules
## (we need to do git things since reset removes git)
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

## Actually install cmdstanr including dependent r packages and cmdstan
Rscript -e "options(error=quit, status = 1); require(remotes, lib='$R_LIBRARY_DIRECTORY'); install_version('GADMTools', version = '3.8-1', repos = 'http://cran.us.r-project.org', lib='$R_LIBRARY_DIRECTORY')" 

echo "DONE"
