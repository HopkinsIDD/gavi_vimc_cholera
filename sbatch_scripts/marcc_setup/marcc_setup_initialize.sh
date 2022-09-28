#!/bin/bash
#SBATCH --job-name setup
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
mkdir -p $R_LIBRARY_DIRECTORY
cd $HOME/rlibs/gcm
rm -r *

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
mkdir -p $R_LIBRARY_DIRECTORY \
 && cd $CHOLERA_DIRECTORY \
 && Rscript -e "options(error=quit, status = 1); 
                install.packages('remotes', lib='$R_LIBRARY_DIRECTORY'); 
                require(remotes, lib='$R_LIBRARY_DIRECTORY'); 
                install_version('sf', version = '1.0.8', repos = 'http://cran.us.r-project.org', lib='$R_LIBRARY_DIRECTORY', dependencies = TRUE); 
                install_version('GADMTools', version = '3.9.1', repos = 'http://cran.us.r-project.org', lib='$R_LIBRARY_DIRECTORY');
                install_version('exactextractr', version = '0.9.0', repos = 'http://cran.us.r-project.org', lib='$R_LIBRARY_DIRECTORY');  
                install_version('raster', version = '3.4.13', repos = 'http://cran.us.r-project.org', lib='$R_LIBRARY_DIRECTORY'); 
                install_version('Rcpp', version = '1.0.9', repos = 'http://cran.us.r-project.org', lib='$R_LIBRARY_DIRECTORY'); 
                install_version('terra', version = '1.4.22', repos = 'http://cran.us.r-project.org', lib='$R_LIBRARY_DIRECTORY'); 

                install.packages('drat', lib='$R_LIBRARY_DIRECTORY'); 
                install.packages('roxygen2', lib='$R_LIBRARY_DIRECTORY'); 
                install.packages('data.table', lib='$R_LIBRARY_DIRECTORY'); 
                install.packages('fasterize', lib='$R_LIBRARY_DIRECTORY'); 
                install.packages('truncnorm', lib='$R_LIBRARY_DIRECTORY'); 
                install.packages('MCMCglmm', lib='$R_LIBRARY_DIRECTORY'); 
                install.packages('codetools', lib='$R_LIBRARY_DIRECTORY'); 
                library(drat, lib='$R_LIBRARY_DIRECTORY'); 
                drat:::add('vimc'); 
                install.packages('montagu', lib='$R_LIBRARY_DIRECTORY'); 
                library(montagu, lib='$R_LIBRARY_DIRECTORY'); 
                
                library(desc, lib='$R_LIBRARY_DIRECTORY'); 
                library(pkgload, lib='$R_LIBRARY_DIRECTORY'); 
                roxygen2::roxygenise('$CHOLERA_DIRECTORY/packages/ocvImpact');
                install.packages('$CHOLERA_DIRECTORY/packages/ocvImpact', type='source', repos = NULL, lib='$R_LIBRARY_DIRECTORY'); 
                library(ocvImpact, lib='$R_LIBRARY_DIRECTORY')" 

echo "DONE"
