#!/bin/bash
#SBATCH --job-name setup
#SBATCH --time=6:00:00
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


## Actually install all the 
mkdir -p $R_LIBRARY_DIRECTORY \
 && cd $CHOLERA_DIRECTORY \
 && Rscript -e "options(error=quit, status = 1); 
                if(!require(package = 'remotes', character.only = T, lib='$R_LIBRARY_DIRECTORY')){ 
                    install.packages('remotes', lib='$R_LIBRARY_DIRECTORY') 
                    }; 
                require(remotes, lib='$R_LIBRARY_DIRECTORY'); 
                chooseCRANmirror(ind = 77); 
                #install.packages('netcdf', lib='$R_LIBRARY_DIRECTORY'); 
                #library('netcdf', lib='$R_LIBRARY_DIRECTORY'); 
                install_version('sf', version = '1.0.8', repos = 'http://cran.us.r-project.org', lib='$R_LIBRARY_DIRECTORY')" 

echo "DONE"

# mkdir -p $R_LIBRARY_DIRECTORY \
#  && cd $CHOLERA_DIRECTORY \
#  && Rscript -e "options(error=quit, status = 1); install.packages('remotes', lib='$R_LIBRARY_DIRECTORY')" \
#  && Rscript -e "options(error=quit, status = 1); require(remotes, lib='$R_LIBRARY_DIRECTORY'); install_version('sf', version = '1.0.8', repos = 'http://cran.us.r-project.org', lib='$R_LIBRARY_DIRECTORY', dependencies = TRUE)" \
#  && Rscript -e "options(error=quit, status = 1); require(remotes, lib='$R_LIBRARY_DIRECTORY'); install_version('GADMTools', version = '3.8-1', repos = 'http://cran.us.r-project.org', lib='$R_LIBRARY_DIRECTORY')" \
#  && Rscript -e "options(error=quit, status = 1); install.packages('drat', lib='$R_LIBRARY_DIRECTORY')" \
#  && Rscript -e "options(error=quit, status = 1); install.packages('roxygen2', lib='$R_LIBRARY_DIRECTORY')" \
#  && Rscript -e "options(error=quit, status = 1); install.packages('data.table', lib='$R_LIBRARY_DIRECTORY')" \
#  && Rscript -e "options(error=quit, status = 1); install.packages('exactextractr', lib='$R_LIBRARY_DIRECTORY')" \
#  && Rscript -e "options(error=quit, status = 1); install.packages('fasterize', lib='$R_LIBRARY_DIRECTORY')" \
#  && Rscript -e "options(error=quit, status = 1); install.packages('truncnorm', lib='$R_LIBRARY_DIRECTORY')" \
#  && Rscript -e "options(error=quit, status = 1); install.packages('MCMCglmm', lib='$R_LIBRARY_DIRECTORY')" \
#  && Rscript -e "options(error=quit, status = 1); install.packages('$CHOLERA_DIRECTORY/packages/ocvImpact', type='source', repos = NULL, lib='$R_LIBRARY_DIRECTORY')"

# echo "DONE"
