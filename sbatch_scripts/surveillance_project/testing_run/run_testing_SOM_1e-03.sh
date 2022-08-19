#!/bin/bash
#SBATCH --job-name=SOM_test
#SBATCH --output=log/SOM_testing_run_1e-03.txt
#SBATCH --time=10000:00:00
#SBATCH --ntasks=1
#SBATCH --mem=15G
#SBATCH --nodelist=idmodeling2

echo "Beginning of script"
date

TAXDIR=/home/kaiyuezou/VIMC_Model/surveillance_project/gavi_vimc_cholera
CONFIGDIR=configs/202110gavi-3
RSCRIPT=/opt/R/4.0.3/bin/Rscript

cd $TAXDIR
$RSCRIPT $TAXDIR/scripts/run_model.R -c $TAXDIR/$CONFIGDIR/campaign-default/district-estimate/0.001/SOM_campaign-default_district-estimate_50.yml
$RSCRIPT $TAXDIR/scripts/run_model.R -c $TAXDIR/$CONFIGDIR/campaign-default/global-estimate/0.001/SOM_campaign-default_global-estimate_50.yml
$RSCRIPT $TAXDIR/scripts/run_model.R -c $TAXDIR/$CONFIGDIR/campaign-default/no-estimate/0.001/SOM_campaign-default_no-estimate_50.yml
$RSCRIPT $TAXDIR/scripts/run_model.R -c $TAXDIR/$CONFIGDIR/no-vaccination/district-estimate/0.001/SOM_no-vaccination_district-estimate_50.yml || exit 1

echo "End of script"
date
