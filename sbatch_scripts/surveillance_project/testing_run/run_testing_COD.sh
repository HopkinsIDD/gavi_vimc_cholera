#!/bin/bash
#SBATCH --job-name=COD_test
#SBATCH --output=log/COD_testing_run.txt
#SBATCH --time=10000:00:00
#SBATCH --ntasks=1
#SBATCH --mem=10G
#SBATCH --nodelist=idmodeling2

echo "Beginning of script"
date

TAXDIR=/home/kaiyuezou/VIMC_Model/surveillance_project/gavi_vimc_cholera
CONFIGDIR=configs/202110gavi-3
RSCRIPT=/opt/R/4.0.3/bin/Rscript

cd $TAXDIR
$RSCRIPT $TAXDIR/scripts/run_model.R -c $TAXDIR/$CONFIGDIR/campaign-default/district-estimate/COD_campaign-default_district-estimate_50.yml 
$RSCRIPT $TAXDIR/scripts/run_model.R -c $TAXDIR/$CONFIGDIR/campaign-default/global-estimate/COD_campaign-default_global-estimate_50.yml
$RSCRIPT $TAXDIR/scripts/run_model.R -c $TAXDIR/$CONFIGDIR/campaign-default/no-estimate/COD_campaign-default_no-estimate_50.yml
$RSCRIPT $TAXDIR/scripts/run_model.R -c $TAXDIR/$CONFIGDIR/no-vaccination/district-estimate/COD_no-vaccination_district-estimate_50.yml || exit 1

echo "End of script"
date
