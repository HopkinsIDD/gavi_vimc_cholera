#!/bin/bash
#SBATCH --job-name=rerun_eth
#SBATCH --time=10000:00:00
#SBATCH --ntasks=1
#SBATCH --mem=10G
#SBATCH --output=YEM_debugging_log/eth_rerun.txt
#SBATCH --nodelist=idmodeling2

echo "Beginning of script"
date

TAXDIR=/home/kaiyuezou/VIMC_Model/202110gavi-3/gavi_vimc_cholera
RSCRIPT=/opt/R/4.0.3/bin/Rscript

cd $TAXDIR
$RSCRIPT $TAXDIR/scripts/run_model.R -c $TAXDIR/configs/202110gavi-3/campaign-default/ETH_campaign-default_30.yml
$RSCRIPT $TAXDIR/scripts/run_model.R -c $TAXDIR/configs/202110gavi-3/no-vaccination/ETH_no-vaccination_30.yml

echo "End of script"
date
