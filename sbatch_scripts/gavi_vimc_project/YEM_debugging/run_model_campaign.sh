#!/bin/bash
#SBATCH --job-name=run_model
#SBATCH --time=10000:00:00
#SBATCH --mem=10G
#SBATCH --output=YEM_debugging_log/novac.txt
#SBATCH --nodelist=idmodeling2
#SBATCH -c 4

echo "Beginning of script"
date

TAXDIR=/home/kaiyuezou/VIMC_Model/202110gavi-3/gavi_vimc_cholera
RSCRIPT=/opt/R/4.0.3/bin/Rscript

cd $TAXDIR
$RSCRIPT $TAXDIR/scripts/run_model.R -c $TAXDIR/configs/202110gavi-3/no-vaccination/BGD_no-vaccination_30.yml
$RSCRIPT $TAXDIR/scripts/run_model.R -c $TAXDIR/configs/202110gavi-3/no-vaccination/COD_no-vaccination_30.yml
$RSCRIPT $TAXDIR/scripts/run_model.R -c $TAXDIR/configs/202110gavi-3/no-vaccination/ETH_no-vaccination_30.yml
$RSCRIPT $TAXDIR/scripts/run_model.R -c $TAXDIR/configs/202110gavi-3/no-vaccination/HTI_no-vaccination_30.yml
$RSCRIPT $TAXDIR/scripts/run_model.R -c $TAXDIR/configs/202110gavi-3/no-vaccination/YEM_no-vaccination_30.yml

echo "End of script"
date

