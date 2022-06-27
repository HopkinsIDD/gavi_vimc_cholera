#!/bin/bash
#SBATCH --job-name=novax
#SBATCH --time=10000:00:00
#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --array=0-46%5
#SBATCH --nodelist=idmodeling2

echo "Beginning of script"
date

TAXDIR=/home/kaiyuezou/VIMC_Model/202110gavi-3/gavi_vimc_cholera
CONFIGDIR=configs/202110gavi-3/no-vaccination
RSCRIPT=/opt/R/4.0.3/bin/Rscript

cd $TAXDIR
CONFIGNAMES=($(ls $TAXDIR/$CONFIGDIR | tr ' ' '\n'))
echo $RSCRIPT $TAXDIR/scripts/run_model.R -c $TAXDIR/$CONFIGDIR/${CONFIGNAMES[$SLURM_ARRAY_TASK_ID]}
$RSCRIPT $TAXDIR/scripts/run_model.R -c $TAXDIR/$CONFIGDIR/${CONFIGNAMES[$SLURM_ARRAY_TASK_ID]} || exit 1

echo "End of script"
date
