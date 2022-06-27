#!/bin/bash
#SBATCH --job-name=run_model
#SBATCH --time=10000:00:00
#SBATCH --ntasks=1
#SBATCH --mem=10G
#SBATCH --array=0-3%4
#SBATCH --nodelist=idmodeling2

echo "Beginning of script"
date

TAXDIR=/home/kaiyuezou/VIMC_Model/202110gavi-3/gavi_vimc_cholera
CONFIGDIR=configs/202110gavi-3/campaign-default/zzz_setting2_new_round2
RSCRIPT=/opt/R/4.0.3/bin/Rscript

cd $TAXDIR
CONFIGNAMES=($(ls $TAXDIR/$CONFIGDIR | tr ' ' '\n'))
echo $RSCRIPT $TAXDIR/scripts/run_model.R -c $TAXDIR/$CONFIGDIR/${CONFIGNAMES[$SLURM_ARRAY_TASK_ID]}
$RSCRIPT $TAXDIR/scripts/run_model.R -c $TAXDIR/$CONFIGDIR/${CONFIGNAMES[$SLURM_ARRAY_TASK_ID]} || exit 1

echo "End of script"
date
