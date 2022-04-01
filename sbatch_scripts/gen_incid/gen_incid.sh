#!/bin/bash
#SBATCH --job-name=gen_incid
#SBATCH --time=10000:00:00
#SBATCH --mem=20G
#SBATCH --array=0-46%1
#SBATCH --nodelist=idmodeling2
#SBATCH -c 4

echo "Beginning of script"
date

TAXDIR=/home/kaiyuezou/VIMC_Model/202110gavi-3/gavi_vimc_cholera
CONFIGDIR=configs/202110gavi-3/campaign-default
RSCRIPT=/opt/R/4.0.3/bin/Rscript

cd $TAXDIR
CONFIGNAMES=($(ls $TAXDIR/$CONFIGDIR | tr ' ' '\n'))
echo $RSCRIPT $TAXDIR/scripts/run_country_incid_crop.R -c $TAXDIR/$CONFIGDIR/${CONFIGNAMES[$SLURM_ARRAY_TASK_ID]}
$RSCRIPT $TAXDIR/scripts/run_country_incid_crop.R -c $TAXDIR/$CONFIGDIR/${CONFIGNAMES[$SLURM_ARRAY_TASK_ID]} || exit 1

echo "End of script"
date
