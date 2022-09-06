#!/bin/bash
#SBATCH --job-name=gen_incid
#SBATCH --output=sbatch_logs/pipeline_run_%A_%a_%u.log
#SBATCH --time=10000:00:00
#SBATCH --mem=10G
#SBATCH --array=0-43%5
#SBATCH --nodelist=idmodeling2

echo "Beginning of script"
date

TAXDIR=/home/.../gavi_vimc_cholera
CONFIGDIR=configs/202110gavi-3/campaign-default/district-estimate/0.001
RSCRIPT=/opt/R/4.0.3/bin/Rscript

cd $TAXDIR
CONFIGNAMES=($(ls $TAXDIR/$CONFIGDIR | tr ' ' '\n'))
echo $RSCRIPT $TAXDIR/scripts/run_country_incid_crop.R -c $TAXDIR/$CONFIGDIR/${CONFIGNAMES[$SLURM_ARRAY_TASK_ID]}
$RSCRIPT $TAXDIR/scripts/run_country_incid_crop.R -c $TAXDIR/$CONFIGDIR/${CONFIGNAMES[$SLURM_ARRAY_TASK_ID]} || exit 1

echo "End of script"
date
