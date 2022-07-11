#!/bin/bash
#SBATCH --job-name=diag
#SBATCH --output=log/Diagnostic_report_output_log.txt
#SBATCH --time=10000:00:00
#SBATCH --mem=20G
#SBATCH --nodelist=idmodeling2
#SBATCH -c 4

echo "Beginning of script"
date
TAXDIR=/home/kaiyuezou/VIMC_Model/surveillance_project/gavi_vimc_cholera
RSCRIPT=/opt/R/4.0.3/bin/Rscript

cd $TAXDIR
$RSCRIPT -e "rmarkdown::render('diagnostics/surveillance_project/diagnostic_report.Rmd', params=list(args = 'myarg'))" || exit 1

echo "End of script"
date
