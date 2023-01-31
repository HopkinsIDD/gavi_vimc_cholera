#!/bin/bash
#SBATCH --job-name=gen_incid
#SBATCH --time=2-23:59:00
#SBATCH --mem=10G
#SBATCH --array=0-43%10
#SBATCH --partition=defq
#SBATCH --account=aazman1


echo "Beginning of script"
date
pwd
echo "End of script"
date
