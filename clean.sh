#!/bin/bash
#SBATCH --account=SSCM030364
#SBATCH --job-name=clean
#SBATCH --partition=veryshort
#SBATCH --nodes=20
#SBATCH --ntasks-per-node=1
#SBATCH --time=1:00:00 #SBATCH --mem=19000

module add languages/r/4.3.1
echo "Beginning of my job steps" date
Rscript code/data_clean.r date
echo "End of my job steps"