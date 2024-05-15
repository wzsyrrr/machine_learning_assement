#!/bin/bash
#SBATCH --account=SSCM030364
#SBATCH --job-name=analysis
#SBATCH --partition=veryshort
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=6:00:00 #SBATCH --mem=15000

module add languages/r/4.3.1
echo "Beginning of my job steps" date
Rscript code/analysis.r date
echo "End of my job steps"