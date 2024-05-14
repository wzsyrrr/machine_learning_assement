#!/bin/bash
#SBATCH --account=SSCM030364
#SBATCH --job-name=model
#SBATCH --partition=veryshort
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=6:00:00 #SBATCH --mem=19000

module add languages/r/4.3.1
echo "Beginning of my job steps" date
Rscript code/model.r date
echo "End of my job steps"