#!/bin/bash
#SBATCH --account=SSCM030364
#SBATCH --job-name=analysis
#SBATCH --partition=veryshort
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1:00:00 #SBATCH --mem=19000

module load python
pip3 install gdown
echo "start downloading clinical"
gdown 1ivYpFh5a2cRbQ42qbZq3vmM45-_xjmNk
echo "downloaded clinical.txt"

echo "start downloading methylation"
gdown 15qGOW8_hLYeGPpzuJhMVBBy-VcPxke6Y
echo "downloaded methylation"

echo "start downloading mirna"
gdown 1K1y92ZOpgt72Q2LhYYDvjoLGI8FZSPp1
echo "downloaded mirna"

echo "start downloading mrna"
gdown 1ZICrKLKZjCavu7eXFedq6bCAmQSetIa9
echo "downloaded mrna"

echo "start downloading mutations"
gdown 1CdlxjsDfx3rC6BDNUZuVG9fx-iEU2RV8
echo "downloaded mutations"

echo "start downloading protein"
gdown 148pVzaCiWZLyilDUI09-p12lOu64mDDV
echo "downloaded protein"

mv clinical.txt data_raw
mv methylation.txt data_raw
mv mirna.txt data_raw
mv mrna.txt data_raw
mv mutations.txt data_raw
mv protein.txt data_raw