#!/bin/bash
#SBATCH --job-name=LOGO-predictions

set -euo pipefail

module purge
module load R/4.4.2-foss-2024a

i="$1"

  echo "computing LOGO predictions for run ${i}"
  Rscript ./src/R/LOGO_CV_prediction.R "$i"