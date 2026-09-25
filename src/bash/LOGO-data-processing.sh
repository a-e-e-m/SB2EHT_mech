#!/bin/bash
#SBATCH --job-name=LOGO-data-processing

set -euo pipefail

module purge
module load R/4.4.2-foss-2024a

N_T="$1"
MODELS_LOGO_FILE="$2"

for k in $(seq 1 "$N_T"); do
  echo "Processing data for LOGO k=${k}"
  Rscript ./src/R/data_processing-job_LOGO.R "$k" "$MODELS_LOGO_FILE"
done