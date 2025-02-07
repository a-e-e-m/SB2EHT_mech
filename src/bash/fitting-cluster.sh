#!/bin/bash
#SBATCH --job-name=run-stan
#SBATCH --mem-per-cpu=8G

###########################################
# Script to execute rscript-job.R 
# with input taken from row ${1} of file 'models2run.csv'
###########################################
# SBATCH --qos=1day
# SBATCH --time=12:00:00
# SBATCH --qos=6hours
# SBATCH --time=06:00:00

module purge 
module load R/4.1.0-foss-2018b

echo "RUNNING TASK ${1}"
echo "Command -----"
cmd="Rscript ./src/R/fitting-job.R ${1} && Rscript ./src/R/samplingdiagnostics-job.R $1"
echo $cmd
echo "Log -----"
eval $cmd
