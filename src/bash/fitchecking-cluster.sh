#!/bin/bash
#SBATCH --job-name=run-stan
#SBATCH --mem-per-cpu=16G
#SBATCH --qos=6hours
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1

###########################################
# Script to execute rscript-job.R 
# with input taken from row ${1} of file 'models2run.csv'
###########################################

module purge 
module load R/4.1.0-foss-2018b

echo "RUNNING TASK ${1}"
echo "Command -----"
cmd="Rscript ./src/R/fitchecking-job.R ${1}"
echo $cmd
echo "Log -----"
eval $cmd
