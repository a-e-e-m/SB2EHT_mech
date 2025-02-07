#!/bin/bash
# A bash script to run a R script 
Rscript ./src/R/fitting-job.R ${1} && Rscript ./src/R/samplingdiagnostics-job.R $1
