# ===========================
# rscript-job.R
# 
# prepares data and runs stan according to models2run.csv
#
# Parameter passed to this script: integer indicating row of models2run.csv
# ===========================


## config
################################
library(tidyverse)
library(reshape2)
library(lubridate)

# source script containing needed functions
functionsfolder <- file.path('./src/R/functions')
files.source <- list.files(functionsfolder)
invisible(sapply(files.source, function(x) source(paste0(functionsfolder, "/", x))))

# get model index from argument passed to this script
args = commandArgs(trailingOnly=TRUE)
i <- as.numeric(args[1])

# read in job specification from models2run
models2run <- read.csv2(file='./models2run.csv', stringsAsFactors=FALSE)
run <- as.character(models2run$run[i])
LOO <- as.numeric(models2run$LOO[i])
exactLOO <- as.numeric(models2run$exactLOO[i])
#NODATA <- as.numeric(models2run$NODATA[i])

# print model specification to output and to file in specific folder
models2run[i,]
write.csv(models2run[i,], file = paste0("./fitting/", run, "/", "job.csv"))

sink(file = paste0("./fitting/", run, "/", "stanfit_log.txt"), split = TRUE)

## fit model including validation
#################################

result <- validfit(model = as.character(models2run$model[i]),
                   init_function = as.character(models2run$init_function[i]),
                   with_feeding <- as.logical(models2run$with_feeding[i]),
                   BA_only <- as.logical(models2run$BA_only[i]),
                   add_data = as.character(models2run$add_data[i]),
                   LOO = as.numeric(models2run$LOO[i]),
                   exactLOO = as.numeric(models2run$exactLOO[i]),
                   pareto_k_threshold = as.numeric(models2run$pareto_k_threshold[i]),
                   LOO_CV = as.numeric(models2run$LOO_CV[i]),
                   LOO_CV_groups = as.character(models2run$LOO_CV_groups[i]),
                    chains = as.numeric(models2run$chains[i]),
                    warmup = as.numeric(models2run$warmup[i]),
                    iter = as.numeric(models2run$iter[i]),
                    run = as.character(models2run$run[i]),
                    control_list_name = as.character(models2run$control_list_name[i])
                    )

sink()

								
