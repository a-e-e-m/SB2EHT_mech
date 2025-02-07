## script to run models according to rows in the file models2run.csv

## NOTE: If you work on the .csv file with libreOffice calc, choose always delimeter ';', 
## string delimeter '"' and when opening 'Format quoted fields as text'

library("lubridate")
library("dplyr")

# define which models to run
###############################
# read file with all model runs
models2run <- read.csv(file = './models2run.csv', sep = ';', stringsAsFactors = FALSE)
  
# get arguments passed to this script if called from the console with Rscript
args = commandArgs(trailingOnly=TRUE)
j <- as.integer(args[1])

# # OR select one line from models2run.csv
# j <- 8

# # if needed: Install required packages 
# system('Rscript inst_pack.R')


# create folder structure
###############################
# stops and throwns error if folder already exists.
if (dir.exists(file.path('fitting', models2run$run[j]))){
  stop("Folder for run with the same name already exists!")
} else {
  dir.create(file.path('fitting', models2run$run[j]), showWarnings = TRUE)
}

## all commands to copy into console
 cmd <- paste0("Rscript ./src/R/data_processing-job.R ", j, " && Rscript ./src/R/estimating-job.R ", j, " && Rscript ./src/R/fitting-job.R ", j,  " && Rscript ./src/R/samplingdiagnostics-job.R ", j, " && Rscript ./src/R/mech_model_specific/fitchecking-job_mech_model.R ", j)

# run from Rstudio
system(cmd)


