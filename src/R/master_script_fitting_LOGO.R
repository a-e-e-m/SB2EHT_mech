## script to prepare LOGO for a specified row in the file models2run.csv
# take 'run' from models2run.csv and create models2run_LOGO_run.csv with one row per LOGO iteration

# NOTE: This only works if there is an existing folder with data for the run (wrt to models2run.csv)!

## NOTE: If you work on the .csv file with libreOffice calc, choose always delimeter ';', 
## string delimeter '"' and when opening 'Format quoted fields as text'
## use this to copy folder to cluster: rsync -vzre ssh ./fitting/dose_response_1_7_rerun_CvforEHTmortality_adaptdelta_moreESS denadr00@transfer12.scicore.unibas.ch:/scicore/home/chitnis/denadr00/lAIRepmi/BA2EHT/fitting/
## and this one to copy from cluster to local: rsync -vzre ssh denadr00@transfer12.scicore.unibas.ch:/scicore/home/chitnis/denadr00/lAIRepmi/BA2EHT/fitting/dose_response_1_7_rerun_CvforEHTmortality_adaptdelta* ./fitting/

library("dplyr")
library("lubridate")


# define which models to run
###############################
# read file with all model runs
models2run <- read.csv(file = './models2run.csv', sep = ';', stringsAsFactors = FALSE)

# # # get arguments passed to this script
# args = commandArgs(trailingOnly=TRUE)
# j <- as.integer(args[1])
# where <- as.character(args[2])
# if (where=="cluster"){cluster <- TRUE} else if (where=="local"){cluster<- FALSE}

# OR set manually
# select one line from models2run.csv
j <- 4
# cluster <- FALSE
cluster <- TRUE

# # if needed: Install required packages 
# source("./src/R/packages_cluster.R")


# function for runtime in sec to qos
qos <- function(runtime){
  case_when(runtime <= 30*60 ~ "30min", 
            runtime <= 6*60*60 ~ "6hours", 
            runtime <= 6*60*60 ~ "1day"
  )
}

#function to capture slurm job ID
submit <- function(cmd) {
  cat("\nSubmitting:\n", cmd, "\n")
  jobid <- system(cmd, intern = TRUE)
  jobid <- trimws(jobid[1])
  cat("Submitted job:", jobid, "\n")
  jobid
}

# log folder with dynamic name
logdir <- file.path(
  "fitting",
  models2run$run[j],
  "slurm_logs_LOGO"
)
dir.create(logdir, recursive = TRUE, showWarnings = FALSE)




# get number of groups (nT)
# nT equals the number of model runs for LOGO-CV
trials <- readRDS(file.path(".", "fitting", models2run$run[j],  "trials.rds"))
nT <- trials$nT 

# load or create models2run_LOGO_run 
file_models2run_LOGO_run <- file.path(".", "models2run4LOGO", paste0("models2run_logo_", models2run$run[j], ".csv"))

if (file.exists(file_models2run_LOGO_run)){
  models2run_LOGO_run <- read.csv2(file = file_models2run_LOGO_run, stringsAsFactors = FALSE)
} else {
  models2run_LOGO_run <- tidyr::expand_grid(models2run[j,], LOGO_ID = seq(1,nT) ) |>
  rename(run_primary = run) |>
  mutate(run = paste0(run_primary, "_LOGO_", LOGO_ID),
         take_data_from_run = run_primary) |>
  relocate(LOGO_ID, run, .after = run_primary)

write.csv2(models2run_LOGO_run, file = file_models2run_LOGO_run, row.names=FALSE)
}


# create a folder per LOGO_ID
###############################
for (k in seq(1,nT)){
  # stops and throwns error if folder already exists.
  if (dir.exists(file.path('fitting', models2run_LOGO_run$run[k]))){
    stop("Folder for run with the same name already exists!")
  } else {
    dir.create(file.path('fitting', models2run_LOGO_run$run[k]), showWarnings = TRUE)
  }
}




# process data and stan control input 
###############################
# CAREFUL: Make sure to call data_processing-job_LOGO.R here, calling data_processing-job.R also works but would put all data in the folder, so no LOGO-CV!!! 
# Check this by looking at the intervention sample size numbers in data_summary.csv!!!
# 1st argument (k) gives index of run within LOGO-CV, 2nd argument (file_models2run_LOGO_run) gives models2run_LOGO_run (with k looping through the rows in there)
if (cluster){
  data_job <- submit(paste(
    "sbatch --parsable",
    paste0("--cpus-per-task=", "1"),
    paste0("--qos=", "30min"),
    paste0("--time=", "00:15:00"),
    "--job-name=LOGO-data",
    paste0("--output=", shQuote(file.path(logdir, "%x_%j.out"))),
    paste0("--error=",  shQuote(file.path(logdir, "%x_%j.err"))),
    "src/bash/LOGO-data-processing.sh",
    nT,
    shQuote(file_models2run_LOGO_run)
  ))
} else {
  for (k in seq(1,nT)){
    cmd <- paste0("Rscript ./src/R/data_processing-job_LOGO.R ", k, " ", file_models2run_LOGO_run)
    system(cmd)
  }
}



# # data visualisation and summary
# ###############################
# cmd <- paste0("Rscript ./src/R/data_visualisation_rich.R ", k)
# system(cmd)


# point estimate / optimiser
###############################
if (cluster){
  opt_job <- submit(paste(
    "sbatch --parsable",
    paste0("--dependency=afterok:", data_job),
    paste0("--cpus-per-task=", "1"),
    paste0("--qos=", "30min"),
    paste0("--time=", "00:30:00"),
    "--job-name=LOGO-opt",
    paste0("--output=", shQuote(file.path(logdir, "%x_%j.out"))),
    paste0("--error=",  shQuote(file.path(logdir, "%x_%j.err"))),
    "src/bash/LOGO-estimating.sh",
    nT,
    shQuote(file_models2run_LOGO_run)
  ))
} else {
  for (k in seq(1,nT)){
    cmd <- paste0("Rscript ./src/R/estimating-job.R ", k, " ", file_models2run_LOGO_run)
    system(cmd)
  }
}


# # check (actual vs predicted) point estimate / optimiser
# ###############################
# for (k in seq(1,nT)){
#   cmd <- paste0("Rscript ./src/R/estchecking-job.R ", k, " ", file_models2run_LOGO_run)
#   system(cmd)
# }


# copying files to cluster
###############################
# files to sync
if (!cluster){
  files_to_sync <- c(
  file.path("fitting", models2run_LOGO_run$run),
  file_models2run_LOGO_run,
  file.path("fitting", models2run$run[j])
  )
  
  # ensure directories have trailing "/"
  is_dir <- dir.exists(files_to_sync)
  files_to_sync[is_dir] <- paste0(files_to_sync[is_dir], "/")
  
  # write txt file listing files to sync 
  writeLines(unique(files_to_sync), "files_to_sync.txt")
  
  # and run the bash script ./src/bash/copy__files_to_sync__to_cluster.sh with argument "./files_to_sync.txt" or run here
  system2("bash", c("./src/bash/copy__files_to_sync__to_cluster.sh", "files_to_sync.txt"))
}



# running stan
###############################
if (cluster){
  N_cpu <- max(models2run_LOGO_run$chains, na.rm = TRUE)   # safe upper bound
  runtime <- "00:30:00"
  qos <- case_when(hms(runtime) <= minutes(30) ~ "30min", 
                   hms(runtime) <= hours(6) ~ "6hours", 
                   hms(runtime) <= days(1) ~ "1day",
                   hms(runtime) <= days(7) ~ "1week"
  )
  
  sample_job <- submit(paste(
    "sbatch --parsable",
    paste0("--dependency=afterok:", opt_job),
    paste0("--array=1-", nT),
    paste0("--cpus-per-task=", N_cpu),
    paste0("--qos=", qos),
    paste0("--time=", runtime),
    "--job-name=LOGO-stan",
    paste0("--output=", shQuote(file.path(logdir, "%x_%A_%a.out"))),
    paste0("--error=",  shQuote(file.path(logdir, "%x_%A_%a.err"))),
    "src/bash/fitting-cluster.sh",
    shQuote(file_models2run_LOGO_run)
  ))
} else {
  sbatch_cmd <- paste0(
    "sbatch",
    " --array=1-", nT,
    " --cpus-per-task=", N_cpu,
    " --qos=", qos,
    " --time=", runtime,
    " --job-name=run-stan",
    " --output=./fitting/%x_%A_%a.out",
    " src/bash/fitting-cluster.sh ",
    file_models2run_LOGO_run
  )
  cat(sbatch_cmd, "\n")
}


# copying files to local
###############################
if (!cluster){
  for (k in seq(1,nT)){
    cmd <- paste0("rsync -vzre ssh denadr00@transfer12.scicore.unibas.ch:/scicore/home/chitnis/denadr00/lAIRepmi/BA2EHT/fitting/", models2run_LOGO_run$run[k], "/", " ./fitting/", models2run_LOGO_run$run[k], "/")
    system(cmd)
  }
  
  # or use this in the terminal to sync all at once
  cmd <- paste0("rsync -vzre ssh denadr00@transfer12.scicore.unibas.ch:/scicore/home/chitnis/denadr00/lAIRepmi/BA2EHT/fitting/", models2run$run[j], "_LOGO_*", " ./fitting/")
  cat(cmd, "\n")
}



# stan sampling diagnostics
###############################
if (cluster){
  diag_job <- submit(paste(
    "sbatch --parsable",
    paste0("--dependency=afterok:", sample_job),
    paste0("--cpus-per-task=", "1"),
    paste0("--qos=", "30min"),
    paste0("--time=", "00:10:00"),
    "--job-name=logo-diag",
    paste0("--output=", shQuote(file.path(logdir, "%x_%j.out"))),
    paste0("--error=",  shQuote(file.path(logdir, "%x_%j.err"))),
    "src/bash/LOGO-diagnostics.sh",
    nT,
    shQuote(file_models2run_LOGO_run)
  ))
}else{
  for (k in seq(1,nT)){
    cmd <- paste0("Rscript ./src/R/samplingdiagnostics-job.R ", k, " ", file_models2run_LOGO_run)
    system(cmd)
  }
}



# computing LOGO predictions
###############################
if (cluster){
  pred_job <- submit(paste(
    "sbatch --parsable",
    paste0("--dependency=afterok:", diag_job),
    paste0("--cpus-per-task=", "1"),
    paste0("--mem=", "16G"),
    paste0("--qos=", "30min"),
    paste0("--time=", "00:20:00"),
    "--job-name=logo-pred",
    paste0("--output=", shQuote(file.path(logdir, "%x_%j.out"))),
    paste0("--error=",  shQuote(file.path(logdir, "%x_%j.err"))),
    "src/bash/LOGO-predictions.sh",
    j
  ))
}else{
    cmd <- paste0("Rscript ./src/R/LOGO_CV_prediction.R ", j)
    system(cmd)
}







