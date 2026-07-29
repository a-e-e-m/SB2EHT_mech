# script to prepare data for leave-one-group-out cross validation ('LOGO-CV')
# ===========================
# Parameter passed to this script: integer indicating row of models2run_logo_run and file name (including relative path) of models2run_logo_run
# ===========================
# This script processes the data for a specific model run as defined in models2run.
# It copies the data from another run.
# Otherwise it calls the functions BA_data_processing and EHT_data_processing which read in, check and correct/allign the data if necessary. 
# Subsequently, the data is selected according to the run specification, and variables are aligned across BA and EHT.
# Finally, the BA and EHT data is paired by adding group_number as id in both BA and EHT data.
#
# Notes: 
# on data input side I use the term data_set, once the data is processed I use the term trial.

## config
################################
library(tidyverse)
library(reshape2)
library(lubridate)

# source script containing needed functions
functionsfolder <- file.path('./src/R/functions')
files.source <- list.files(functionsfolder)
invisible(sapply(files.source, function(x) source(paste0(functionsfolder, "/", x))))


## get instructions
################################
# get model index from argument passed to this script
args = commandArgs(trailingOnly=TRUE)
i <- as.numeric(args[1])
if (length(args) < 2){
  stop("No csv file in models2run4LOGO was specified in call to ./src/R/data_processing-job_LOGO.R.")
} else {
  file_models2run <- args[2]
}


# read in job specification from models2run
models2run <- read.csv2(file=file_models2run, stringsAsFactors=FALSE)

# define variables from this
run <- as.character(models2run$run[i])
take_data_from_run <- as.character(models2run$take_data_from_run[i])
if (str_detect(models2run$data_sets[i], "all")){
  list_data_sets <- readxl::read_xlsx(file.path("data", "list_data_sets.xlsx"), trim_ws = FALSE)
  data_sets <- sort(unique(as.integer(list_data_sets$Trial_code)))
} else {
  list_data_sets <- readxl::read_xlsx(file.path("data", "list_data_sets.xlsx"), trim_ws = FALSE)
  data_sets <- eval(parse(text=paste0('c(', models2run$data_sets[i], ")")))
}
bioassay_type <- unlist(strsplit(models2run$bioassay_type[i], split=", "))
with_feeding <- as.logical(models2run$with_feeding[i])
BA_only <- as.logical(models2run$BA_only[i])
single_pyrethroid_input <- as.logical(models2run$BA_single_pyrethroid[i])
net_gen_input <- unlist(strsplit(models2run$net_gen[i], split=", "))
net_status_input <-  unlist(strsplit(models2run$net_status[i], split=", "))
file_outcome_timepoints <- as.character(models2run$outcome_timepoint[i])
group_by_vars <- unlist(strsplit(models2run$group_by[i], split=", "))
sumup_species <- as.logical(models2run$sumup_species[i])
manual_match_functions <- unlist(strsplit(models2run$manual_match_functions[i], split=", "))
include_controls <- as.logical(models2run$include_controls[i])
CV_samples_sizes <- eval(parse(text=paste0('c(', models2run$CV_samples_sizes[i], ")")))
LOGO_ID <- as.numeric(models2run$LOGO_ID[i])

# print model specification to output and to file in specific folder
models2run[i,]
write.csv(models2run[i,], file = paste0("./fitting/", run, "/", "job.csv"))



# take data from specified run and exclude EHT intervention arm data for LOGO_ID assay pair
if (nchar(take_data_from_run) > 0){
  file.copy(file.path("fitting", take_data_from_run, "trials.rds"), file.path("fitting", run, "trials.rds"), overwrite = TRUE)
  B_all <- readRDS(file.path("fitting", take_data_from_run,  "B_all.rds"))
  saveRDS(B_all, file.path("fitting", run, "B_all.rds"))
  
  # for runs including EHT mortality data
  if (!BA_only){
    H_all <- readRDS(file.path("fitting", take_data_from_run,  "H_all.rds"))
    
    H_all <- H_all |>
      filter(!(!control & group_number == LOGO_ID))
    
    saveRDS(H_all, file.path("fitting", run, "H_all.rds"))
  }
  
  # for runs including EHT feeding data
  if (with_feeding){
    H_f_all <- readRDS(file.path("fitting", take_data_from_run,  "H_f_all.rds"))
    
    H_f_all <- H_f_all |>
      filter(!(!control & group_number == LOGO_ID))
    
    saveRDS(H_f_all, file.path("fitting", run, "H_f_all.rds"))
  }
  
  
  # output table summarising both BA and EHT data
    B_all_sum <- B_all |>
      group_by(across(all_of(c(setdiff(group_by_vars, "insecticide"), "group_number", "Trial_code", "int_dose_available", "test_type")))) |> 
      summarise(
        "Insecticides tested in resistance bio assay" = toString(unique(insecticide[!control])),
        "Number of replicates for intervention resistance bio assays" = sum(!control),
        "Total number of mosquitoes in intervention resistance bio assays" = sum(N_b[!control]),
        "Number of replicates for control resistance bio assays" = sum(control),
        "Total number of mosquitoes in control resistance bio assays" = sum(N_b[control]),
        .groups = "keep"
      ) |>
      rename("Resistance bio assay site" = site) |>
      arrange(group_number)
    
    H_all_sum <- H_all |>
      group_by(across(all_of(c(setdiff(group_by_vars, "insecticide"), "group_number", "Trial_code", "hut_type", "outcome_timepoint")))) |> 
      summarise(
        "Net products" = toString(unique(net_product[!control])),
        "Insecticides in nets" = toString(unique(insecticide[!control])),
        "Number of replicates for intervention arm of EHT" = sum(!control),
        "Total number of mosquitoes in intervention arm of EHT" = sum(N_h[!control]),
        "Number of replicates for control arm of EHT" = sum(control),
        "Total number of mosquitoes in control arm of EHT" = sum(N_h[control]),
        .groups = "keep"
      ) |>
      rename("EHT site" = site) |>
      arrange(group_number, country, "bio assay site", "EHT site" )     
    
    publication_lookup <- list_data_sets |>
      filter(!is.na(Trial_code)) |>
      select(c("Trial_code", "Senior Author", "PubYear", "Title", "Publication/source", "Link")) |>
      mutate(
        Trial_code = as.integer(Trial_code)
      )
    
    summary_table <- B_all_sum |>
      full_join(H_all_sum) |>
      left_join(
        publication_lookup,
        by = c("Trial_code")
      ) |>
      ungroup() |>
      # select(!Trial_code) |>
      mutate(
        Reference = NA, 
        int_dose_available = ifelse(int_dose_available == 1, "yes", "no")
      ) |>
      rename(
        "Intensity dose resistance bio assay" = int_dose_available,
        "Group number" = group_number,
        Country = country, 
        Year = year
      ) |>
      relocate("Group number", "Trial_code", Reference,
               "Senior Author", "PubYear", "Title", "Publication/source", "Link",
               Country, "Resistance bio assay site", "EHT site", Year,  
               "test_type",
               "Intensity dose resistance bio assay",
               "Insecticides tested in resistance bio assay", 
               "hut_type", "outcome_timepoint",
               "Insecticides in nets", "Net products",
               "Number of replicates for intervention resistance bio assays", "Total number of mosquitoes in intervention resistance bio assays",
               "Number of replicates for control resistance bio assays", "Total number of mosquitoes in control resistance bio assays",
               
               "Number of replicates for intervention arm of EHT","Total number of mosquitoes in intervention arm of EHT", 
               "Number of replicates for control arm of EHT", "Total number of mosquitoes in control arm of EHT"
      )
    
    write.csv(summary_table, file = file.path("fitting", run, "data_summary.csv"), row.names=FALSE)
  
} else {
  stop("Error: No run to take data from is specified!")
}