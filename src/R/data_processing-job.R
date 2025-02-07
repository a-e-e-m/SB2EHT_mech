# ===========================
# Parameter passed to this script: integer indicating row of models2run.csv
# ===========================
# This script processes the data for a specific model run as defined in models2run.
# If specified, it copies the data from another run.
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

# read in job specification from models2run
models2run <- read.csv2(file='./models2run.csv', stringsAsFactors=FALSE)

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

# print model specification to output and to file in specific folder
models2run[i,]
write.csv(models2run[i,], file = paste0("./fitting/", run, "/", "job.csv"))


# if existing run to take data from is specified, copy data and skip rest of script
if (nchar(take_data_from_run) > 0){
  file.copy(file.path("fitting", take_data_from_run, "B.rds"), file.path("fitting", run, "B.rds"), overwrite = TRUE)
  file.copy(file.path("fitting", take_data_from_run, "B_all.rds"), file.path("fitting", run, "B_all.rds"), overwrite = TRUE)
  file.copy(file.path("fitting", take_data_from_run, "data_summary.csv"), file.path("fitting", run, "data_summary.csv"), overwrite = TRUE)
  file.copy(file.path("fitting", take_data_from_run, "trials.rds"), file.path("fitting", run, "trials.rds"), overwrite = TRUE)
  
  if (!BA_only){
    file.copy(file.path("fitting", take_data_from_run, "H.rds"), file.path("fitting", run, "H.rds"), overwrite = TRUE)
    file.copy(file.path("fitting", take_data_from_run, "H_all.rds"), file.path("fitting", run, "H_all.rds"), overwrite = TRUE)
  }
  
  
  if (with_feeding){
    file.copy(file.path("fitting", take_data_from_run, "H_f.rds"), file.path("fitting", run, "H_f.rds"), overwrite = TRUE)
    file.copy(file.path("fitting", take_data_from_run, "H_f_all.rds"), file.path("fitting", run, "H_f_all.rds"), overwrite = TRUE)
  }
  
  file.copy(file.path("fitting", take_data_from_run,  paste0("group_CVsample_5",".rds")), file.path("fitting", run,  paste0("group_CVsample_5",".rds")), overwrite = TRUE )
  file.copy(file.path("fitting", take_data_from_run,  paste0("group_CVsample_8",".rds")), file.path("fitting", run,  paste0("group_CVsample_8",".rds")), overwrite = TRUE )

} else {
  ## process data
  #############################
  # write output to a txt file
  sink(file.path("./fitting", run, "data_processing.txt"))
  
    # process BA data
    list2env(BA_data_processing(data_sets = data_sets, file_list_data_sets = "data/list_data_sets.xlsx"), envir = .GlobalEnv)
    BA_unfiltered <- BA
    
    # process EHT data
    list2env(EHT_data_processing(data_sets = data_sets, file_list_data_sets = "data/list_data_sets.xlsx"), envir = .GlobalEnv)
    EHT_unfiltered <- EHT
    
    ## filter EHT data according to instructions (run specification)
      EHT <- EHT_unfiltered %>%
        filter(net_gen %in% net_gen_input | control,
               net_status %in% net_status_input | control) 
      
      if (!include_controls){
        EHT <- EHT %>%
          filter(!control)   
      }
      
      # select with respect to what time point outcomes are taken
      outcome_timepoint_selection <- read.csv2(file = file.path("data", file_outcome_timepoints))
      EHT <- EHT %>%
        semi_join(select(outcome_timepoint_selection, -notes))
      
      #print excluded EHT data to output
      EHT_excluded_by_specs <- dplyr::setdiff(EHT_unfiltered, EHT) %>% group_by(across( c("Trial_code", "treat_name", "control", 
                                                                                          "net_product", "net_gen", "insecticide", "net_status", "N_washes", "outcome_timepoint", 
                                                                                          "species", "country", "site", "hut_type", "holes", 
                                                                                          "year"))) %>% 
                                                                       summarise(n_datapoints = n())
  
      EHT_excluded_by_specs_short <- EHT_excluded_by_specs %>% group_by(across( c("Trial_code", "control", 
                                                                                  "net_product", "net_gen", "insecticide", 
                                                                                  "country", "site", "hut_type",
                                                                                  "year"))) %>% 
                                                               summarise(n_net_status = n_distinct(net_status), 
                                                                         n_N_washes = n_distinct(N_washes), 
                                                                         n_outcome_timepoint = n_distinct(outcome_timepoint), 
                                                                         n_species = n_distinct(species), 
                                                                         n_holes = n_distinct(holes), 
                                                                         n_datapoints = sum(n_datapoints))
      
      EHT_excluded_by_specs_veryshort <- EHT_excluded_by_specs %>% group_by(Trial_code) %>% 
                                                                  summarise(n_control = n_distinct(control),
                                                                            n_net_product = n_distinct(net_product),
                                                                            n_net_gen = n_distinct(net_gen),
                                                                            n_insecticide = n_distinct(insecticide),
                                                                            n_country = n_distinct(country),
                                                                            n_site = n_distinct(site),
                                                                            n_hut_type = n_distinct(hut_type),
                                                                            n_year = n_distinct(year),
                                                                            n_net_status = n_distinct(net_status), 
                                                                            n_N_washes = n_distinct(N_washes), 
                                                                            n_outcome_timepoint = n_distinct(outcome_timepoint), 
                                                                            n_species = n_distinct(species), 
                                                                            n_holes = n_distinct(holes), 
                                                                            n_datapoints = sum(n_datapoints))
      
      
      print("The following EHT data was excluded by run specification (summarised with n_somevariable meaning 'number of distinct values for variable 'somevariable').")
      print.data.frame(EHT_excluded_by_specs_veryshort)
      
      # sumup species if sumup_species set to TRUE (and if no grouping by species)
      # WARNING: If there are separate data points that cannot be distinguished based on covars,
      # then this operation pools data points unintentionally!
      if ( !('species' %in% group_by_vars) & sumup_species) {
        EHT_test <- EHT %>%
          group_by(across(setdiff(covars_EHT, "species"))) %>%
          summarise(across(all_of(outcomes_EHT), ~ sum(.x))) |>
          ungroup()
      }    
      
      
    ## filter BA data according to instructions (run specification)
      BA <- BA_unfiltered 
      
      if (!include_controls){
        BA <- BA %>%
          filter(!control)   
      }
      
      if (single_pyrethroid_input){
        BA <- BA %>%
          filter(single_pyrethroid | control)   
      }
      
      if ("disc" %in% bioassay_type & !("int" %in% bioassay_type)){
        BA <- BA %>%
          filter(is_disc_dose | control) 
      } 
    
      #print excluded BA data to output
      BA_excluded_by_specs <- dplyr::setdiff(BA_unfiltered, BA) %>% group_by(across( c("Trial_code", "control", 
                                                                                       "test_type", "insecticide", "single_pyrethroid", "BA_total_known", 
                                                                                       "strain", "country", "site", "year",
                                                                                       "species",
                                                                                       "conc_perc", "conc_mug", "is_disc_dose", "times_disc_dose"))) %>% 
        summarise(n_datapoints = n())
      
      BA_excluded_by_specs_short <- BA_excluded_by_specs %>% group_by(across( c("Trial_code", "control", 
                                                                                "test_type", "insecticide", "single_pyrethroid", "BA_total_known", 
                                                                                "strain", "country", "site", "year",
                                                                                "conc_perc", "conc_mug", "is_disc_dose", "times_disc_dose"))) %>% 
        summarise(n_species = n_distinct(species), 
                  n_datapoints = sum(n_datapoints))
      
      BA_excluded_by_specs_veryshort <- BA_excluded_by_specs %>% group_by(Trial_code) %>% 
        summarise(n_control = n_distinct(control),
                  n_test_type = n_distinct(test_type),
                  n_insecticide = n_distinct(insecticide),
                  n_strain = n_distinct(strain),
                  n_site = n_distinct(site), 
                  n_year = n_distinct(year), 
                  n_conc_perc = n_distinct(conc_perc), 
                  n_conc_mug = n_distinct(conc_mug), 
                  n_is_disc_dose = n_distinct(is_disc_dose), 
                  n_times_disc_dose = n_distinct(times_disc_dose), 
                  n_species = n_distinct(species), 
                  n_datapoints = sum(n_datapoints))
      
      
      print("The following BA data was excluded by run specification (summarised with n_somevariable meaning 'number of distinct values for variable 'somevariable').")
      print.data.frame(BA_excluded_by_specs_veryshort)
     
      # sumup species if sumup_species set to TRUE (and if no grouping by species)
      # WARNING: If there are separate data points that cannot be distinguished based on covars,
      # then this operation pools data points unintentionally!
      # AND FOR BA DATA THIS IS LIKELY TO BE THE CASE!
      if ( !('species' %in% group_by_vars) & sumup_species) {
        BA <- BA %>%
          group_by(across(setdiff(covars_BA, "species"))) %>%
          summarise(across(all_of(outcomes_BA), ~ sum(.x))) |>
          ungroup()
      }
      
    
    ## align site names across BA and EHT
      # if necessary mutate site names on EHT side
      EHT <- EHT %>%
       mutate(site = case_when(site == "Mwagaggala" ~ "Mwagagala",
                               TRUE ~ site))
      
      # check for remaining ambiguities in site names
      # gives list of similar site names according to specified similarity score 
      # by the similarity metric 'Optimal string aligment, (restricted Damerau-Levenshtein distance)'
      BA_site_names <- BA %>% 
        select(site, country) %>% 
        unique() 
      EHT_site_names <- EHT %>% 
        select(site, country) %>% 
        unique() 
      
      similarity_index <- 0.8
      
      all_site_names_simmatrix <- stringdist::stringsimmatrix(BA_site_names$site, EHT_site_names$site) # compute pairwise similarity score
      all_site_names_simmatrix <- all_site_names_simmatrix > similarity_index & all_site_names_simmatrix != 1
      all_site_names_sim <- which(all_site_names_simmatrix, arr.ind = TRUE)
      all_site_names_simpair <- tibble( cbind(
        cbind( BA_site_names$site[all_site_names_sim[ ,1]], BA_site_names$country[all_site_names_sim[ ,1]]), 
        cbind( EHT_site_names$site[all_site_names_sim[ ,2]], EHT_site_names$country[all_site_names_sim[ ,2]])))
      
      if (nrow(all_site_names_simpair) != 0){
        writeLines(paste0("WARNING: There are remaining ambiguities in the site names between BA and EHT side."))
        print.data.frame(all_site_names_simpair)
        writeLines("\n")
      }
      
      
    ## align country names across BA and EHT
      # if necessary mutate country names on EHT side
      EHT <- EHT %>%
        mutate(country = case_when(country == "Burkina Faso" ~ "BurkinaFaso",
                                TRUE ~ country))
      
      # check for remaining ambiguities in country names
      # gives list of similar country names according to specified similarity score 
      # by the similarity metric 'Optimal string aligment, (restricted Damerau-Levenshtein distance)'
      BA_countrys <- unique(BA$country)
      EHT_countrys <- unique(EHT$country)
      similarity_index <- 0.8
      
      all_country_names_simmatrix <- stringdist::stringsimmatrix(BA_countrys, EHT_countrys) # compute pairwise similarity score
      all_country_names_simmatrix <- all_country_names_simmatrix > similarity_index & all_country_names_simmatrix != 1
      all_country_names_sim <- which(all_country_names_simmatrix, arr.ind = TRUE)
      all_country_names_simpair <- tibble(cbind(BA_countrys[all_country_names_sim[ ,1]], EHT_countrys[all_country_names_sim[ ,2]]))
      
      if (nrow(all_country_names_simpair) != 0){
        writeLines(paste0("WARNING: There are remaining ambiguities in the country names between BA and EHT side."))
        print.data.frame(all_country_names_simpair)
        writeLines("\n")
      }
    
      
    ## align insecticide names across BA and EHT
      # if necessary mutate insecticide names on EHT side
      # EHT <- EHT %>%
      #   mutate(insecticide = case_when(insecticide == "" ~ "",
      #                           TRUE ~ insecticide))
      
      # check for remaining ambiguities in insecticide names
      # gives list of similar insecticide names according to specified similarity score 
      # by the similarity metric 'Optimal string aligment, (restricted Damerau-Levenshtein distance)'
      BA_insecticides <- unique(BA$insecticide)
      EHT_insecticides <- unique(EHT$insecticide)
      similarity_index <- 0.8
      
      all_insecticide_names_simmatrix <- stringdist::stringsimmatrix(BA_insecticides, EHT_insecticides) # compute pairwise similarity score
      all_insecticide_names_simmatrix <- all_insecticide_names_simmatrix > similarity_index & all_insecticide_names_simmatrix != 1
      all_insecticide_names_sim <- which(all_insecticide_names_simmatrix, arr.ind = TRUE)
      all_insecticide_names_simpair <- tibble(cbind(BA_insecticides[all_insecticide_names_sim[ ,1]], EHT_insecticides[all_insecticide_names_sim[ ,2]]))
      
      if (nrow(all_insecticide_names_simpair) != 0){
        writeLines(paste0("WARNING: There are remaining ambiguities in the insecticide names between BA and EHT side."))
        print.data.frame(all_insecticide_names_simpair)
        writeLines("\n")
      }
      
    
    
    ## MATCHING BA AND EHT PAIRS (THIS EXCLUDES CONTROLS, added again later if wanted)
    ## strategy is to assign unique group IDs for each BA and EHT data point
    ## 'group' generally refers to df containing defining group characteristics, 'group_number' refers to a column of IDs in that df
    ## NOTE: It is needed to keep separate variables BA_groups_matched and EHT_groups_matched,
    ## as the manual matches match pairs with non-identical characteristics.
      # matching by group_by_vars input

   # if only int dose data sets are to be selected, keep only BA_groups for which any(times_disc_dose > 1) TRUE   
      if ("int" %in% bioassay_type & !("disc" %in% bioassay_type)){
        BA_groups <- BA |> 
          filter(!control) %>%
          group_by(across(all_of(group_by_vars))) |> 
          filter(any(times_disc_dose > 1)) |>
          select(all_of(group_by_vars)) %>%
          unique()
      } else {
        BA_groups <- BA %>%
          filter(!control) %>%
          select(all_of(group_by_vars)) %>%
          unique()
      }
    
      EHT_groups <- EHT %>%
        filter(!control) %>%
        select(all_of(group_by_vars)) %>%
        unique()
      
      groups_common <- intersect(BA_groups, EHT_groups) %>%
        ungroup() %>%
        mutate(group_number = row_number())
      
      BA_groups_matched <- groups_common
      BA_groups_unmatched <- dplyr::setdiff(BA_groups, select(BA_groups_matched, !group_number)) 
      
      EHT_groups_matched <- groups_common
      EHT_groups_unmatched <- dplyr::setdiff(EHT_groups, select(EHT_groups_matched, !group_number))
        
      # add manual matches
      for (kk in 1: length(manual_match_functions)){
        
        f <- get(manual_match_functions[kk])
        # this sends the return values of the manualmatch function to the global environment,
        # being BA_groups_newmatch, EHT_groups_newmatch, BA_groups_unmatched, EHT_groups_unmatched
        list2env(f(BA_groups = BA_groups_unmatched, EHT_groups = EHT_groups_unmatched, N_groups_previous = max(BA_groups_matched$group_number)), envir = .GlobalEnv)
        
        # append matched groups
        BA_groups_matched <- rbind(BA_groups_matched, BA_groups_newmatch)
        EHT_groups_matched <- rbind(EHT_groups_matched, EHT_groups_newmatch)
      }
      
      print("The matched groups on the BA side are:")
      print.data.frame(BA_groups_matched  %>% arrange(Trial_code))
      print("The matched groups on the EHT side are:")
      print.data.frame(EHT_groups_matched  %>% arrange(Trial_code))
      
      # print non-matched and thus excluded groups to output
      BA_groups_still_unmatched <- dplyr::setdiff(BA_groups, select(BA_groups_matched, !group_number)) %>% arrange(Trial_code)
      EHT_groups_still_unmatched <- dplyr::setdiff(EHT_groups, select(EHT_groups_matched, !group_number)) %>% arrange(Trial_code)
      
      print("The following BA data groups remained unmatched and were therefore excluded from the processed data.")
      print.data.frame(BA_groups_still_unmatched)
      
      print("The following EHT data groups remained unmatched and were therefore excluded from the processed data.")
      print.data.frame(EHT_groups_still_unmatched)
      
      # add group numbers to df
      BA_matched <- BA %>%
        filter(!control) %>%
        inner_join(BA_groups_matched, by = group_by_vars)
      
      EHT_matched <- EHT %>%
        filter(!control) %>%
        inner_join(EHT_groups_matched, by = group_by_vars)
      
      
    ## ADD CONTROLS AGAIN IF WANTED
      # BE CAREFULL, THIS DUPLICATES THE CONTROL DATA IF THERE ARE MULTIPLE GROUPS PER Trial_code-site PAIR
      # SINCE IT ADDS THE CONTROL TO EACH GROUP OF A GIVEN Trial_code-site PAIR.
      ## Hmm, what group_number do the controls get?
      if (include_controls) {
        BA_groups_controls <- BA_groups_matched %>% 
        select(!insecticide) %>%
          unique()
        
        EHT_groups_controls <- EHT_groups_matched %>% 
          select(!insecticide) %>%
          unique()
  
        BA_control <- BA %>%
          filter(control) %>%
          inner_join(BA_groups_controls, by = setdiff(group_by_vars, "insecticide"))
        BA_matched <- bind_rows(BA_matched, BA_control)
        
        EHT_control <- EHT %>%
          filter(control) %>%
          inner_join(EHT_groups_controls, by = setdiff(group_by_vars, "insecticide"))
        EHT_matched <- bind_rows(EHT_matched, EHT_control)
        
        # groups with control
        BA_groups_withcontrols <- semi_join(BA_groups_controls, BA_control, by = setdiff(group_by_vars, "insecticide"))
        EHT_groups_withcontrols <- semi_join(EHT_groups_controls, EHT_control, by = setdiff(group_by_vars, "insecticide"))
      }
    
      
    # Data summaries
      print("Number of BA data points per group:")
      print.data.frame(BA_matched %>% group_by(group_number) %>% summarise(number_of_datapoints = n()))
      print("Number of EHT data points per group:")
      print.data.frame(EHT_matched %>% group_by(group_number) %>% summarise(number_of_datapoints = n()))
      
    
    ## prepare data for stan (select variables)
      # rename 
      EHT_matched <- EHT_matched %>% 
        rename(N_h = total,
               DF = bf_dead,
               DU = unf_dead,
               AF = bf_live,
               AU = unf_live,
               D_h = tot_dead,
               A = tot_live,
               U = tot_unf,
               F = tot_blf)  
      
      # rename
      BA_matched <- BA_matched %>% 
        rename(N_b = BA_total,
               D_b = BA_total_mort)
      
    # BA data
      # with all covariates
      T_bint_indices <- BA_matched |> filter(times_disc_dose > 1) |> pull(group_number) |> unique() |> sort()
      T_bdisc_indices <- sort( setdiff( unique(BA_matched$group_number), T_bint_indices ) )
      
      B_all <- BA_matched %>% 
        mutate(treat_b = treat,
               T_b = group_number,
               int_dose_available = as.integer(T_b %in% T_bint_indices)) 
      
      # only covariates relevant for stan
      B <- B_all %>%
        select(Trial_code, T_b, treat_b, N_b, D_b, country, site, year, insecticide, times_disc_dose, int_dose_available)
    
    # EHT data, all data points, i.e. at least part breakdown 
      # with all covariates
      H_all <- EHT_matched %>% 
        mutate(treat_h = treat,
               T_h = group_number,
               hut_type_h = as.numeric(factor(hut_type))) 
      
      # only covariates relevant for stan
      H <- H_all %>%
        select(Trial_code, T_h, treat_h, N_h, D_h)
      
    # EHT data, only data points with full breakdown
      # with all covariates
      H_f_all <- EHT_matched %>% 
        mutate(treat_h_f = treat,
               T_h_f = group_number,
               hut_type_h_f = as.numeric(factor(hut_type))) %>%
        mutate(H_ID = row_number()) %>% # adds index of same replicate in H
        drop_na(A, AF)
      
      # only covariates relevant for stan  
      H_f <- H_f_all %>%  
        select(Trial_code, T_h_f, H_ID, treat_h_f, A, AF)
  
    
    # number of data sets 
    nT <- length(unique(H$T_h))
    nT_f <- length(unique(H_f$T_h_f))
    nT_bint <- length(T_bint_indices)
    
    # number of hut types in H data set
    n_hut_type <- length(unique(H_all$hut_type))
    
    trials <- list(nT = nT, nT_f = nT_f, T_f_indices = sort(unique(H_f_all$group_number)), nT_bint = nT_bint, T_bint_indices = T_bint_indices, T_bdisc_indices = T_bdisc_indices, n_hut_type = n_hut_type)
  
  sink()
  
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
  
  
  ## save processed data as rds in respective run folder
  #############################
  saveRDS(B, file.path("fitting", run, "B.rds"))
  saveRDS(B_all, file.path("fitting", run, "B_all.rds"))
  
  if (!BA_only){
    saveRDS(H, file.path("fitting", run, "H.rds"))
    saveRDS(H_all, file.path("fitting", run, "H_all.rds"))
  }
  
  if (with_feeding){
    saveRDS(H_f, file.path("fitting", run, "H_f.rds"))
    saveRDS(H_f_all, file.path("fitting", run, "H_f_all.rds"))
  }
  
  saveRDS(trials, file.path("fitting", run, "trials.rds"))
  
  
  # randomly select groups to perform leave-one-out cross-validation
  #############################
  # weight the groups by the sample size, i.e. the sum of the totals, in the EHT treatment data, while summing up N_h and A.
  # I rely only on the EHT so that I can pick the same 5 groups for both disc and int version.
  
  H_sum <- H |>
    filter(treat_h == 1) |>
    group_by(T_h) |>
    summarise(N_h = sum(N_h))
  
  H_f_sum <- H_f |>
    filter(treat_h_f == 1) |>
    group_by(T_h_f) |>
    summarise(A = sum(A))
  
  EHT_treat_sum <- left_join(H_sum, H_f_sum, by = join_by(T_h == T_h_f)) |>
    mutate(total = N_h + if_else(is.na(A), 0, A))
  
  if (setequal(CV_samples_sizes, c(5,8))){
    if (with_feeding){
      set.seed(111)
      group_CVsample_5 <- sort(sample(EHT_treat_sum$T_h, 5, replace = F, prob = EHT_treat_sum$total))
      group_CVsample_8 <- sort(sample(EHT_treat_sum$T_h, 8, replace = F, prob = EHT_treat_sum$total))
    } else {
      set.seed(111)
      group_CVsample_5 <- sort(sample(EHT_treat_sum$T_h, 5, replace = F, prob = EHT_treat_sum$N_h))
      group_CVsample_8 <- sort(sample(EHT_treat_sum$T_h, 8, replace = F, prob = EHT_treat_sum$N_h))
    }
    saveRDS(group_CVsample_5, file = file.path("fitting", run,  paste0("group_CVsample_5",".rds")))
    saveRDS(group_CVsample_8, file = file.path("fitting", run,  paste0("group_CVsample_8",".rds")))
  }
} # ends if (!is.na(take_data_from_run))
