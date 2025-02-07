BA_data_processing <- function(data_sets, file_list_data_sets){
  # data_sets needs to be an array of integers
  # file_list_data_sets is a string containing the path of the xlsx containing information on data sets, including where the files are
  
  library(readxl)
  library(tidyverse)
  library(reshape2)
  library(lubridate)
  library(stringr)
  
  list_data_sets <- readxl::read_xlsx(file_list_data_sets, trim_ws = FALSE)

  remaining_data_sets <- data_sets
  
  covars_BA <- c("Trial_code", "control", "treat",
                 "test_type", "insecticide", "single_pyrethroid", "BA_total_known", 
                 "strain", "country", "site", "tube", "date", "year",
                 "species",
                 "conc_perc", "conc_mug", "is_disc_dose", "times_disc_dose"
                 )
  outcomes_BA <- c("BA_total", 
               "kd15", "kd30", "kd45", "kd60", 
               "mort_24h", "mort_48h", "mort_72h",
               "BA_total_mort", "BA_perc_mort"
               )
  
  vars_BA <- c(covars_BA, outcomes_BA)
  
  
  # start with disc BA data from Rebecca's data set
  B_df<- readxl::read_xlsx(file.path("data/data_discBA", "Database_A1_with unpublished_split.xlsx"),
                                sheet = "Data") %>% 
  filter(Trial_code %in% remaining_data_sets) %>%
    rename(test_type = Bioassay_type,
           insecticide = Bioassay_insecticide,
           BA_total = Bioassay_N_tested,
           BA_total_mort = Bioassay_N_dead,
           BA_perc_mort = Bioassay_mort,
           country = Country,
           site = Site,
           year = Est_Start) %>%
    filter( !(is.na(BA_total_mort) & is.na(BA_perc_mort)) ) %>% # removes empty lines that remained from removing BA duplicates 
    mutate(control = FALSE,
           conc_perc = as.numeric(conc),
           conc_mug = NA,
           strain = NA,
           site = str_remove_all(site, "\\_(.*)"),
           is_disc_dose = TRUE,
           BA_total_known = as.logical(Bioassay_N_known),
           KD = NA, 
           mort_24h = BA_total_mort, # I assume this was always 24h!
           mort_48h = NA,
           mort_72h = NA) %>%
    select(any_of(vars_BA))

  
  #additional information found in articles
  B_df <- B_df %>%
    mutate(test_type = case_when(Trial_code == 26 ~ "WHO", # found in article
                                 TRUE ~ test_type)
           )
  
  # data corrections
  B_df <- B_df %>%
    mutate(conc_perc = case_when(Trial_code == 32 & insecticide == "permethrin" ~ 0.75 , # assumed data entry error
                                 TRUE ~ conc_perc)
    )
  
  remaining_data_sets <- setdiff(remaining_data_sets, B_df$Trial_code)
  
  ## Intensity BA EHT data
  
  # Ngufor data 
  Trial_code_i <- 201
  if (Trial_code_i %in% remaining_data_sets){
    B_df_add <- readxl::read_xlsx(file.path(select(filter(list_data_sets, Trial_code==Trial_code_i), path2BAfile), 
                                            select(filter(list_data_sets, Trial_code==Trial_code_i), BA_file)),
                                  sheet = paste0(select(filter(list_data_sets, Trial_code==Trial_code_i), BA_sheet))) %>% 
      rename(insecticide = Treatment,
             BA_total = N,
             strain = Strain,
             date = Date,
             kd15 = "15", 
             kd30 = "30", 
             kd45 = "45", 
             kd60 = "60", 
             mort_24h = "24 h",
             mort_48h = "48 h",
             mort_72h = "72 h") %>%
      mutate(Trial_code = Trial_code_i,
             test_type = "CDC", # check this!!!
             control = str_detect(insecticide, "Control|control"), 
             BA_total_known = TRUE,
             BA_total_mort = pmax(mort_24h, mort_48h, mort_72h, na.rm = TRUE),
             BA_perc_mort = NA,
             country = "Benin",
             site = "Cove", # check this!!!
             conc_perc = NA,
             conc_mug = as.numeric(str_remove(Dose, "µg"))
             ) %>%
      select(any_of(vars_BA))
    
    B_df <- B_df %>% bind_rows(B_df_add)
    rm(B_df_add)
    remaining_data_sets <- setdiff(remaining_data_sets, B_df$Trial_code)
  }
  
  Trial_code_i <- 202
  if (Trial_code_i %in% remaining_data_sets){
    B_df_add <- readxl::read_xlsx(file.path(select(filter(list_data_sets, Trial_code==Trial_code_i), path2BAfile), 
                                            select(filter(list_data_sets, Trial_code==Trial_code_i), BA_file)),
                                  sheet = paste0(select(filter(list_data_sets, Trial_code==Trial_code_i), BA_sheet))) %>% 
      rename(insecticide = "Treatment",
             BA_total = "N",
             strain = Strain,
             date = Date,
             kd15 = "15", 
             kd30 = "30", 
             kd45 = "45", 
             kd60 = "60",
             mort_24h = "24 h",
             mort_48h = "48 h",
             mort_72h = "72 h") %>%
      mutate(Trial_code = Trial_code_i,
             test_type = "CDC", # check this!!!
             control = str_detect(insecticide, "Control|control"), 
             BA_total_known = TRUE,
             BA_total_mort = pmax(mort_24h, mort_48h, mort_72h, na.rm = TRUE),
             BA_perc_mort = NA,
             country = "Benin",
             site = "Cove", # check this!!!
             conc_perc = NA,
             conc_mug = as.numeric(str_remove(Dose, "µg"))
             ) %>%
      select(any_of(vars_BA))
    
    B_df <- B_df %>% bind_rows(B_df_add)
    rm(B_df_add)
  }  

  
  # Antoine Sanou's data
  # maybe consistency check Dea_24h + Alive_24h = Nb_tested?
  Trial_code_i <- 203
  if (Trial_code_i %in% remaining_data_sets){
    B_df_add<- readxl::read_xlsx(file.path(select(filter(list_data_sets, Trial_code==Trial_code_i), path2BAfile), 
                                            select(filter(list_data_sets, Trial_code==Trial_code_i), BA_file)),
                                  sheet = paste0(select(filter(list_data_sets, Trial_code==Trial_code_i), BA_sheet))) %>% 
      rename(insecticide = insecticide,
             BA_total = Nb_tested,
             species = Species,
             date = Date_test,
             conc_perc = Concentration,
             kd60 = KD60, 
             mort_24h = "Dea_24h",
             site = Village) %>%
      mutate(Trial_code = Trial_code_i,
             test_type = "WHO", # check this!!!
             conc_mug = NA,
             control = str_detect(insecticide, "Control|control"), 
             BA_total_known = if_else(!is.na(BA_total), TRUE, FALSE),
             BA_total_mort = mort_24h,
             BA_perc_mort = NA,
             country = "BurkinaFaso",
             mort_48h = NA,
             mort_72h = NA,
             year = year(date)) %>%  # Year in data was sometimes missing, maybe check consistency
      select(any_of(vars_BA)) 
    
    # exclude Lab data
    B_df_add <- B_df_add %>%
      filter(site != "Lab")
    
    B_df <- B_df %>% bind_rows(B_df_add)
    rm(B_df_add)
    remaining_data_sets <- setdiff(remaining_data_sets, B_df$Trial_code)
  }
  
  
  # Hyacinthe Toe's data
  # maybe consistency check Dead + Aliveh = Total tested?
  if (204 %in% remaining_data_sets | 205 %in% remaining_data_sets){
    B_df_add<- readxl::read_xlsx(file.path(select(filter(list_data_sets, Trial_code==204), path2BAfile), # note: Trial_code 204 and 205 in same excel sheet
                                           select(filter(list_data_sets, Trial_code==204), BA_file)),
                                 sheet = paste0(select(filter(list_data_sets, Trial_code==204), BA_sheet)),
                                 range = cell_cols("A:J")) %>%  # note: this excludes column with summaries over replicates
      rename(insecticide = Insecticide,
             BA_total = "Total tested",
             species = Specie,
             date = Date,
             strain = Strain,
             tube = "Tube N",
             kd60 = "KD 60", 
             mort_24h = "Dead", # check it really is 24h!!!
             site = Strain
             ) %>% 
      mutate(Trial_code = case_when(str_detect(site, "Tengrela") ~ 204,
                                    str_detect(site, "VK5") ~ 205),
             test_type = "WHO", # check!!
             is_disc_dose = TRUE, # I assumed disc dose
             conc_perc = as.numeric(Dosage),
             control = FALSE, # check this, especially for VK5 where I don't have the concentrations
             BA_total_known = if_else(!is.na(BA_total), TRUE, FALSE),
             BA_total_mort = mort_24h,
             inconsistent = Alive + BA_total_mort != BA_total,
             BA_perc_mort = NA,
             country = "BurkinaFaso",
             mort_48h = NA,
             mort_72h = NA,
             year = year(date)) %>%
      select(any_of(vars_BA)) 
    
    # reduce in case either only 204 or only 205 was in data_sets
    if (!(204 %in% remaining_data_sets)){
      B_df_add <- B_df_add %>% filter(Trial_code != 204)
    } 
    if (!(205 %in% remaining_data_sets)){
      B_df_add <- B_df_add %>% filter(Trial_code != 205)
    } 
    
    B_df <- B_df %>% bind_rows(B_df_add)
    rm(B_df_add)
    remaining_data_sets <- setdiff(remaining_data_sets, B_df$Trial_code)
  }
  
  
  # Emidi's data from Tanzania
  Trial_code_i <- 208
  if (Trial_code_i %in% remaining_data_sets){
      B_df_add<- readxl::read_xlsx(file.path(select(filter(list_data_sets, Trial_code==Trial_code_i), path2BAfile), 
                                             select(filter(list_data_sets, Trial_code==Trial_code_i), BA_file)),
                                   sheet = paste0(select(filter(list_data_sets, Trial_code==Trial_code_i), BA_sheet)),
                                   range = cell_cols("A:O"),  # note: this excludes column with summaries over replicates
                                   col_types = c("date", rep("guess", 14))) %>% 
        rename(BA_total = `Total exposed`,
               kd60 = `1 hr KD`,
               conc_perc = Dose,
               date = Date,
               site = Site,
               mort_24h = `24h Dead`
               ) %>%
        filter(!is.na(mort_24h) & !is.na(BA_total)) %>% # remove empty rows
        mutate(insecticide = Treatment,
               Trial_code = Trial_code_i,
               test_type = case_when(str_detect(`Expo. Method`, regex("WHO", ignore_case = TRUE)) ~ "WHO"),
               conc_mug = NA,
               conc_perc = as.numeric(conc_perc) * 100, # guessed so 
               control = str_detect(Treatment, "Control|control"), 
               BA_total_known = if_else(!is.na(BA_total), TRUE, FALSE),
               BA_perc_mort = NA,
               country = "Tanzania",
               mort_48h = NA,
               mort_72h = NA,
               year = 2021) %>%
        select(any_of(vars_BA)) 
      
      # data correction
      B_df_add$site[c(35,39)] <- "Igunga"
      B_df_add <- B_df_add %>% 
        mutate(insecticide = case_when(str_detect(insecticide, regex("delta", ignore_case = TRUE)) ~ "deltamethrin",
                                       TRUE ~ insecticide))
    
    B_df <- B_df %>% bind_rows(B_df_add)
    rm(B_df_add)
    remaining_data_sets <- setdiff(remaining_data_sets, B_df$Trial_code)
  }
  
  
  # Luc Djogbenou's data
  Trial_code_i <- 206
  if (Trial_code_i %in% remaining_data_sets){
    B_df_add<- readxl::read_xlsx(file.path(select(filter(list_data_sets, Trial_code==Trial_code_i), path2BAfile), 
                                           select(filter(list_data_sets, Trial_code==Trial_code_i), BA_file)),
                                 sheet = paste0(select(filter(list_data_sets, Trial_code==Trial_code_i), BA_sheet))) %>% 
      rename(species = Strain_Name,
             date = Test_Date,
             country = Test_Country,
             site = Test_Site,
             tube = Tube_number
             ) %>%
      mutate(Trial_code = Trial_code_i,
             insecticide = Treatment_Name,
             test_type = case_when(str_detect(Test_Name, "WHO") ~ "WHO"),
             conc_mug = NA,
             conc_perc = as.numeric(Treatment_dose),
             control = str_detect(Treatment_Name, "Untreated"), 
             year = year(date)
             )  
      
    B_df_add <- B_df_add %>% 
      group_by(across(intersect(covars_BA, names(B_df_add)))) %>%
      summarise(  
        BA_total = n(),
        mort_24h = sum(`24h_dead`),
        kd60 = sum(`1h_KD`)
      ) %>%
      mutate(
        mort_48h = NA,
        mort_72h = NA,          
        BA_total_mort = mort_24h,
        BA_total_known = TRUE,
        BA_perc_mort = NA,
      ) %>% 
      select(any_of(vars_BA)) 
    
    # data correction, probably to be changed later by new data version
    B_df_add <- B_df_add %>%
      mutate(conc_perc = case_when(insecticide == "Deltamethrine" ~ 0.05, # check!!
                                   insecticide == "Permethrine" ~ 0.75, # check!!
                                   TRUE ~ conc_perc))
    
    B_df <- B_df %>% bind_rows(B_df_add)
    rm(B_df_add)
    remaining_data_sets <- setdiff(remaining_data_sets, B_df$Trial_code)
  }
  
  
  # Protopopoff's data 
  Trial_code_i <- 209
  if (Trial_code_i %in% remaining_data_sets){
    B_df_add<- readxl::read_xlsx(file.path(select(filter(list_data_sets, Trial_code==Trial_code_i), path2BAfile), 
                                           select(filter(list_data_sets, Trial_code==Trial_code_i), BA_file)),
                                 sheet = paste0(select(filter(list_data_sets, Trial_code==Trial_code_i), BA_sheet))) %>% 
      rename(insecticide = Insecticide_name,
             species = Species,
             times_disc_dose = Concetration, # sic!
             BA_total = Tot_mosq_exposed,
             mort_24h = `24hrsD`,
             mort_48h = `48hrsD`,
             mort_72h = `72hrsD` 
             ) %>%
      mutate(Trial_code = Trial_code_i,
             insecticide = case_when(str_detect(insecticide, "PPF") ~ "Pyriproxyfen",
                                        TRUE ~ insecticide),
             test_type = "CDC", # according to email by Jackline Martin 21.04.2022
             tube = as.numeric(Bottle),
             conc_mug = NA,
             conc_perc = NA,
             control = case_when(str_detect(insecticide, regex("Acetone", ignore_case = TRUE)) ~ TRUE,  # according to email by Jackline Martin 21.04.2022
                                 TRUE ~ FALSE), 
             times_disc_dose = case_when(control ~ 0,
                                         TRUE ~ times_disc_dose),
             BA_total_known = if_else(!is.na(BA_total), TRUE, FALSE),
             BA_perc_mort = NA,
             country = "Tanzania",
             site = "Welamasonga",  # according to email by Jackline Martin 21.04.2022
             BA_total_mort = mort_24h + mort_48h + mort_72h,
             year = case_when(`Study year` == 1 ~ 2020, # from email with Jackline Martin
                              `Study year` == 2 ~ 2021, # from email with Jackline Martin
                              `Study year` == 3 ~ 2022) # from email with Jackline Martin
             ) %>%
      select(any_of(vars_BA)) %>%
      filter(species != 'kisumu' ) # lab strain, email with Jackline Martin 2022-11-15
    
    B_df <- B_df %>% bind_rows(B_df_add)
    rm(B_df_add)
    remaining_data_sets <- setdiff(remaining_data_sets, B_df$Trial_code)
  }
  
  
  
  
  ## after all data sets have been included
  #########################################
  
  # classify and unify insecticides
  pyrethroids <- c("alphacypermethrin", "lambdacyhalothrin", "permethrin", "pyrethroid", "deltamethrin", "bifenthrin")  
  B_df <- B_df %>%
    mutate(
      insecticide = str_to_lower(insecticide),
      insecticide = case_when(control ~ "none", 
                                 TRUE ~ insecticide),
      single_pyrethroid = case_when(str_detect(insecticide, regex("pirimiphos-methyl|bendiocarb|fenitrothion|pyriproxyfen|chlorfenapyr|pbo|\\+", ignore_case = TRUE) ) ~ FALSE,
                                    str_detect(insecticide, regex(paste0(pyrethroids, collapse = "|"), ignore_case = TRUE) ) ~ TRUE),
      insecticide = case_when(single_pyrethroid & str_detect(insecticide, "alpha(|-)cyper") ~ "alphacypermethrin",
                                single_pyrethroid & str_detect(insecticide, "lambda(|-)cyhalothrin") ~ "lambdacyhalothrin",
                                single_pyrethroid & str_detect(insecticide, "deltamethrin") ~ "deltamethrin",
                                single_pyrethroid & str_detect(insecticide, "permethrin") ~ "permethrin",
                                TRUE ~ insecticide),
      #manual corrections
      single_pyrethroid = case_when(str_detect(insecticide, "(pbo synergism)") ~ TRUE,
                                    str_detect(insecticide, "alphacypermthrine") ~ TRUE,
                                    TRUE ~ single_pyrethroid),
      insecticide = case_when(str_detect(insecticide, "(pbo synergism)") ~ "deltamethrin",
                                 str_detect(insecticide, "alphacypermthrine") ~ "alphacypermethrin",
                                 TRUE ~ insecticide)
    )
  
  

  # compute missing concentration information
  # if variable times_disc_dose not there yet
  if (!("times_disc_dose" %in% names(B_df))){
    B_df <- B_df |>
      mutate(times_disc_dose = NA)
  }
  B_df <- B_df %>%
    mutate(is_disc_dose = if_else(is.na(as.logical(is_disc_dose)), 
                             case_when(control ~ NA,
                                       insecticide == "alphacypermethrin" & (conc_perc == 0.05 | conc_mug == 12.5 | times_disc_dose == 1) ~ TRUE,
                                       insecticide == "lambdacyhalothrin" & (conc_perc == 0.05 | conc_mug == 12.5 | times_disc_dose == 1) ~ TRUE,
                                       insecticide == "deltamethrin" & (conc_perc == 0.05 | conc_mug == 12.5 | times_disc_dose == 1) ~ TRUE,
                                       insecticide == "permethrin" & (conc_perc == 0.75 | conc_mug == 21.5 | times_disc_dose == 1) ~ TRUE,
                                       TRUE ~ FALSE
                                       ),
                             as.logical(is_disc_dose)),
           times_disc_dose = if_else(is.na(as.numeric(times_disc_dose)), 
                                     case_when(control ~ 0,
                                               is_disc_dose ~ 1,
                                               insecticide == "alphacypermethrin" & test_type == "WHO" ~ conc_perc / 0.05,
                                               insecticide == "alphacypermethrin" & test_type == "CDC" ~ conc_mug / 12.5,
                                               insecticide == "lambdacyhalothrin" & test_type == "WHO" ~ conc_perc / 0.05,
                                               insecticide == "lambdacyhalothrin" & test_type == "CDC" ~ conc_mug / 12.5,
                                               insecticide == "deltamethrin" & test_type == "WHO" ~ conc_perc / 0.05,
                                               insecticide == "deltamethrin" & test_type == "CDC" ~ conc_mug / 12.5,
                                               insecticide == "permethrin" & test_type == "WHO" ~ conc_perc / 0.75,
                                               insecticide == "permethrin" & test_type == "CDC" ~ conc_mug / 21.5,
                                               TRUE ~ NA_real_
                                     ),
                                     as.numeric(times_disc_dose)),
           #compute conc_perc if only times_disc_dose given
           conc_perc = if_else(is.na(conc_perc), 
                               case_when(control ~ 0, 
                                         test_type == "WHO" & insecticide == "alphacypermethrin" ~ times_disc_dose * 0.05,
                                         test_type == "WHO" & insecticide == "lambdacyhalothrin" ~ times_disc_dose * 0.05,
                                         test_type == "WHO" & insecticide == "deltamethrin" ~ times_disc_dose * 0.05,
                                         test_type == "WHO" & insecticide == "permethrin" ~ times_disc_dose * 0.75,
                                         TRUE ~ conc_perc),
                               conc_perc),
           #compute conc_mug if only times_disc_dose given
           conc_mug = if_else(is.na(conc_mug), 
                              case_when(control ~ 0, 
                                        test_type == "CDC" & insecticide == "alphacypermethrin" ~ times_disc_dose * 12.5,
                                        test_type == "CDC" & insecticide == "lambdacyhalothrin" ~ times_disc_dose * 12.5,
                                        test_type == "CDC" & insecticide == "deltamethrin" ~ times_disc_dose * 12.5,
                                        test_type == "CDC" & insecticide == "permethrin" ~ times_disc_dose * 21.5,
                                        TRUE ~ conc_mug),
                              conc_mug)
    )
  
  
  # add years
  B_df <- B_df %>%
    mutate(year = as.integer(year),
           year = if_else(is.na(year), 
                          as.integer(year(date)),
                          as.integer(year))
    )
  
  
  # unify covariates
  B_df <- B_df %>%
    mutate(test_type = case_when(str_detect(test_type, regex("WHO", ignore_case = TRUE)) ~ "WHO",
                                 str_detect(test_type, regex("CDC", ignore_case = TRUE)) ~ "CDC",
                                 Trial_code == 25 ~ "WHO", # guessed from the fact that Trial_code 25 is from a WHO report
                                 Trial_code == 36 ~ "WHO" , # guessed since concentration percentage was given
                                 TRUE ~ test_type),
           species = case_when(str_detect(species, regex("gamb(|.*)s(|.)s(|.)", ignore_case = TRUE)) ~ "gambiae ss",
                               str_detect(species, regex("gamb", ignore_case = TRUE)) ~ "gambiae sl",
                               str_detect(species, regex("fune", ignore_case = TRUE)) ~ "funestus",
                               TRUE ~ species)
    )
  
  
  # align country names
  B_df <- B_df %>%
    mutate(country = case_when(country == "Burkina Faso" ~ "BurkinaFaso",
                               TRUE ~ country) 
           )
  
  # check for remaining ambiguities in country names
  # gives list of similar site names according to specified similarity score 
  # by the similarity metric 'Optimal string aligment, (restricted Damerau-Levenshtein distance)'
  all_country_names <- unique(B_df$country) 
  similarity_index <- 0.8
  
  all_country_names_simmatrix <- stringdist::stringsimmatrix(all_country_names, all_country_names) # compute pairwise similarity score
  all_country_names_simmatrix[lower.tri(all_country_names_simmatrix, diag = TRUE)] <- 0 # set lower triangular matrix including diag to 0 
  all_country_names_simmatrix <- all_country_names_simmatrix > similarity_index
  all_country_names_sim <- which(all_country_names_simmatrix, arr.ind = TRUE)
  all_country_names_simpair <- tibble( cbind(all_country_names[all_country_names_sim[ ,1]], all_country_names[all_country_names_sim[ ,2]]))
  
  if (nrow(all_country_names_simpair) != 0){
    writeLines(paste0("WARNING: There are remaining ambiguities in the country names."))
    print.data.frame(all_country_names_simpair)
    writeLines("\n")
  }
  
  
  # align site names
  B_df <- B_df %>%
    mutate(site = case_when(str_detect(site, regex("soumousso", ignore_case = TRUE)) ~ "Soumousso",
                            TRUE ~ site))
  
  # check for remaining ambiguities in site names
  # gives list of similar site names according to specified similarity score 
  # by the similarity metric 'Optimal string aligment, (restricted Damerau-Levenshtein distance)'
  all_site_names <- B_df %>% 
    select(site, country) %>% 
    unique() 
  
  similarity_index <- 0.8
  
  all_site_names_edit <- str_replace_all(all_site_names$site, " village", "")
  all_site_names_simmatrix <- stringdist::stringsimmatrix(all_site_names_edit, all_site_names_edit) # compute pairwise similarity score
  all_site_names_simmatrix[lower.tri(all_site_names_simmatrix, diag = TRUE)] <- 0 # set lower triangular matrix including diag to 0 
  all_site_names_simmatrix <- all_site_names_simmatrix > similarity_index
  all_site_names_sim <- which(all_site_names_simmatrix, arr.ind = TRUE)
  all_site_names_simpair <- tibble( cbind(
    cbind( all_site_names$site[all_site_names_sim[ ,1]], all_site_names$country[all_site_names_sim[ ,1]]), 
    cbind( all_site_names$site[all_site_names_sim[ ,2]], all_site_names$country[all_site_names_sim[ ,2]])))
  
  if (nrow(all_site_names_simpair) != 0){
    writeLines(paste0("WARNING: There are remaining ambiguities in the site names."))
    print.data.frame(all_site_names_simpair)
    writeLines("\n")
  }
  
  # quick fix for BA_total_mort
  B_df <- B_df %>% 
    mutate(BA_total_mort = case_when(is.na(BA_total_mort) ~ mort_24h,
                                     TRUE ~ BA_total_mort)
    )
  
  ## omit inconsistent data points  
  B_df <- B_df %>% 
    mutate(data_inconsistent = case_when(conc_perc == 0 & !control ~ TRUE,
                                         conc_mug == 0 & !control ~ TRUE,
                                         BA_total_known & is.na(BA_total) ~ TRUE,
                                         TRUE ~ FALSE))

  if (TRUE %in% B_df$data_inconsistent){
    print("The following BA data contains inconsistent information and is removed.")
    print.data.frame(B_df %>% filter(data_inconsistent == TRUE))
  }
  
  B_df <- B_df %>%
    filter(!data_inconsistent)
  
  
  ## omit data points with missing key data
  B_df <- B_df %>% 
    mutate(data_missing = case_when(is.na(mort_24h) & is.na(mort_48h) & is.na(mort_72h) & is.na(BA_total_mort) & is.na(BA_perc_mort) ~ TRUE,
                                    is.na(BA_total_mort) ~ TRUE,
                                    is.na(BA_total) ~ TRUE,
                                    TRUE ~ FALSE))
  
  if (TRUE %in% B_df$data_missing){
    print("The following BA data points missed key data and are removed.")
    print.data.frame(B_df %>% filter(data_missing == TRUE))
  }
  
  B_df <- B_df %>%
    filter(!data_missing)
  

  # define new variable 'treat' from 'control'
  B_df <- B_df %>%
    mutate(treat = as.integer(!as.logical(control)))

  
  return <-list(BA = B_df, covars_BA = covars_BA, outcomes_BA = outcomes_BA, vars_BA = vars_BA)
}
