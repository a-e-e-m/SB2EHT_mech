EHT_data_processing <- function(data_sets, file_list_data_sets){
  # data_sets needs to be an array of integers
  # file_list_data_sets is the filename of the xlsx containing information on data sets, including where the files are
  
  library(readxl)
  library(tidyverse)
  library(reshape2)
  library(lubridate)
  library(stringr)
  library(haven)
  
  
  list_data_sets <- readxl::read_xlsx(file_list_data_sets, trim_ws = FALSE)
  
  # for Trial_code to ID conversion for complete breakdown data from Rebecca
    lookup_data_sets_complete_breakdown <- readxl::read_xlsx(file.path("data/data_discBA", "Complete_Breakdown.xlsx"),
                                                           sheet = "Studies") %>% 
    select(Trial_code, ID) %>%
    filter(Trial_code %in% data_sets)
  
  remaining_data_sets <- data_sets
  
  covars_EHT <- c("Trial_code", "treat_name", "control", "treat",
              "net_product", "net_gen", "insecticide", "net_status", "N_washes", "outcome_timepoint", 
              "species", "country", "site", "hut_type", "holes", 
              "date", "year",
              "round", "week", "day", "hut", "sleeper",
              "compartment")
  
  
  
  #WARNING: gravid coutns are dealt with in the end, don't add them to other categories in the section loading the individual data sets!
  # totals are computed in the end if not given
  outcomes_EHT <- c("total", "unf_live", "unf_dead", "bf_live", "bf_dead", "gravid_live", "gravid_dead",
                    "tot_unf", "tot_blf", "tot_live", "tot_dead")
  
  vars_EHT <- c(covars_EHT, outcomes_EHT)
  
  #initialise data frame
  H_df <- data.frame(
    # covars
      "Trial_code" = integer(), 
      "treat_name" = character(), 
      "control" = logical(), 
      "treat" = integer(),
      "net_product" = character(), 
      "net_gen" = character(), 
      "insecticide" = character(), 
      "net_status" = character(), 
      "N_washes" = integer(), 
      "outcome_timepoint" = character(), 
      "species" = character(), 
      "country" = character(), 
      "site" = character(), 
      "hut_type" = character(), 
      "holes" = character(), 
      "date" = as.Date(character()), 
      "year" = integer(),
      "round" = integer(), 
      "week" = integer(), 
      "day" = integer(), 
      "hut" = integer(), 
      "sleeper" = integer(),
      "compartment" = character(),
    # outcomes
      "total" = double(), 
      "unf_live" = double(), 
      "unf_dead" = double(), 
      "bf_live" = double(), 
      "bf_dead" = double(), 
      "gravid_live" = double(), 
      "gravid_dead" = double(),
      "tot_unf" = double(), 
      "tot_blf" = double(), 
      "tot_live" = double(), 
      "tot_dead" = double()
  )
  
  ## NOTE: Generally, net_gen, insecticide and N_washes are to be guessed from treat_name after loading all data sources
  
  ## start with complete breakdown disc BA data
  
    # Nets24 (may contain species, but actually not, so skipped that column)
    H_df_add <- readxl::read_xlsx(paste0("data/data_discBA/", "Complete_Breakdown.xlsx"),
                              sheet = "Nets24", col_types = c(rep("guess", 31), "text")) %>% 
      filter(ID %in% lookup_data_sets_complete_breakdown$ID) %>%
      left_join(lookup_data_sets_complete_breakdown, by = "ID") %>%
      rename(treat_name = treat) %>%    
      mutate(control = str_detect(treat_name, "control"),
             outcome_timepoint = "24h"
             ) %>%
      select(any_of(vars_EHT)) %>%  
      mutate(Trial_code = as.double(Trial_code),
             hut = as.double(hut),
             sleeper = as.double(sleeper))
  
    H_df <- H_df %>% bind_rows(H_df_add)
    rm(H_df_add)
    
    # Nets72
    # ??careful, these also need to be loaded as summaries (?) for 24h mortality 
    # ??hence, no remaining_data_set after that block.
    H_df_add <- readxl::read_xlsx(paste0("data/data_discBA/", "Complete_Breakdown.xlsx"),
                                  sheet = "Nets72") %>% 
      filter(ID %in% lookup_data_sets_complete_breakdown$ID) %>%
      left_join(lookup_data_sets_complete_breakdown, by = "ID") %>%
      rename(treat_name = treat) %>%    
      mutate(control = str_detect(treat_name, "control"),
             outcome_timepoint = "72h"
             ) %>%
      select(any_of(vars_EHT)) %>%  
      mutate(Trial_code = as.double(Trial_code),
             hut = as.double(hut),
             sleeper = as.double(sleeper))
    
    H_df <- H_df %>% bind_rows(H_df_add)
    rm(H_df_add)
    
    # load Rebecca's BA-EHT-matched data set to add covariates 
    compl_breakd_covar <- readxl::read_xlsx(file.path("data/data_discBA", "Database_A1_with unpublished_split.xlsx"),
                                            sheet = "Data") %>% 
      filter(Trial_code %in% remaining_data_sets) %>%
      rename(treat_name = Net,
             insecticide = Net_insecticide,
             country = Country, 
             site = Site,
             year = Est_Start, 
             hut_type = Hut, 
             species = Species) %>%
      filter( !(is.na(N_total) & is.na(N_dead)) ) %>%
      select(Trial_code, treat_name, insecticide, country, site, year, hut_type, species) %>%
      mutate(site = str_remove_all(site, "\\_(.*)"),) %>%
      unique()
    
    # add country, site, year, hut_type
    # I only join by Trial_code since I checked that country, site, year and hut_type 
    # are univariate per Trial_code
    H_df <- H_df %>%
      select(!c(country, site, year, hut_type)) %>% # I checked theu are all NA so far
      left_join(compl_breakd_covar %>% select(!c(species, treat_name, insecticide)) %>% unique(), by = c("Trial_code"))
    
    lookup_data_sets_complete_breakdown_remaining <- lookup_data_sets_complete_breakdown %>% # not sure what I need this for anymore. check!!
      filter(!(Trial_code %in% H_df$Trial_code))
    remaining_data_sets <- setdiff(remaining_data_sets, H_df %>% filter(outcome_timepoint == "24h") %>% pull(Trial_code) )
    
    #add species manually, according to Rebecca's table
    H_df <- H_df %>%
      mutate(species = case_when(Trial_code == 13 ~ "gambiae sl",
                                 Trial_code == 14 ~ "gambiae sl",
                                 TRUE ~ species))
  
  
  # EHT without full breakdown, treatment arms
  H_df_add <- readxl::read_xlsx(file.path("data/data_discBA", "Database_A1_with unpublished_split.xlsx"),
                                sheet = "Data") %>% 
    filter(Trial_code %in% remaining_data_sets) %>%
    mutate(control = FALSE,
           outcome_timepoint = "24h") %>%
    rename(treat_name = Net,
           total = N_total,
           tot_dead = N_dead,
           insecticide = Net_insecticide,
           net_gen = Inter_Type,
           hut_type = Hut,
           country = Country,
           site = Site,
           year = Est_Start,
           species = Species) %>%
    mutate(site = str_remove_all(site, "\\_(.*)"),
           N_washes = 0) %>%
    select(any_of(vars_EHT))  %>%
    filter( !is.na(tot_dead) )
  
  H_df <- H_df %>% bind_rows(H_df_add)
  rm(H_df_add)
  
  # EHT without full breakdown, control arms 
  H_df_add <- readxl::read_xlsx(file.path("data/data_discBA", "Database_A1_with unpublished_split.xlsx"),
                                sheet = "Data") %>% 
    filter(Trial_code %in% remaining_data_sets) %>%
    mutate(control = TRUE,
           outcome_timepoint = "24h") %>%
    rename(treat_name = Control_type,
           total = Control_total,
           tot_dead = Control_dead,
           hut_type = Hut,
           country = Country,
           site = Site,
           year = Est_Start,
           species = Species) %>%
    mutate(site = str_remove_all(site, "\\_(.*)"),) %>%
    select(any_of(vars_EHT)) %>%
    filter( !is.na(tot_dead) )
  
  H_df <- H_df %>% bind_rows(H_df_add)
  rm(H_df_add)
  remaining_data_sets <- setdiff(remaining_data_sets, H_df$Trial_code)
  
  
  # Olivier's data sets
  
  
  # Intensity BA EHT data
  # note: for each study / data set there is a separat code snippet in ./src/R/EHT_data_processing_snippets
  # each snippet appends its data (usually called H_df_add) to H_df and then should delete H_df_add and other intermediate files
  
  # Ngufor's data
  Trial_code_i <- 201
  if (Trial_code_i %in% remaining_data_sets){
    source('./src/R/EHT_data_processing_snippets/201_Ngufor.R', local=TRUE)
  }
  
  Trial_code_i <- 202
  if (Trial_code_i %in% remaining_data_sets){
    source('./src/R/EHT_data_processing_snippets/202_Ngufor.R', local=TRUE)
  }  
  
  # Antoine Sanou's data  
  Trial_code_i <- 203
  if (Trial_code_i %in% remaining_data_sets){
    source('./src/R/EHT_data_processing_snippets/203_Sanou.R', local=TRUE)
  }
  
  
  # Hyacinthe Toe's data  
  Trial_code_i <- 204 
  if (Trial_code_i %in% remaining_data_sets){
    source('./src/R/EHT_data_processing_snippets/204_Toe.R', local=TRUE)
  }
  
  Trial_code_i <- 205
  if (Trial_code_i %in% remaining_data_sets){
    source('./src/R/EHT_data_processing_snippets/205_Toe.R', local=TRUE)
  }
  
  
  # Basiliana Emidi's data from Tanzania
  Trial_code_i <- 208
  if (Trial_code_i %in% remaining_data_sets){
    source('./src/R/EHT_data_processing_snippets/208_Emidi.R', local=TRUE)
  }
  
  
  # Djogbenou's data from Benin
  Trial_code_i <- 206
  if (Trial_code_i %in% remaining_data_sets){
    source('./src/R/EHT_data_processing_snippets/206_Djogbenou.R', local=TRUE)
  }
  
  
  # Natasha Protopopoff's data
  Trial_code_i <- 209
  if (Trial_code_i %in% remaining_data_sets){
    source('./src/R/EHT_data_processing_snippets/209_Protopopoff.R', local=TRUE)
  }
  
  
  # Sarah Moore's data from Tiassale, Burkina Faso
  Trial_code_i <- 210
  if (Trial_code_i %in% remaining_data_sets){
    source('./src/R/EHT_data_processing_snippets/210_Moore.R', local=TRUE)
  }
  
  
  
  ## after all data sets have been included
  #########################################
  
  # fill in missing variables based on information given in treat_name
  # NOTE: case_when evaluates with conditions in the given order, 
  # hence substrings of other strings need to come later,
  # e.g. first 'Interceptor G2' and only then 'Interceptor'
  H_df <- H_df %>%
    mutate(net_product = if_else(is.na(net_product), 
                                    case_when(str_detect(treat_name, "Interceptor G2|IntG2") ~ "InterceptorG2", 
                                              str_detect(treat_name, "Interceptor|Int") ~ "Interceptor",
                                              str_detect(treat_name, "Royal Guard") ~ "RoyalGuard",
                                              str_detect(treat_name, "PermaNet 3.0|Perm3.0") ~ "PermaNet3",
                                              str_detect(treat_name, "PermaNet 2.0|Perm2.0") ~ "PermaNet2",
                                              str_detect(treat_name, "DawaPlus2.0") ~ "DawaPlus2",
                                              str_detect(treat_name, "OlysetPlus") ~ "OlysetPlus",
                                              str_detect(treat_name, "Olyset") ~ "Olyset",
                                              str_detect(treat_name, "chlorfenapyr_CTN") ~ "chlorfenapyr_CTN",
                                              str_detect(treat_name, "alphacy_CTN") ~ "alphacy_CTN",
                                              str_detect(treat_name, "PermaNet P191") ~ "PermaNetP191",
                                              control ~ "untreated",
                                              TRUE ~ treat_name),
                                 net_product),
      net_gen = if_else(is.na(net_gen), 
                           case_when(str_detect(net_product, "InterceptorG2") ~ "LLING2", 
                                     str_detect(net_product, "Interceptor") ~ "LLIN",
                                     str_detect(net_product, "RoyalGuard") ~ "LLING2",
                                     str_detect(net_product, "PermaNet3") ~ "LLING2",
                                     str_detect(net_product, "PermaNet2") ~ "LLIN",
                                     str_detect(net_product, "PermaNetP191") ~ "LLING2",
                                     str_detect(net_product, "DawaPlus2") ~ "LLIN",
                                     str_detect(net_product, "OlysetPlus") ~ "LLING2",
                                     str_detect(net_product, "Olyset") ~ "LLIN",
                                     str_detect(net_product, "RoyalSentry") ~ "LLIN",
                                     str_detect(net_product, "chlorfenapyr_CTN") ~ "hand_treated",
                                     str_detect(net_product, "alphacy_CTN") ~ "hand_treated",
                                     control ~ "untreated"), 
                           net_gen),
      insecticide = if_else(is.na(insecticide), 
                                case_when(str_detect(net_product, "InterceptorG2") ~ "alphacypermethrin+chlorfenapyr",
                                          str_detect(net_product, "Interceptor") ~ "alphacypermethrin",
                                          str_detect(net_product, "RoyalGuard") ~ "alphacypermethrin+pyriproxifen",
                                          str_detect(net_product, "PermaNet3") ~ "Deltamethrin+PBO",
                                          str_detect(net_product, "PermaNet2") ~ "deltamethrin",
                                          str_detect(net_product, "DawaPlus2") ~ "deltamethrin",
                                          str_detect(net_product, "OlysetPlus") ~ "permethrin+PBO",
                                          str_detect(net_product, "Olyset") ~ "permethrin",
                                          str_detect(net_product, "RoyalSentry") ~ "alphacypermethrin",
                                          str_detect(net_product, "chlorfenapyr_CTN") ~ "chlorfenapyr",
                                          str_detect(net_product, "alphacy_CTN") ~ "alpha-cypermethrin",
                                          str_detect(net_product, "PermaNetP191") ~ "deltamethrin+chlorfenapyr", # according to Tom
                                          control ~ "none"),
                                insecticide),
      insecticide = str_replace(insecticide, "alpha(|-)cypermethrin", "alphacypermethrin"),
      insecticide = str_replace(insecticide, "lambda(|-)cyhalothrin", "lambdacyhalothrin"),
      N_washes = as.integer(N_washes),
      N_washes = if_else(is.na(N_washes), 
                         # check this!!! make more clever to detect any number of washes (probably needs regex)
                         # check whether assuming 0 as follows is OK
                         # needs a lot of care here: 20W needs to come before 0W!!!
                         case_when(str_detect(treat_name, "20W|20x") ~ 20L,
                                   str_detect(treat_name, "0W|0x") ~ 0L,
                                   str_detect(treat_name, "15W|15x") ~ 15L,
                                   #!str_detect(treat_name, "[0-9]++W") & control != TRUE ~ 0L, # this assumes unwashed if no wash indication in treat_name
                                   net_status == "unwashed" ~ 0L,
                                   TRUE ~ NA_integer_),
                         N_washes),
      net_status = if_else(is.na(net_status), 
                           case_when(control ~ "control",
                                     N_washes == 0 ~ "unwashed",
                                     N_washes > 0 ~ "washed",
                                     TRUE ~ NA_character_), 
                           net_status),
      net_status = str_to_lower(net_status)
    )
    
    
    
    # align covariates
    H_df <- H_df %>%
    mutate(hut_type = case_when(str_detect(hut_type, regex("west", ignore_case = TRUE)) ~ "West",
                                str_detect(hut_type, regex("east", ignore_case = TRUE)) ~ "East",
                                str_detect(hut_type, regex("ifakara", ignore_case = TRUE)) ~ "Ifakara",
                                TRUE ~ hut_type),
           species = case_when(str_detect(species, regex("gamb(|.*)s(|.)s(|.)", ignore_case = TRUE)) & !str_detect(species, regex("NON", ignore_case = TRUE)) ~ "gambiae ss",
                               str_detect(species, regex("gamb", ignore_case = TRUE)) & !str_detect(species, regex("NON", ignore_case = TRUE)) ~ "gambiae sl",
                               str_detect(species, regex("arab", ignore_case = TRUE)) & !str_detect(species, regex("NON", ignore_case = TRUE)) ~ "arabiensis",
                               str_detect(species, regex("fune", ignore_case = TRUE)) & !str_detect(species, regex("NON", ignore_case = TRUE)) ~ "funestus",
                               str_detect(species, regex("nili", ignore_case = TRUE)) & !str_detect(species, regex("NON", ignore_case = TRUE)) ~ "nili",
                               str_detect(species, regex("An_NON_gamb_sl", ignore_case = TRUE))~ "NON(gamb_sl)",
                               str_detect(species, regex("An_NON_gambsl_fune_nili", ignore_case = TRUE))~ "NON(gamb_sl, fune, nili)",
                               TRUE ~ species)
    )

  
  
  
  
  
  # align country names
  # currently all fine!
   H_df <- H_df %>%
     mutate(country = case_when(country == "BurkinaFaso" ~ "Burkina Faso",
                                country == "CotedIvoire" ~ "Cote d'Ivoire",
                                TRUE ~ country) )
  
  # check for remaining ambiguities in country names
  # gives list of similar country names according to specified similarity score 
  # by the similarity metric 'Optimal string aligment, (restricted Damerau-Levenshtein distance)'
  all_country_names <- unique(H_df$country) 
  similarity_index <- 0.8
  
  all_country_names_simmatrix <- stringdist::stringsimmatrix(all_country_names, all_country_names) # compute pairwise similarity score
  all_country_names_simmatrix[lower.tri(all_country_names_simmatrix, diag = TRUE)] <- 0 # set lower triangular matrix including diag to 0 
  all_country_names_simmatrix <- all_country_names_simmatrix > similarity_index
  all_country_names_sim <- which(all_country_names_simmatrix, arr.ind = TRUE)
  all_country_names_simpair <- tibble( cbind(all_country_names[all_country_names_sim[ ,1]], all_country_names[all_country_names_sim[ ,2]]))
  
  if (nrow(all_country_names_simpair) != 0){
    writeLines(paste0("WARNING: There are remaining ambiguities in the country unit names withing some countries."))
    print.data.frame(all_country_names_simpair)
    writeLines("\n")
  }
  
  
  
  # align site names
  # currently all fine!
  # H_df <- H_df %>%
  #   mutate(site = case_when(site == "" ~ "") )
  
  # check for remaining ambiguities in site names
  # gives list of similar site names according to specified similarity score 
  # by the similarity metric 'Optimal string alignment, (restricted Damerau-Levenshtein distance)'
  all_site_names <- H_df %>% 
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
    writeLines(paste0("WARNING: There are remaining ambiguities in the site unit names."))
    print.data.frame(all_site_names_simpair)
    writeLines("\n")
  }
  
  
  # deal with gravid ones 
  # assume that gravid ones are unfed and that they were not yet counted as unfed
  # NOTE: For Trial_code 208, outcome_timepoint 24h and 72h, gravid ones have already been counted as unfed (but no problem here as gravid_live and gravid_dead were set to NA there)
  H_df <- H_df %>%
    mutate(unf_live = if_else(!is.na(gravid_live), unf_live + gravid_live, unf_live),
           unf_dead = if_else(!is.na(gravid_dead), unf_dead + gravid_dead, unf_dead),
           # and correct tot_unf if gravid ones were given and if tot_unf was given
           tot_unf = if_else( !is.na(gravid_live) & !is.na(tot_unf), tot_unf + gravid_live, tot_unf),
           tot_unf = if_else( !is.na(gravid_dead) & !is.na(tot_unf), tot_unf + gravid_dead, tot_unf),
           ) %>%
    # and omit gravid counts to prevent confusion 
    select(!c(gravid_live, gravid_dead))
  
  outcomes_EHT <- setdiff(outcomes_EHT, c("gravid_live", "gravid_dead")) 
  
  
  # compute totals
  H_df <- H_df %>%
    mutate(tot_unf = if_else(is.na(tot_unf), unf_live + unf_dead, tot_unf),
           tot_blf = if_else(is.na(tot_blf), bf_live + bf_dead, tot_blf),
           tot_live = if_else(is.na(tot_live), unf_live + bf_live, tot_live),
           tot_dead = if_else(is.na(tot_dead), unf_dead + bf_dead, tot_dead),
           total = if_else(is.na(total), if_else( tot_live + tot_dead == tot_unf + tot_blf, tot_live + tot_dead, NA_real_), total),
           # compute tot_live as total - tot_dead for cases where split-off counts not given (delayed mortality without unfed/fed)
           tot_live = if_else(is.na(tot_live), total - tot_dead, tot_live),
    )
  
  
  
  # data checks and omit inconsistent data points
    H_df <- H_df %>% 
    mutate(data_inconsistent = total != unf_live + unf_dead + bf_live + bf_dead | 
                               total != tot_unf + tot_blf | 
                               total != tot_live + tot_dead |
                               total < tot_unf |
                               total < tot_blf |
                               total < tot_live |
                               total < tot_dead |
                               tot_live != unf_live + bf_live |
                               tot_dead != unf_dead + bf_dead |
                               tot_unf != unf_live + unf_dead |
                               tot_blf != bf_live + bf_dead |
                               if_any(any_of(outcomes_EHT), ~ .  < 0)
          )
  
  if (TRUE %in% H_df$data_inconsistent){
    print("The following EHT data contains inconsistent information and is removed.")
    print.data.frame(H_df %>% filter(data_inconsistent == TRUE))
  }
  
  H_df <- H_df %>% 
    filter(!(data_inconsistent) | is.na(data_inconsistent)) %>%
    select(-data_inconsistent)

  
  # check for data points missing essential information and omit those

  
  
  ## omit data points with missing key data
    H_df <- H_df %>% 
      mutate(data_missing = case_when(is.na(total) ~ TRUE,
                                    is.na(tot_dead) ~ TRUE,
                                    TRUE ~ FALSE))
  
  if (TRUE %in% H_df$data_missing){
    print("The following EHT data points missed key data and are removed.")
    print.data.frame(H_df %>% filter(data_missing == TRUE))
  }
  
  H_df <- H_df %>%
    filter(!data_missing) |>
    select(-data_missing)
  
  
  # define new variable 'treat' from 'control'
  H_df <- H_df %>%
    mutate(treat = as.integer(!as.logical(control)))
  
  
  return <- list(EHT = H_df, covars_EHT = covars_EHT, outcomes_EHT = outcomes_EHT, vars_EHT = vars_EHT)
}

