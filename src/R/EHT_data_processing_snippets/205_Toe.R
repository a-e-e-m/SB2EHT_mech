H_df_add<- readxl::read_xlsx(file.path(select(filter(list_data_sets, Trial_code==Trial_code_i), path2EHTfile), 
                                       select(filter(list_data_sets, Trial_code==Trial_code_i), EHT_file)),
                             sheet = paste0(select(filter(list_data_sets, Trial_code==Trial_code_i), EHT_sheet))) %>%
  rename(treat_name = "Treatment",
         site = Locality,
         sleeper = Sleeper,
         hut = Hut,
         day = Day,
         date = Date) %>%
  mutate(Trial_code = Trial_code_i, 
         control = str_detect(treat_name, "Control"), 
         # check nets with author!!! put here instead of further down since maybe too specific
         net_product = case_when(str_detect(treat_name, "PermaNet 3.0") ~ "PermaNet3",
                                 str_detect(treat_name, "PermaNet 2.0") ~ "PermaNet2",
                                 str_detect(treat_name, "OlysetPlus") ~ "OlysetPlus",
                                 str_detect(treat_name, "Olyset") ~ "Olyset",
                                 str_detect(treat_name, "DawaPlus") ~ "DawaPlus2", # I checked with publication that this is actually DawaPlus2.0
                                 control ~ "untreated"), 
         net_status = "unwashed", # from Publication
         N_washes = 0L, # from Publication
         hut_type = "west", # from Publication
         holes = "WHO_6",
         country = "BurkinaFaso",
         outcome_timepoint = "24h", # guessed as no InterceptorG2, to check!!!
         year = year(date),
         day = as.integer(day),
         hut = as.integer(hut),
         sleeper = as.integer(str_remove(sleeper, "S")),
         species = "gambiae sl", # according to Rebecca's table (Trial_code 22 and 23)
         unf_live = case_when(str_detect(statut, "Unfed") ~ Alive,
                              TRUE ~ 0),
         bf_live = case_when(str_detect(statut, "Blood fed") ~ Alive,
                             TRUE ~ 0),
         unf_dead = case_when(str_detect(statut, "Unfed") ~ `Total dead`,
                              TRUE ~ 0),
         bf_dead = case_when(str_detect(statut, "Blood fed") ~ `Total dead`,
                             TRUE ~ 0),             
         tot_unf = case_when(str_detect(statut, "Unfed") ~ `Total collected`,
                             TRUE ~ 0),
         tot_blf = case_when(str_detect(statut, "Blood fed") ~ `Total collected`,
                             TRUE ~ 0),
         gravid_live = case_when(str_detect(statut, "Gravid") ~ Alive,
                                 TRUE ~ 0),
         gravid_dead = case_when(str_detect(statut, "Gravid") ~ `Total dead`,
                                 TRUE ~ 0),
         tot_live = Alive,
         tot_dead = `Total dead`,
         total = `Total collected`
  )

H_df_add <- H_df_add %>%
  group_by(across(intersect(covars_EHT, names(H_df_add)))) %>%
  summarise( unf_live = sum(unf_live),
             bf_live = sum(bf_live),
             unf_dead = sum(unf_dead),
             bf_dead = sum(bf_dead),
             gravid_live = sum(gravid_live),
             gravid_dead = sum(gravid_dead),
             tot_unf = sum(tot_unf),
             tot_blf = sum(tot_blf),
             tot_live = sum(tot_live),
             tot_dead = sum(tot_dead),
             total = sum(total)
  ) %>%
  ungroup() %>%
  select(any_of(vars_EHT))


H_df <- H_df %>% bind_rows(H_df_add)
rm(H_df_add)
remaining_data_sets <- setdiff(remaining_data_sets, H_df$Trial_code)