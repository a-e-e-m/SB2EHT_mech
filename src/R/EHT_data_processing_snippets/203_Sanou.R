H_df_add<- readxl::read_xlsx(file.path(select(filter(list_data_sets, Trial_code==Trial_code_i), path2EHTfile), 
                                       select(filter(list_data_sets, Trial_code==Trial_code_i), EHT_file)),
                             sheet = paste0(select(filter(list_data_sets, Trial_code==Trial_code_i), EHT_sheet)))
# needed since strange character after Village!
names(H_df_add)[str_detect(names(H_df_add), "Village")] <- "Village"

H_df_add <- H_df_add %>%
  rename(treat_name = "Treatment-C",
         site = Village,
         sleeper = Sleeper,
         hut = Case,
         day = Day,
         round = Round,
         week = Week,
         date = Date) %>%
  mutate(Trial_code = Trial_code_i, 
         control = str_detect(treat_name, "UTN"), 
         # check all of this with authors!!! put here instead of further down since maybe too specific
         net_product = case_when(str_detect(treat_name, "IG2") ~ "InterceptorG2", 
                                 str_detect(treat_name, "IG1") ~ "Interceptor",
                                 str_detect(treat_name, "RG") ~ "RoyalGuard",
                                 str_detect(treat_name, "P3") ~ "PermaNet3",
                                 str_detect(treat_name, "P2") ~ "PermaNet2",
                                 str_detect(treat_name, "OSP") ~ "OlysetPlus",
                                 str_detect(treat_name, "UTN") ~ NA_character_),
         N_washes = case_when(str_detect(treat_name, "NU") ~ 0L,
                              TRUE ~ NA_integer_),
         net_status = case_when(str_detect(treat_name, "NU") ~ "unwashed",
                                str_detect(treat_name, "NW") ~ "washed",
                                str_detect(treat_name, regex("age", ignore_case = TRUE)) ~ "aged", # pool all age categories
                                TRUE ~ NA_character_),
         hut_type = "West", # according to email with Antoine
         country = "BurkinaFaso",
         year = year(date),
         week = as.integer(str_remove(week, "W")),
         day = as.integer(str_remove(day, "d")),
         hut = as.integer(str_remove(hut, "C")),
         sleeper = as.integer(str_remove(sleeper, "S")))

# summarise over different areas (room, net, veranda) and 'collect' wanted outcomes (outcome_cat) by reshaping data frame    
H_df_add_cat <- H_df_add %>%
  mutate(outcome_cat = case_when(str_detect(Physiology, "Unfed") & str_detect(Status, "Alive") ~ "unf_live",
                                 str_detect(Physiology, "Unfed") & str_detect(Status, "Dead") ~ "unf_dead",
                                 str_detect(Physiology, "Fed|Partially.fed") & str_detect(Status, "Alive") ~ "bf_live",
                                 str_detect(Physiology, "Fed|Partially.fed") & str_detect(Status, "Dead") ~ "bf_dead",
                                 str_detect(Physiology, "Gravid") & str_detect(Status, "Alive") ~ "gravid_live",
                                 str_detect(Physiology, "Gravid") & str_detect(Status, "Dead") ~ "gravid_dead")) %>%
  select(!c(Status, Physiology)) %>%
  filter(!is.na(outcome_cat))    


# outcome_timepoint == morning
# gambiae sl  
H_df_add_sum_gambsl <- H_df_add_cat %>%
  group_by(across(c(intersect(covars_EHT, names(H_df_add_cat)), "outcome_cat"))) %>%
  dplyr::summarise(An_gamb_sl = sum(`An. gambiae sl`)) %>%
  ungroup() %>%
  pivot_wider(names_from = outcome_cat, values_from = An_gamb_sl) %>%
  mutate(outcome_timepoint = "morning",
         species = "An_gamb_sl"
  ) %>%
  select(any_of(vars_EHT))

H_df <- H_df %>% bind_rows(H_df_add_sum_gambsl)
rm(H_df_add_sum_gambsl)

# non gambiae sl  
H_df_add_sum_NONgambsl <- H_df_add_cat %>%
  group_by(across(c(intersect(covars_EHT, names(H_df_add_cat)), "outcome_cat"))) %>%
  dplyr::summarise(An_NON_gamb_sl = sum(Other.Anophelines)) %>%
  ungroup() %>%
  pivot_wider(names_from = outcome_cat, values_from = An_NON_gamb_sl) %>%
  mutate(outcome_timepoint = "morning",
         species = "An_NON_gamb_sl"
  ) %>%
  select(any_of(vars_EHT))

H_df <- H_df %>% bind_rows(H_df_add_sum_NONgambsl)
rm(H_df_add_sum_NONgambsl, H_df_add_cat)

# additional mortality for outcome_timepoint == 24h
# WARNING: Only An_gamb_sl
H_df_add_aliveminus_24h <- H_df_add %>%
  mutate(`An. gambiae sl` = if_else(Status == "Alive" & Tot.gamb.follow != 0, # for safety, `24h_dead` should otherwise be 0 anyway, second clause needed to avoid NaN
                                    round(`An. gambiae sl` - `24h_dead` * `An. gambiae sl` / Tot.gamb.follow),
                                    `An. gambiae sl`)
  )

H_df_add_deadplus_24h <- H_df_add %>%
  filter(Status == "Alive" & Tot.gamb.follow != 0) %>% # for safety, should otherwise be 0 anyway
  mutate(`An. gambiae sl` = round(`24h_dead` * `An. gambiae sl` / Tot.gamb.follow),
         Status = "Dead")


H_df_add_sum_gambsl_24h <- rbind(H_df_add_aliveminus_24h, H_df_add_deadplus_24h)  %>%
  mutate(outcome_cat = case_when(str_detect(Physiology, "Unfed") & str_detect(Status, "Alive") ~ "unf_live",
                                 str_detect(Physiology, "Unfed") & str_detect(Status, "Dead") ~ "unf_dead",
                                 str_detect(Physiology, "Fed|Partially.fed") & str_detect(Status, "Alive") ~ "bf_live",
                                 str_detect(Physiology, "Fed|Partially.fed") & str_detect(Status, "Dead") ~ "bf_dead",
                                 str_detect(Physiology, "Gravid") & str_detect(Status, "Alive") ~ "gravid_live",
                                 str_detect(Physiology, "Gravid") & str_detect(Status, "Dead") ~ "gravid_dead")) %>%
  select(!c(Status, Physiology)) %>%
  filter(!is.na(outcome_cat))  %>%
  group_by(across(c(intersect(covars_EHT, names(H_df_add)), "outcome_cat"))) %>%
  dplyr::summarise(An_gamb_sl = sum(`An. gambiae sl`)) %>%
  ungroup() %>%
  pivot_wider(names_from = outcome_cat, values_from = An_gamb_sl) %>%
  mutate(outcome_timepoint = "24h",
         species = "An_gamb_sl",
  ) %>%
  select(any_of(vars_EHT))


H_df <- H_df %>% bind_rows(H_df_add_sum_gambsl_24h)
rm(H_df_add_sum_gambsl_24h, H_df_add_aliveminus_24h, H_df_add_deadplus_24h)


# additional mortality for outcome_timepoint == 48h
# WARNING: Only An_gamb_sl
H_df_add_aliveminus_48h <- H_df_add %>%
  mutate(`An. gambiae sl` = if_else(Status == "Alive" & Tot.gamb.follow != 0, # for safety, `72h_dead` should otherwise be 0 anyway, second clause needed to avoid NaN
                                    round(`An. gambiae sl` - (`24h_dead` + `48H_dead`) * `An. gambiae sl` / Tot.gamb.follow),
                                    `An. gambiae sl`)
  )

H_df_add_deadplus_48h <- H_df_add %>%
  filter(Status == "Alive" & Tot.gamb.follow != 0) %>% # for safety, should otherwise be 0 anyway
  mutate(`An. gambiae sl` = round((`24h_dead` + `48H_dead`) * `An. gambiae sl` / Tot.gamb.follow),
         Status = "Dead")


H_df_add_sum_gambsl_48h <- rbind(H_df_add_aliveminus_48h, H_df_add_deadplus_48h)  %>%
  mutate(outcome_cat = case_when(str_detect(Physiology, "Unfed") & str_detect(Status, "Alive") ~ "unf_live",
                                 str_detect(Physiology, "Unfed") & str_detect(Status, "Dead") ~ "unf_dead",
                                 str_detect(Physiology, "Fed|Partially.fed") & str_detect(Status, "Alive") ~ "bf_live",
                                 str_detect(Physiology, "Fed|Partially.fed") & str_detect(Status, "Dead") ~ "bf_dead",
                                 str_detect(Physiology, "Gravid") & str_detect(Status, "Alive") ~ "gravid_live",
                                 str_detect(Physiology, "Gravid") & str_detect(Status, "Dead") ~ "gravid_dead")) %>%
  select(!c(Status, Physiology)) %>%
  filter(!is.na(outcome_cat))  %>%
  group_by(across(c(intersect(covars_EHT, names(H_df_add)), "outcome_cat"))) %>%
  dplyr::summarise(An_gamb_sl = sum(`An. gambiae sl`)) %>%
  ungroup() %>%
  pivot_wider(names_from = outcome_cat, values_from = An_gamb_sl) %>%
  mutate(outcome_timepoint = "48h",
         species = "An_gamb_sl",
  ) %>%
  select(any_of(vars_EHT))

H_df <- H_df %>% bind_rows(H_df_add_sum_gambsl_48h)
rm(H_df_add_sum_gambsl_48h, H_df_add_aliveminus_48h, H_df_add_deadplus_48h)


# additional mortality for outcome_timepoint == 72h
# WARNING: Only An_gamb_sl
H_df_add_aliveminus_72h <- H_df_add %>%
  mutate(`An. gambiae sl` = if_else(Status == "Alive" & Tot.gamb.follow != 0, # for safety, `72h_dead` should otherwise be 0 anyway, second clause needed to avoid NaN
                                    round(`An. gambiae sl` - (`24h_dead` + `48H_dead` + `72h_dead`) * `An. gambiae sl` / Tot.gamb.follow),
                                    `An. gambiae sl`)
  )

H_df_add_deadplus_72h <- H_df_add %>%
  filter(Status == "Alive"& Tot.gamb.follow != 0) %>% # for safety, should otherwise be 0 anyway
  mutate(`An. gambiae sl` = round((`24h_dead` + `48H_dead` + `72h_dead`) * `An. gambiae sl` / Tot.gamb.follow),
         Status = "Dead")


H_df_add_sum_gambsl_72h <- rbind(H_df_add_aliveminus_72h, H_df_add_deadplus_72h)  %>%
  mutate(outcome_cat = case_when(str_detect(Physiology, "Unfed") & str_detect(Status, "Alive") ~ "unf_live",
                                 str_detect(Physiology, "Unfed") & str_detect(Status, "Dead") ~ "unf_dead",
                                 str_detect(Physiology, "Fed|Partially.fed") & str_detect(Status, "Alive") ~ "bf_live",
                                 str_detect(Physiology, "Fed|Partially.fed") & str_detect(Status, "Dead") ~ "bf_dead",
                                 str_detect(Physiology, "Gravid") & str_detect(Status, "Alive") ~ "gravid_live",
                                 str_detect(Physiology, "Gravid") & str_detect(Status, "Dead") ~ "gravid_dead")) %>%
  select(!c(Status, Physiology)) %>%
  filter(!is.na(outcome_cat))  %>%
  group_by(across(c(intersect(covars_EHT, names(H_df_add)), "outcome_cat"))) %>%
  dplyr::summarise(An_gamb_sl = sum(`An. gambiae sl`)) %>%
  ungroup() %>%
  pivot_wider(names_from = outcome_cat, values_from = An_gamb_sl) %>%
  mutate(outcome_timepoint = "72h",
         species = "An_gamb_sl",
  ) %>%
  select(any_of(vars_EHT))

# check 72h with column 'Cum_72H' before adding
###############################################
H_df_add_aliveminus_72h_cum <- H_df_add %>%
  mutate(`An. gambiae sl` = if_else(Status == "Alive" & Tot.gamb.follow != 0, # for safety, `72h_dead` should otherwise be 0 anyway, second clause needed to avoid NaN
                                    round(`An. gambiae sl` - Cum_72H * `An. gambiae sl` / Tot.gamb.follow),
                                    `An. gambiae sl`)
  )

H_df_add_deadplus_72h_cum <- H_df_add %>%
  filter(Status == "Alive" & Tot.gamb.follow != 0) %>% 
  mutate(`An. gambiae sl` = round(Cum_72H * `An. gambiae sl` / Tot.gamb.follow),
         Status = "Dead")


H_df_add_sum_gambsl_72h_cum <- rbind(H_df_add_aliveminus_72h_cum, H_df_add_deadplus_72h_cum)  %>%
  mutate(outcome_cat = case_when(str_detect(Physiology, "Unfed") & str_detect(Status, "Alive") ~ "unf_live",
                                 str_detect(Physiology, "Unfed") & str_detect(Status, "Dead") ~ "unf_dead",
                                 str_detect(Physiology, "Fed|Partially.fed") & str_detect(Status, "Alive") ~ "bf_live",
                                 str_detect(Physiology, "Fed|Partially.fed") & str_detect(Status, "Dead") ~ "bf_dead",
                                 str_detect(Physiology, "Gravid") & str_detect(Status, "Alive") ~ "gravid_live",
                                 str_detect(Physiology, "Gravid") & str_detect(Status, "Dead") ~ "gravid_dead")) %>%
  select(!c(Status, Physiology)) %>%
  filter(!is.na(outcome_cat))  %>%
  group_by(across(c(intersect(covars_EHT, names(H_df_add)), "outcome_cat"))) %>%
  dplyr::summarise(An_gamb_sl = sum(`An. gambiae sl`)) %>%
  ungroup() %>%
  pivot_wider(names_from = outcome_cat, values_from = An_gamb_sl) %>%
  mutate(outcome_timepoint = "72h",
         species = "An_gamb_sl",
  ) %>%
  select(any_of(vars_EHT))

equal_check_72h_cum72h <- all_equal(H_df_add_sum_gambsl_72h, H_df_add_sum_gambsl_72h_cum)
print(paste0("For Sanou's data, a comparison between incremental 24h, 48h, 72h mortalities with cumulative mortality (column 'Cum_72H) is performed. Are they equal?: ", equal_check_72h_cum72h))
#################################

H_df <- H_df %>% bind_rows(H_df_add_sum_gambsl_72h)
rm(H_df_add, H_df_add_sum_gambsl_72h, H_df_add_aliveminus_72h, H_df_add_deadplus_72h, equal_check_72h_cum72h)    


remaining_data_sets <- setdiff(remaining_data_sets, H_df$Trial_code)