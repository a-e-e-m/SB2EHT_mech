# This matches EHT trials with only one arm to the corresponding BA data, irrespective of the insecticide within the pyrethroid class.
# I.e. this pools all pyrethroids on BA side if corresponding EHT has only one arm, irrespective of the insecticide used in that net (given it is a single pyrethroid).



manualmatch_insecticides_loose <- function(BA_groups, EHT_groups, N_groups_previous){
# input: non-automatically matched BA and EHT groups, last group_number from automatic match
# output: manually matched BA and EHT groups, with group_numbers assigned (consecutive from last automatic match)
#         remaining unmatched BA and EHT groups 
  
  EHT_groups <- EHT_groups %>% 
    mutate(group_sign = case_when(Trial_code == 13 & site == "Cove" & year == 2012 & insecticide == "alphacypermethrin" ~ "a",
                                  Trial_code == 15 & site == "Cove" & year == 2015 & insecticide == "alphacypermethrin" ~ "b",
                                  Trial_code == 36 & site == "LowerMoshi" & year == 2014 & insecticide == "alphacypermethrin" ~ "c",
                                  Trial_code == 35 & site == "Muheza" & year == 2015 & insecticide == "alphacypermethrin" ~ "d",
                                  Trial_code == 25 & site == "Lupiro" & year == 2015 & insecticide == "alphacypermethrin" ~ "e",
                                  Trial_code == 26 & site == "Muheza" & year == 2015 & insecticide == "alphacypermethrin" ~ "f",
                                  Trial_code == 27 & site == "Yaokoffikro" & year == 2000 & insecticide == "deltamethrin" ~ "g",
                                  Trial_code == 33 & site == "VDK" & year == 2014 & insecticide == "alphacypermethrin" ~ "h",
                                  Trial_code == 34 & site == "Cove" & year == 2015 & insecticide == "alphacypermethrin" ~ "i",
                                  Trial_code == 30 & site == "WaliKunda" & year == 1986 & insecticide == "cypermethrin" ~ "j",
                                  Trial_code == 30 & site == "WaliKunda" & year == 1986 & insecticide == "deltamethrin" ~ "j",
                                  Trial_code == 30 & site == "WaliKunda" & year == 1986 & insecticide == "lambdacyhalothrin" ~ "j",
                                  Trial_code == 30 & site == "WaliKunda" & year == 1986 & insecticide == "permethrin" ~ "j"
                                  )
           )
  
  
  BA_groups <- BA_groups %>% 
    mutate(group_sign = case_when(Trial_code == 13 & site == "Cove" & year == 2012 & insecticide == "permethrin" ~ "a",
                                  Trial_code == 13 & site == "Cove" & year == 2012 & insecticide == "deltamethrin" ~ "a",
                                  Trial_code == 15 & site == "Cove" & year == 2015 & insecticide == "permethrin" ~ "b",
                                  Trial_code == 15 & site == "Cove" & year == 2015 & insecticide == "deltamethrin" ~ "b",
                                  Trial_code == 36 & site == "LowerMoshi" & year == 2014 & insecticide == "permethrin" ~ "c",
                                  Trial_code == 35 & site == "Muheza" & year == 2015 & insecticide == "permethrin" ~ "d",
                                  Trial_code == 25 & site == "Lupiro" & year == 2015 & insecticide == "permethrin" ~ "e",
                                  Trial_code == 25 & site == "Lupiro" & year == 2015 & insecticide == "deltamethrin" ~ "e",
                                  Trial_code == 25 & site == "Lupiro" & year == 2015 & insecticide == "lambdacyhalothrin" ~ "e",
                                  Trial_code == 26 & site == "Muheza" & year == 2015 & insecticide == "permethrin" ~ "f",
                                  Trial_code == 27 & site == "Yaokoffikro" & year == 2000 & insecticide == "pyrethroid" ~ "g",
                                  Trial_code == 33 & site == "VDK" & year == 2014 & insecticide == "permethrin" ~ "h",
                                  Trial_code == 33 & site == "VDK" & year == 2014 & insecticide == "deltamethrin" ~ "h",
                                  Trial_code == 34 & site == "Cove" & year == 2015 & insecticide == "permethrin" ~ "i",
                                  Trial_code == 34 & site == "Cove" & year == 2015 & insecticide == "deltamethrin" ~ "i",
                                  Trial_code == 30 & site == "WaliKunda" & year == 1986 & insecticide == "pyrethroid" ~ "j"
                                  )
    )

  # get consistent group_numbers, consecutive after last present, from group_sign
  group_sign_arr <- setdiff( lubridate::intersect( unique(EHT_groups$group_sign), unique(BA_groups$group_sign)), NA)
  group_number_arr <- seq_along(group_sign_arr) + N_groups_previous
  
  EHT_groups$group_number <- group_number_arr[ match(EHT_groups$group_sign, group_sign_arr) ]
  EHT_groups <- EHT_groups %>% select(!group_sign)
  BA_groups$group_number <- group_number_arr[ match(BA_groups$group_sign, group_sign_arr) ]
  BA_groups <- BA_groups %>% select(!group_sign)
  
  # get dfs for matched and still unmatched groups
  EHT_groups_matched <- EHT_groups %>%
    filter(!is.na(group_number))
  
  EHT_groups_unmatched <- EHT_groups %>%
    filter(is.na(group_number))  
  
  BA_groups_matched <- BA_groups %>%
    filter(!is.na(group_number))
  
  BA_groups_unmatched <- BA_groups %>%
    filter(is.na(group_number))
  
  
  return(list(BA_groups_newmatch = BA_groups_matched, EHT_groups_newmatch = EHT_groups_matched, BA_groups_unmatched = BA_groups_unmatched, EHT_groups_unmatched = EHT_groups_unmatched))
}