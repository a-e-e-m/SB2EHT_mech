# This adds manual matches for site-year

manualmatch_site_year <- function(BA_groups, EHT_groups, N_groups_previous){
  # input: non-automatically matched BA and EHT groups, last group_number from automatic match
  # output: manually matched BA and EHT groups, with group_numbers assigned (consecutive from last automatic match)
  #         remaining unmatched BA and EHT groups 
  
  EHT_groups <- EHT_groups %>% 
    mutate(group_sign = case_when(Trial_code == 206 & site == "Ganhoua" & year == 2021 & insecticide == "alphacypermethrin" ~ "j"
                                  )
    )
  
  
  BA_groups <- BA_groups %>% 
    mutate(group_sign = case_when(Trial_code == 206 & site == "Zakpota" & year == 2021 & insecticide == "alphacypermethrin" ~ "j"
                                  )
    )
  
  # get consistent group_numbers, consecutive after last present, from group_sign
  group_sign_arr <- setdiff(lubridate::intersect( unique(EHT_groups$group_sign), unique(BA_groups$group_sign)), NA)
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