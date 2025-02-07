# function to load and process model predictions for publication plot of dose response curves for mech model paper

# output : 
# prob_D_b
# prob_D_b_summary
# list_dfs_doseresponse
# list_dfs_doseresponse_summary



modelpred_4_pubplot <- function(i){

# config and prepare
###############################
###############################

#libraries
library("HDInterval")
library("tidyr")
library("dplyr")
library("reshape2")
library("rstan")
library("stringr")
library("data.table")

# source script containing needed functions
functionsfolder <- file.path('./src/R/functions')
files.source <- list.files(functionsfolder)
invisible(sapply(files.source, function(x) source(paste0(functionsfolder, "/", x))))


models2run <- read.csv(file = './models2run.csv', sep = ';', stringsAsFactors = FALSE)


# loading specifications
#############
with_feeding <- as.logical(models2run$with_feeding[i])
BA_only <- as.logical(models2run$BA_only[i])
dir_i <- paste0("fitting/", models2run$run[i]) # folder
trials <- readRDS(file.path(dir_i, 'trials.rds'))


# loading real data
#############
B_all <- readRDS(file.path(dir_i, 'B_all.rds'))


# loading model fit
#############
fit_i <- readRDS(file.path(dir_i, 'fit.rds'))

# take posterior sub-sample of size param_ss if MCMC sample is larger
param_s <- length(rstan::extract(fit_i, pars = "sigma_v", permuted = TRUE)[["sigma_v"]])
param_ss <- min(2000, param_s)
index_sample <- sort(sample(x = 1:param_ss, size = param_ss, replace = FALSE))

# load fitted model parameters, subsample if necessary
mu_d <- rstan::extract(fit_i, pars = "mu_d", permuted = TRUE)[["mu_d"]][index_sample,]
sigma_d <- rstan::extract(fit_i, pars = "sigma_d", permuted = TRUE)[["sigma_d"]][index_sample,]
sigma_v <- rstan::extract(fit_i, pars = "sigma_v", permuted = TRUE)[["sigma_v"]][index_sample]
# switch needed because p_b_control sometimes group specific and sometimes not
if (length(dim(rstan::extract(fit_i, pars = "p_b_control", permuted = TRUE)[["p_b_control"]])) == 1){
  p_b_control <- rstan::extract(fit_i, pars = "p_b_control", permuted = TRUE)[["p_b_control"]][index_sample]
} else {
  p_b_control <- rstan::extract(fit_i, pars = "p_b_control", permuted = TRUE)[["p_b_control"]][index_sample,]
}



# load model predictions (probability) per data point
# inv_logit transform if necessary
inv_logit <- function(x){
  1 / ( 1 + exp( -x ))
}

if ("prob_D_b" %in% fit_i@model_pars){
  prob_D_b <-  rstan::extract(fit_i, pars = "prob_D_b")[["prob_D_b"]][index_sample,]
} else if ("problogit_D_b" %in% fit_i@model_pars){
  prob_D_b <- inv_logit( rstan::extract(fit_i, pars = "problogit_D_b")[["problogit_D_b"]][index_sample,])
} else {
  prob_D_b <- numeric(0)
}


# data processing 
#################

# summaries data and add estimates, and plot for BA 
# connect model estimates of death probabilty and real data, summarise
# NOTE: this used the prob estimates from the stan output directly, not using the other parameters
col_D_b <- paste("V", seq(1, nrow(B_all),1), sep = "")
prob_D_b <- as.data.table(prob_D_b)[, iter := .I]
prob_D_b <- melt(prob_D_b, measure = col_D_b, variable.name = "data_point", value.name = c("prob"))
prob_D_b <- prob_D_b[, `:=`("data_point" = as.integer(str_remove(data_point, "V")))]

B_all <- as.data.table(B_all)[, data_point := .I]
setkey(B_all, ON = data_point)
setkey(prob_D_b, ON = data_point)
prob_D_b <- B_all[prob_D_b, nomatch = 0] 

prob_D_b <- prob_D_b |>
  left_join(group_number_conversion) |>
  filter(group_number_conversion < 8) |>
  select(!group_number) |>
  rename(group_number = group_number_conversion)

prob_D_b_summary <- prob_D_b[
  ,
  list(
    "pred_median" = median(prob),
    "pred_hdi_l" = hdi(prob)[["lower"]],
    "pred_hdi_r" = hdi(prob)[["upper"]]
  ),
  .(group_number, control, treat, times_disc_dose)
]


B_all_summary <- B_all |>
  group_by(group_number, control, treat, times_disc_dose, int_dose_available) |>
  # be warned: this is not fully general but probably covers all group_by variants wanted
  summarise(D_b = sum(D_b),
            N_b = sum(N_b),
            Trial_code = paste(unique(Trial_code), collapse = ', '),
            country = paste(unique(country), collapse = ', '),
            site = paste(unique(site), collapse = ', '),
            year = paste(unique(year), collapse = ', '),
            insecticide = paste(unique(insecticide), collapse = ', '),
            test_type = paste(unique(test_type), collapse = ', '),
            .groups = "keep"
  ) |>
  mutate(insecticide = case_when( str_detect(insecticide, ",") ~ "pyrethroid",
                                  TRUE ~ insecticide),
         MLE = D_b / N_b, # note this is equivalent to the mode of beta(D_b + 1, N_b - D_b +1)
         q025 = qbeta(0.025, shape1 = D_b + 1, shape2 = N_b - D_b + 1),
         q975 = qbeta(0.975, shape1 = D_b + 1, shape2 = N_b - D_b + 1),
         median = qbeta(0.5, shape1 = D_b + 1, shape2 = N_b - D_b + 1),
         hdi_l = hdi(qbeta, 0.95, shape1 = D_b + 1, shape2 = N_b - D_b + 1)[[1]],
         hdi_r = hdi(qbeta, 0.95, shape1 = D_b + 1, shape2 = N_b - D_b + 1)[[2]]
  ) |> 
  left_join(prob_D_b_summary, by = c("group_number", "control", "treat", "times_disc_dose")) 


# dose response curves per group with facets (only for groups with intensity dose data)
###################
# compute predictions for dose-response
times_disc_dose <- seq(0, 15, 0.05)
list_dfs_doseresponse <- vector(mode = "list", length = trials[["nT"]])
list_dfs_doseresponse_summary <- vector(mode = "list", length = trials[["nT"]])

# lookup for group_number for int_dose groups
# contains quick hack to exclude data from Tanzania
lookup_group_number <- B_all_summary |>
  filter(int_dose_available == 1, country == "BurkinaFaso") |>
  ungroup() |>
  filter(insecticide != "none") |>
  select(group_number, country, site, year, insecticide) |>
  unique()

# loop over int_dose groups only
for (ii in seq_along(sort(unique(lookup_group_number$group_number)))) {
  #switch needed because p_b_control sometimes group specific and sometimes not
  if (length(dim(p_b_control)) == 1){
    list_dfs_doseresponse[[ii]] <- as.data.table(expand_grid(tibble(mu_d = mu_d[,lookup_group_number$group_number[ii]], sigma_d = sigma_d[,lookup_group_number$group_number[ii]], p_b_control = p_b_control, sigma_v, iter = 1:length(sigma_v)), times_disc_dose))
  }else{
    list_dfs_doseresponse[[ii]] <- as.data.table(expand_grid(tibble(mu_d = mu_d[,lookup_group_number$group_number[ii]], sigma_d = sigma_d[,lookup_group_number$group_number[ii]], p_b_control = p_b_control[,lookup_group_number$group_number[ii]], sigma_v, iter = 1:length(sigma_v)), times_disc_dose))
  }
  list_dfs_doseresponse[[ii]] <- list_dfs_doseresponse[[ii]][, p := p_b_control + ( 1 - p_b_control) * pnorm( (log(times_disc_dose) - mu_d ) / sqrt(sigma_v^2 + sigma_d^2 ))]
  list_dfs_doseresponse[[ii]] <- list_dfs_doseresponse[[ii]][, group_number := lookup_group_number$group_number[ii]]
  list_dfs_doseresponse[[ii]] <- tibble(list_dfs_doseresponse[[ii]]) |>
    left_join(lookup_group_number, by = join_by(group_number))
  list_dfs_doseresponse_summary[[ii]] <- tibble(list_dfs_doseresponse[[ii]]) |> 
    group_by(group_number, times_disc_dose) |>
    summarise(p_median = median(p), .groups = "keep") |>
    left_join(lookup_group_number, by = join_by(group_number))
}

return(list("prob_D_b" = prob_D_b, 
            "prob_D_b_summary" = prob_D_b_summary, 
            "list_dfs_doseresponse" = list_dfs_doseresponse, 
            "list_dfs_doseresponse_summary" = list_dfs_doseresponse_summary)
       )
}