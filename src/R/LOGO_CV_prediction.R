# script to predict EHT mortality from SB data with LOGO fit
# produces and saves predictions and elpd estimates in fitting/dir_i
# produces and saves actual VS predicted plots in fitting/dir_i
# script takes integer as input, referring to a row in models2run, can also be set inside script


## config
################################
library(tidyverse)
library(tidybayes)
library(ggplot2)
library(extraDistr)

# source script containing needed functions
functionsfolder <- file.path('./src/R/functions')
files.source <- list.files(functionsfolder)
invisible(sapply(files.source, function(x) source(paste0(functionsfolder, "/", x))))


# user defined functions
################################
logit <- function(p) {
  log(p / (1 - p))
} 

inv_logit <- function(x){
  1 / ( 1 + exp( -x ))
}

safe_log_dbinom <- function(x, size, prob) {
  lchoose(size, x) +
    x * log(prob-1) +
    (size - x) * log1p(-prob)  # log1p(z) = log(1 + z), better near 0
}

log_mean_exp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m))) - log(length(x))
}

richardscurve_C1_A0_K1 <- function( x, B, M, nu){
  (1 + exp(- B * (x - M) ) ) ^ (-1 / nu) 
}

# loading specifications
################################
# read file with all model runs
models2run <- read.csv(file = './models2run.csv', sep = ';', stringsAsFactors = FALSE)

# select one line from models2run.csv
# or get index from argument passed to this script
args = commandArgs(trailingOnly=TRUE)
i <- as.numeric(args[1])

with_feeding <- as.logical(models2run$with_feeding[i])
BA_only <- as.logical(models2run$BA_only[i])
dir_i <- paste0("fitting/", models2run$run[i]) # folder
trials <- readRDS(file.path(dir_i, 'trials.rds'))
T_bint_indices <- trials[["T_bint_indices"]]
nT_int <- length(T_bint_indices)


# loading real data
################################
B_all <- readRDS(file.path(dir_i, 'B_all.rds')) |>
  mutate(country = case_when(country == "BurkinaFaso" ~ "Burkina Faso",
                             TRUE ~ country))

if (!BA_only){
  H_all <- readRDS(file.path(dir_i, 'H_all.rds')) |>
    mutate(country = case_when(country == "BurkinaFaso" ~ "Burkina Faso",
                               TRUE ~ country),
           H_row_id = row_number())
  
}else{H_all <- tibble()}

if (with_feeding){
  H_f_all <- readRDS(file.path(dir_i, 'H_f_all.rds')) |>
    mutate(country = case_when(country == "BurkinaFaso" ~ "Burkina Faso",
                               TRUE ~ country),
           H_row_id = row_number())
  
}else{H_f_all <- tibble()}


# loading LOGO model fits
################################
# get models2run_LOGO_run
file_models2run_LOGO_run <- file.path(".", "models2run4LOGO", paste0("models2run_logo_", models2run$run[i], ".csv"))
models2run_LOGO_run <- read.csv2(file_models2run_LOGO_run)

# collect samples for LOGO
predictions <- data.frame()
elpd <- data.frame()
for (k in seq(1, nT_int)){
  dir_k <- file.path("fitting", models2run_LOGO_run$run[k]) 
  fit_k <- readRDS(file.path(dir_k, 'fit.rds'))
  
  source(file.path("./src/R/model_pred_snipets", paste0(models2run$model[i], ".R")), local=TRUE)
  
  predictions <- bind_rows(predictions, draws_k_spread)
  elpd <- bind_rows(elpd, elpd_k)
}

# for models with feeding data
if (with_feeding){
  # check whether there is a mismatch between the group_numbers with feeding prob predictions and group_numbers declared to have feeding data in trials
  if (!all.equal(sort(unique(pull(filter(predictions, !is.na(prob_AF)), group_number))), trials$T_f_indices[trials$T_f_indices != 0 ])){
    message("WARNING: Mismatch between group_numbers with feeding prob prediction summaries and group_numbers declared to have feeding data in trials. Possibly some feeding probability predictions produced NAs.")
  }
}

# save predictions to run folder
saveRDS(predictions, file = file.path(dir_i, "predictions.rds"))

# save elpd to run folder
saveRDS(elpd, file = file.path(dir_i, "elpd.rds"))


# summarise
  # model predictions
  predictions_summary <- predictions |>
    group_by(across( any_of( c("group_number", "secondai_ID")) ) ) |>
    summarise(
      pred_median = median(prob_D_h),
      pred_q025 = quantile(prob_D_h, probs = 0.025), 
      pred_q975 = quantile(prob_D_h, probs = 0.975)
    )
  
  # also summarise feeding prob if applicable
  if (with_feeding) {
    predictions_summary <- predictions |>
      group_by(across( any_of( c("group_number", "secondai_ID")) ) ) |>
      summarise(
        pred_f_median = median(prob_AF, na.rm = TRUE),
        pred_f_q025 = quantile(prob_AF, probs = 0.025, na.rm = TRUE), 
        pred_f_q975 = quantile(prob_AF, probs = 0.975, na.rm = TRUE)
      ) |>
      right_join(predictions_summary)
  }
  
  saveRDS(predictions_summary, file = file.path(dir_i, "predictions_summary.rds"))

  # data
  # NOTE: separate summary since data in predictions is multiplied (by draws sample size)
  data_summary <- H_all |>
    filter(!control) |>
    group_by(group_number, secondai_ID) |>
    summarise(
      D_h = sum(D_h),
      N_h = sum(N_h),
      A_h = sum(A),
      AF_h = sum(AF),
      Trial_code = paste(unique(Trial_code), collapse = ', '),
      country = paste(unique(country), collapse = ', '),
      site = paste(unique(site), collapse = ', '),
      year = paste(unique(year), collapse = ', '),
      insecticide = paste(unique(insecticide), collapse = ', '),
      hut_type = paste(unique(hut_type), collapse = ', '),
      .groups = "keep"
    ) |>
    mutate(
      insecticide = case_when(str_detect(insecticide, ",") ~ "pyrethroid",
                              TRUE ~ insecticide),
      MLE = D_h / N_h, # note this is equivalent to the mode of beta(D_h + 1, N_h - D_h +1)
      q025 = qbeta(0.025, shape1 = D_h + 1, shape2 = N_h - D_h + 1),
      q975 = qbeta(0.975, shape1 = D_h + 1, shape2 = N_h - D_h + 1),
      MLE_f = AF_h / A_h, # note this is equivalent to the mode of beta(AF_h + 1, A_h - AF_h +1)
      q025_f = qbeta(0.025, shape1 = AF_h + 1, shape2 = A_h - AF_h + 1),
      q975_f = qbeta(0.975, shape1 = AF_h + 1, shape2 = A_h - AF_h + 1)
    )
  
  # merge both
  aVSp_summary <- data_summary |>
    left_join(predictions_summary)
  

# plotting
################################
# fix aes scales for all plots
shapes4insecticides <- c("none"= "plus", "pyrethroid" = "circle filled", "permethrin" = "square filled", 'deltamethrin' = "diamond filled", 'alphacypermethrin' = "triangle filled", 'lambdacyhalothrin' = "triangle down filled")
shapes4huttype <- c('East' = 3, 'West' = 4, 'Ifakara' = 8)
allcountries_sort = sort(unique(B_all$country))

# plot mortality
p_H <- ggplot(aVSp_summary) +
  scale_shape_manual(values = shapes4insecticides) +
  # scale_color_manual(values = colours4country) +
  geom_abline(intercept = 0, slope = 1, color = 'grey') +
  geom_point(aes(x = MLE, y = pred_median, colour = factor(country, levels = allcountries_sort), fill = factor(country, levels = allcountries_sort), shape = insecticide)) +
  geom_linerange(aes(y = pred_median, xmin = q025, xmax = q975, colour = factor(country, levels = allcountries_sort))) +
  geom_linerange(aes(x = MLE, ymin = pred_q025, ymax = pred_q975, colour = factor(country, levels = allcountries_sort))) +
  geom_smooth(aes(x = MLE, y = pred_median), alpha = 0.6, colour = "black", method="lm", se=FALSE) +
  xlim(0,1) + ylim(0,1) +
  theme(legend.position="right", text=element_text(size=9), aspect.ratio=1) +
  ylab("Predicted EHT mortality") + xlab("Actual EHT mortality")

ggsave(file = file.path(dir_i, "LOGO-CV-actVSpred.png"), plot=p_H, width = 8, height = 8)
saveRDS(p_H, file = file.path(dir_i, "LOGO-CV-actVSpred.rds"))


# plot feeding if applicable
if (with_feeding) {
  p_H_f <- ggplot(aVSp_summary) +
    scale_shape_manual(values = shapes4insecticides) +
    # scale_color_manual(values = colours4country) +
    geom_abline(intercept = 0, slope = 1, color = 'grey') +
    geom_point(aes(x = MLE, y = pred_f_median, colour = factor(country, levels = allcountries_sort), fill = factor(country, levels = allcountries_sort), shape = insecticide)) +
    geom_linerange(aes(y = pred_f_median, xmin = q025, xmax = q975, colour = factor(country, levels = allcountries_sort))) +
    geom_linerange(aes(x = MLE, ymin = pred_f_q025, ymax = pred_f_q975, colour = factor(country, levels = allcountries_sort))) +
    geom_smooth(aes(x = MLE, y = pred_f_median), alpha = 0.6, colour = "black", method="lm", se=FALSE) +
    xlim(0,1) + ylim(0,1) +
    theme(legend.position="right", text=element_text(size=9), aspect.ratio=1) +
    ylab("Predicted EHT feeding out of alive") + xlab("Actual EHT feeding out of alive")
  
  ggsave(file = file.path(dir_i, "LOGO-CV-actVSpred_feeding.png"), plot=p_H_f, width = 8, height = 8)
  saveRDS(p_H_f, file = file.path(dir_i, "LOGO-CV-actVSpred_feeding.rds"))
}



