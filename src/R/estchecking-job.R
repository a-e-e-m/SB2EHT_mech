# config and prepare
###############################
###############################

#libraries
options(mc.cores = 7)
library("loo")
library("bayesplot")
#library("Rfast")
library("HDInterval")
library("tidyr")
library("dplyr")
library("reshape2")
library("posterior")
library("readr")
library("rstan")
library("stringr")
library("patchwork")
library("data.table")
rstan_options(auto_write = TRUE) # set FALSE for cluster
library("ggrepel")

# source script containing needed functions
functionsfolder <- file.path('./src/R/functions')
files.source <- list.files(functionsfolder)
invisible(sapply(files.source, function(x) source(paste0(functionsfolder, "/", x))))


models2run <- read.csv(file = './models2run.csv', sep = ';', stringsAsFactors = FALSE)

# get model index from argument passed to this script
args = commandArgs(trailingOnly=TRUE)
i <- as.numeric(args[1])

# print model specification to output 
models2run[i,]

# postprocessing
###############################
###############################

# loading data
#############
with_feeding <- as.logical(models2run$with_feeding[i])
BA_only <- as.logical(models2run$BA_only[i])
dir_i <- paste0("fitting/", models2run$run[i]) # folder
est_i <- readRDS(file.path(dir_i, 'est.rds'))


# ## posterior predictive checks
# ######################################################

feature_to_color_BA <- "factor(group_number)"
feature_to_color_EHT <- "factor(group_number)"

# load data (B, H, H_f)
B_all <- readRDS(file.path(dir_i, 'B_all.rds'))

if (!BA_only){
  H_all <- readRDS(file.path(dir_i, 'H_all.rds'))
}else{H_all <- tibble()}

if (with_feeding){
  H_f_all <- readRDS(file.path(dir_i, 'H_f_all.rds'))
}else{H_f_all <- tibble()}


if (str_detect(models2run$data.generating.process[i], 'beta_binomial')){
  # extract probabilities 
  param_D_b <- est_i[["par"]][["phi_D_b"]]
  param_D_h <- est_i[["par"]][["phi_D_h"]]
  param_AF <- est_i[["par"]][["phi_AF"]]
  
  # sampling from beta_binomial...
  
} else if (str_detect(models2run$data.generating.process[i], 'binomial')){
  
  inv_logit <- function(x){
    1 / ( 1 + exp( -x ))
  }
  
  # get probabilities, inv_logit transform if necessary
  if (any(str_detect(names(est_i[["par"]]), "prob_D_b" ))){
    prob_D_b <- est_i[["par"]][["prob_D_b"]]
  } else if (any(str_detect(names(est_i[["par"]]), "problogit_D_b" ))){
    prob_D_b <- inv_logit( est_i[["par"]][["problogit_D_b"]] )
  } else {
    prob_D_b <- numeric(0)
  }
  
  if (any(str_detect(names(est_i[["par"]]), "prob_D_h" ))){
    prob_D_h <- est_i[["par"]][["prob_D_h"]]
  } else if (any(str_detect(names(est_i[["par"]]), "problogit_D_h" ))){
    prob_D_h <- inv_logit( est_i[["par"]][["prob_D_h"]] )
  } else {
    prob_D_h <- numeric(0)
  }
  
  if (any(str_detect(names(est_i[["par"]]), "prob_AF" ))){
    prob_AF <- est_i[["par"]][["prob_AF"]]
  } else if (any(str_detect(names(est_i[["par"]]), "problogit_AF" ))){
    prob_AF <- inv_logit( est_i[["par"]][["prob_AF"]] )
  } else {
    prob_AF <- numeric(0)
  }

}


# fix aes scales for all plots
shapes4insecticides <- c("none"= "plus", "pyrethroid" = "circle filled", "permethrin" = "square filled", 'deltamethrin' = "diamond filled", 'alphacypermethrin' = "triangle filled", 'lambdacyhalothrin' = "triangle down filled")
shapes4huttype <- c('East' = 3, 'West' = 4, 'Ifakara' = 8)
allcountries_sort = sort(unique(B_all$country))

#initialise empty plots
p_B <- ggplot()
p_H <- ggplot()
p_H_f <- ggplot()


# summaries data and add estimates, and plot
if (length(prob_D_b) !=0){
  B_all_summary <- B_all |>
    cbind(tibble(prob_D_b)) |>
    group_by(group_number, control, treat, times_disc_dose) |>
    # be warned: this is not fully general but probably covers all group_by variants wanted
    summarise(D_b = sum(D_b),
              N_b = sum(N_b),
              Trial_code = paste(unique(Trial_code), collapse = ', '),
              country = paste(unique(country), collapse = ', '),
              site = paste(unique(site), collapse = ', '),
              year = paste(unique(year), collapse = ', '),
              insecticide = paste(unique(insecticide), collapse = ', '),
              test_type = paste(unique(test_type), collapse = ', '),
              pred = mean(prob_D_b),
              .groups = "keep") |>
    mutate(insecticide = case_when( str_detect(insecticide, ",") ~ "pyrethroid",
                                   TRUE ~ insecticide),
           MLE = D_b / N_b, # ntote this is equivalent to the mode of beta(D_b + 1, N_b - D_b +1)
           q025 = qbeta(0.025, shape1 = D_b + 1, shape2 = N_b - D_b + 1),
           q975 = qbeta(0.975, shape1 = D_b + 1, shape2 = N_b - D_b + 1),
           hdi_l = hdi(qbeta, 0.95, shape1 = D_b + 1, shape2 = N_b - D_b + 1)[[1]],
           hdi_r = hdi(qbeta, 0.95, shape1 = D_b + 1, shape2 = N_b - D_b + 1)[[2]]
           )
  
  p_B <- ggplot(B_all_summary) +
    scale_shape_manual(values = shapes4insecticides) +
    # scale_color_manual(values = colours4country) +
    geom_abline(intercept = 0, slope = 1, color = 'grey') +
    geom_point(aes(x = MLE, y = pred, colour = factor(country, levels = allcountries_sort), fill = factor(country, levels = allcountries_sort), shape = insecticide)) +
    geom_label_repel(data = B_all_summary |> filter(times_disc_dose > 1), aes(x = MLE, y = pred, label = times_disc_dose, colour = country), size = 2, box.padding = 0.2, point.size = 5, max.overlaps = Inf, alpha = 0.8) + # use box.padding > 0 to draw lines from points to text
    geom_linerange(aes(y = pred, xmin = hdi_l, xmax = hdi_r, colour = factor(country, levels = allcountries_sort))) +
    geom_smooth(aes(x = MLE, y = pred), alpha = 0.6, colour = "black", method="lm", se=FALSE) +
    xlim(0,1) + ylim(0,1) +
    #coord_fixed(ratio = 1, xlim = c(0,1), ylim = c(0,1)) +
    facet_wrap( ~ treat, labeller = label_both) +
    theme(legend.position="right", text=element_text(size=9), aspect.ratio=1) +
    ylab("Predicted BA mortality") + xlab("Actual BA mortality")
}
  

if (length(prob_D_h) !=0){  
  H_all_summary <- H_all |>
    cbind(tibble(prob_D_h)) |>
    group_by(group_number, control, treat) |>
    # be warned: this is not fully general but probably covers all group_by variants wanted
    summarise(D_h = sum(D_h),
              N_h = sum(N_h),
              Trial_code = paste(unique(Trial_code), collapse = ', '),
              country = paste(unique(country), collapse = ', '),
              site = paste(unique(site), collapse = ', '),
              year = paste(unique(year), collapse = ', '),
              insecticide = paste(unique(insecticide), collapse = ', '),
              hut_type = paste(unique(hut_type), collapse = ', '),
              pred = mean(prob_D_h),
              .groups = "keep"
    ) |>
    mutate(insecticide = case_when( str_detect(insecticide, ",") ~ "pyrethroid",
                                    TRUE ~ insecticide),
           MLE = D_h / N_h, # note this is equivalent to the mode of beta(D_h + 1, N_h - D_h +1)
           q025 = qbeta(0.025, shape1 = D_h + 1, shape2 = N_h - D_h + 1),
           q975 = qbeta(0.975, shape1 = D_h + 1, shape2 = N_h - D_h + 1),
           hdi_l = hdi(qbeta, 0.95, shape1 = D_h + 1, shape2 = N_h - D_h + 1)[1],
           hdi_r = hdi(qbeta, 0.95, shape1 = D_h + 1, shape2 = N_h - D_h + 1)[2]
    ) 
  
  p_H <- ggplot(H_all_summary) +
    scale_shape_manual(values = shapes4insecticides) +
    # scale_color_manual(values = colours4country) +
    geom_abline(intercept = 0, slope = 1, color = 'grey') +
    geom_point(aes(x = MLE, y = pred, colour = factor(country, levels = allcountries_sort), fill = factor(country, levels = allcountries_sort), shape = insecticide)) +
    geom_linerange(aes(y = pred, xmin = hdi_l, xmax = hdi_r, colour = factor(country, levels = allcountries_sort))) +
    geom_smooth(aes(x = MLE, y = pred), alpha = 0.6, colour = "black", method="lm", se=FALSE) +
    xlim(0,1) + ylim(0,1) +
    facet_wrap( ~ treat, labeller = label_both) +
    theme(legend.position="right", text=element_text(size=9), aspect.ratio=1) +
    ylab("Predicted EHT mortality") + xlab("Actual EHT mortality")
}
  
if (length(prob_AF) !=0){
  H_f_all_summary <- H_f_all |>
    cbind(tibble(prob_AF)) |>
    group_by(group_number, control, treat) |>
    # be warned: this is not fully general but probably covers all group_by variants wanted
    summarise(AF = sum(AF),
              A = sum(A),
              Trial_code = paste(unique(Trial_code), collapse = ', '),
              country = paste(unique(country), collapse = ', '),
              site = paste(unique(site), collapse = ', '),
              year = paste(unique(year), collapse = ', '),
              insecticide = paste(unique(insecticide), collapse = ', '),
              hut_type = paste(unique(hut_type), collapse = ', '),
              pred = mean(prob_AF),
              .groups = "keep"
    ) |>
    mutate(insecticide = case_when( str_detect(insecticide, ",") ~ "pyrethroid",
                                    TRUE ~ insecticide),
           MLE = AF / A, # note this is equivalent to the mode of beta(AF + 1, A - AF +1)
           q025 = qbeta(0.025, shape1 = AF + 1, shape2 = A - AF + 1),
           q975 = qbeta(0.975, shape1 = AF + 1, shape2 = A - AF + 1),
           hdi_l = hdi(qbeta, 0.95, shape1 = AF + 1, shape2 = A - AF + 1)[1],
           hdi_r = hdi(qbeta, 0.95, shape1 = AF + 1, shape2 = A - AF + 1)[2]
    )
  
  p_H_f <- ggplot(H_f_all_summary) +
    scale_shape_manual(values = shapes4insecticides) +
    # scale_color_manual(values = colours4country) +
    geom_abline(intercept = 0, slope = 1, color = 'grey') +
    geom_point(aes(x = MLE, y = pred, colour = factor(country, levels = allcountries_sort), fill = factor(country, levels = allcountries_sort), shape = insecticide)) +
    geom_linerange(aes(y = pred, xmin = hdi_l, xmax = hdi_r, colour = factor(country, levels = allcountries_sort))) +
    geom_smooth(aes(x = MLE, y = pred), alpha = 0.6, colour = "black", method="lm", se=FALSE) +
    xlim(0,1) + ylim(0,1) +
    facet_wrap( ~ treat, labeller = label_both) +
    theme(legend.position="right", text=element_text(size=9), aspect.ratio=1) +
    ylab("Predicted fed-alive out of alive in EHT") + xlab("Actual fed-alive out of alive in EHT")
}
  

# combine plots
plot_comb <- p_B / p_H / p_H_f +
  plot_layout(guides = "collect") # & scale_color_discrete(drop = F)

ggsave(file.path(dir_i, "postpc_comb_pointest.png"), plot_comb, width = 10, height = 10)






