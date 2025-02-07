# this script plots for the mechanistic model with the alpha-re-parameterisation the following plots
# - actual vs predicted reference line type plot for BA and EHT mortality
# - BA dose-response data and prediction 



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
library("patchwork")
library("data.table")
#library("ggrepel")
library("scales")
library("ggplot2")
library("RColorBrewer")

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

# loading specifications
#############
with_feeding <- as.logical(models2run$with_feeding[i])
BA_only <- as.logical(models2run$BA_only[i])
dir_i <- paste0("fitting/", models2run$run[i]) # folder
trials <- readRDS(file.path(dir_i, 'trials.rds'))


# loading real data
#############
B_all <- readRDS(file.path(dir_i, 'B_all.rds')) |>
  mutate(country = case_when(country == "BurkinaFaso" ~ "Burkina Faso",
                             TRUE ~ country))

if (!BA_only){
  H_all <- readRDS(file.path(dir_i, 'H_all.rds')) |>
    mutate(country = case_when(country == "BurkinaFaso" ~ "Burkina Faso",
                               TRUE ~ country))
}else{H_all <- tibble()}

if (with_feeding){
  H_f_all <- readRDS(file.path(dir_i, 'H_f_all.rds')) |>
    mutate(country = case_when(country == "BurkinaFaso" ~ "Burkina Faso",
                               TRUE ~ country))
}else{H_f_all <- tibble()}


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


if (!BA_only){
  mu_x <- rstan::extract(fit_i, pars = "mu_x", permuted = TRUE)[["mu_x"]][index_sample]
  sigma_x <- rstan::extract(fit_i, pars = "sigma_x", permuted = TRUE)[["sigma_x"]][index_sample]
  # switch needed because p_b_control sometimes group specific and sometimes not
  if (length(dim( rstan::extract(fit_i, pars = "p_h_control", permuted = TRUE)[["p_h_control"]])) == 1){
    p_h_control <- rstan::extract(fit_i, pars = "p_h_control", permuted = TRUE)[["p_h_control"]][index_sample]
  } else {
    p_h_control <- rstan::extract(fit_i, pars = "p_h_control", permuted = TRUE)[["p_h_control"]][index_sample,]
  }
}else{H_all <- tibble()}



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
  
  if ("prob_D_h" %in% fit_i@model_pars){
    prob_D_h <-  rstan::extract(fit_i, pars = "prob_D_h")[["prob_D_h"]][index_sample,]
  } else if ("problogit_D_h" %in% fit_i@model_pars){
    prob_D_h <- inv_logit( rstan::extract(fit_i, pars = "problogit_D_h")[["problogit_D_h"]][index_sample,])
  } else {
    prob_D_h <- numeric(0)
  }

  
# data processing and plotting
############################
  
# fix aes scales for all plots
shapes4insecticides <- c("none"= "plus", "pyrethroid" = "circle filled", "permethrin" = "square filled", 'deltamethrin' = "diamond filled", 'alphacypermethrin' = "triangle filled", 'lambdacyhalothrin' = "triangle down filled")
shapes4insecticides_treat_publication <- c(
  'alphacypermethrin' = "circle",
  'deltamethrin' = "triangle",
  "permethrin" = "square",
  "pyrethroid" = "plus"
)
shapes4huttype <- c('East' = 3, 'West' = 4, 'Ifakara' = 8)
allcountries_sort = sort(unique(B_all$country))


# summaries data and add estimates, and plot for BA and EHT
if (length(prob_D_b) !=0){
  # connect model estimates of death probability and real data, sumarise
  # NOTE: this uses the prob estimates from the stan output directly, not using the other parameters
    col_D_b <- paste("V", seq(1, nrow(B_all),1), sep = "")
    prob_D_b <- as.data.table(prob_D_b)[, iter := .I]
    prob_D_b <- melt(prob_D_b, measure.vars = col_D_b, variable.name = "data_point", value.name = c("prob"))
    prob_D_b <- prob_D_b[, `:=`("data_point" = as.integer(str_remove(data_point, "V")))]
    
    B_all <- as.data.table(B_all)[, data_point := .I]
    setkey(B_all, ON = data_point)
    setkey(prob_D_b, ON = data_point)
    prob_D_b <- B_all[prob_D_b, nomatch = 0]
    
    prob_D_b_summary <- prob_D_b[
      ,
      list(
        "pred_median" = median(prob),
        "pred_hdi_l" = hdi(prob)[["lower"]],
        "pred_hdi_r" = hdi(prob)[["upper"]],
        "pred_q025" = quantile(prob, probs = 0.025), 
        "pred_q975" = quantile(prob, probs = 0.975)
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
  
  # reference line actual vs predicted plot
  ###################
    p_B <- ggplot(B_all_summary) +
      scale_fill_brewer(palette = "Dark2") +
      scale_color_brewer(palette = "Dark2") +
      scale_shape_manual(name = 'Insecticide', values = shapes4insecticides) +
      geom_abline(intercept = 0, slope = 1, color = 'grey') +
      geom_smooth(aes(x = MLE, y = pred_median), alpha = 0.6, colour = "black", method="lm", se=FALSE) +
      geom_point(aes(x = MLE, y = pred_median, colour = factor(country, levels = allcountries_sort), shape = factor(insecticide, levels = names(shapes4insecticides_treat_publication)))) +
      geom_linerange(aes(y = pred_median, xmin = q025, xmax = q975, colour = factor(country, levels = allcountries_sort))) +
      geom_linerange(aes(x = MLE, ymin = pred_q025, ymax = pred_q975, colour = factor(country, levels = allcountries_sort))) +
      xlim(0,1) + ylim(0,1) +
      facet_wrap( ~ treat, labeller = label_both) +
      ylab("Predicted SB mortality [Probability]") + xlab("Actual SB mortality [Probability]") +
      scale_x_continuous(labels = percent) +
      scale_y_continuous(labels = percent)+
      ggtitle("Bio assay") + 
      labs(colour = "Country", shape = "Insecticide")
    
    p_B_treat <- ggplot(B_all_summary |> filter(treat == 1)) +
      scale_fill_brewer(palette = "Dark2") +
      scale_color_brewer(palette = "Dark2") +
      scale_shape_manual(name = 'Insecticide', values = shapes4insecticides_treat_publication, drop = F) +
      geom_abline(intercept = 0, slope = 1, color = 'grey') +
      geom_smooth(aes(x = MLE, y = pred_median), alpha = 0.6, colour = "black", method="lm", se=FALSE) +
      geom_point(aes(x = MLE, y = pred_median, colour = factor(country, levels = allcountries_sort), shape = factor(insecticide, levels = names(shapes4insecticides_treat_publication))), show.legend = TRUE) +
      geom_linerange(aes(y = pred_median, xmin = q025, xmax = q975, colour = factor(country, levels = allcountries_sort))) +
      geom_linerange(aes(x = MLE, ymin = pred_q025, ymax = pred_q975, colour = factor(country, levels = allcountries_sort))) +
      ylab("Predicted SB mortality [Probability]") + xlab("Actual SB mortality [Probability]") +
      scale_x_continuous(labels = percent) +
      scale_y_continuous(labels = percent) +
      ggtitle("Bio assay") + 
      labs(colour = "Country", shape = "Insecticide")
    
    # # add labels for times_disc_dose if intensity dose data
    # if (length(unique(B_all_summary$times_disc_dose)) > 1 ){
    #   p_B <- p_B + geom_label_repel(data = B_all_summary |> filter(times_disc_dose > 1), aes(x = MLE, y = pred_median, label = times_disc_dose, colour = country), size = 2, box.padding = 0.2, point.size = 5, max.overlaps = Inf, alpha = 0.8) # use box.padding > 0 to draw lines from points to text
    #   p_B_treat <- p_B_treat + geom_label_repel(data = B_all_summary |> filter(treat == 1, times_disc_dose > 1), aes(x = MLE, y = pred_median, label = times_disc_dose, colour = country), size = 2, box.padding = 0.2, point.size = 5, max.overlaps = Inf, alpha = 0.8) # use box.padding > 0 to draw lines from points to text
    # }
  
  
  # dose response curves per group with facets (only for groups with intensity dose data)
  ###################
    # compute predictions for dose-response
    times_disc_dose <- seq(0, 15, 0.05)
    list_dfs_doseresponse <- vector(mode = "list", length = trials[["nT"]])
    list_dfs_doseresponse_summary <- vector(mode = "list", length = trials[["nT"]])
    
    # lookup for group_number for int_dose groups
    # contains quick hack to exclude data from Tanzania
    lookup_group_number <- B_all_summary |>
      filter(int_dose_available == 1, country == "Burkina Faso") |>
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
    
    # generate index sample for plotting sample of dose-response curves only
    iter_sample <- sample(1:max(list_dfs_doseresponse[[1]]$iter), 500, replace = FALSE)

  
    # initialise plot
    p_B_facet <- ggplot() +
      scale_fill_brewer(palette = "Dark2") +
      scale_color_brewer(palette = "Dark2")
    
    # add dose-response curves for each group
    for (ii in seq_along(sort(unique(lookup_group_number$group_number)))){
        p_B_facet <- p_B_facet +
          geom_line(data = list_dfs_doseresponse[[ii]] |> filter(iter %in% iter_sample), aes(x = times_disc_dose, y = p, group = iter, colour = insecticide), alpha = 0.01) +
          geom_line(data = list_dfs_doseresponse_summary[[ii]], aes(x = times_disc_dose, y = p_median, colour = insecticide), linewidth = 0.5) 
    }
    
    
    # add real data and empirical estimates 
    # as samples for violin plots are of equal size, the areas should match across facets, otherwise check
    # this post for making violin area uniform across all facets: https://stackoverflow.com/questions/47174825/same-area-for-all-violins-independent-of-facets-in-ggplot2
    p_B_facet <- p_B_facet +
      # contains quick hack to exclude data from Tanzania
      geom_pointrange(data = B_all_summary |> filter(int_dose_available == 1, country == "Burkina Faso"), mapping = aes(x = times_disc_dose, y = median, ymin = q025, ymax = q975, colour = insecticide), size = 0.2, alpha = 0.9) +
      geom_jitter(data = B_all |> filter(int_dose_available == 1, country == "Burkina Faso"), aes(x = times_disc_dose, y = D_b / N_b, colour = insecticide), size = 0.5, alpha = 0.9) +
      facet_wrap(~ group_number + site + year , labeller = 
                   labeller(
                     site = ~ paste(.),
                     year = ~ paste(.),
                     group_number = ~ paste(.),
                     .multi_line = FALSE
                   )
                 ) +
      scale_y_continuous(labels = percent) +
      xlab("Dose (as multiples of discriminating dose)") + ylab("Mortality") +
      labs(colour = "Insecticide", title = NULL) 
  
    # save that last plot
    ggsave(file.path(dir_i, "dose-response-curves-BA.png"), p_B_facet, width = sqrt(B_all |> filter(int_dose_available == 1)  |> pull(group_number) |> unique() |> length()) * 2 + 3, height = sqrt(B_all |> filter(int_dose_available == 1) |> pull(group_number) |> unique() |> length()) * 2)
}

if (length(prob_D_h) !=0){    
  # sumarise probability estimates
  # data.table 
  # TODO: Turn this into a function
  # H
  col_D_h <- paste("V", seq(1, nrow(H_all),1), sep = "")
  prob_D_h <- as.data.table(prob_D_h)[, iter := .I]
  prob_D_h <- melt(prob_D_h, measure = col_D_h, variable.name = "data_point", value.name = c("prob"))
  prob_D_h <- prob_D_h[, `:=`("data_point" = as.integer(str_remove(data_point, "V")))]
  
  H_all <- as.data.table(H_all)[, data_point := .I]
  setkey(H_all, ON = data_point)
  setkey(prob_D_h, ON = data_point)
  prob_D_h <- H_all[prob_D_h, nomatch = 0]
  
  prob_D_h_summary <- prob_D_h[
    ,
    list(
      "pred_median" = median(prob),
      "pred_hdi_l" = hdi(prob)[["lower"]],
      "pred_hdi_r" = hdi(prob)[["upper"]],
      "pred_q025" = quantile(prob, probs = 0.025), 
      "pred_q975" = quantile(prob, probs = 0.975)
    ),
    .(group_number, control, treat)
  ]
  
  H_all_summary <- H_all |>
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
              .groups = "keep"
    ) |>
    mutate(insecticide = case_when( str_detect(insecticide, ",") ~ "pyrethroid",
                                    TRUE ~ insecticide),
           MLE = D_h / N_h, # note this is equivalent to the mode of beta(D_h + 1, N_h - D_h +1)
           q025 = qbeta(0.025, shape1 = D_h + 1, shape2 = N_h - D_h + 1),
           q975 = qbeta(0.975, shape1 = D_h + 1, shape2 = N_h - D_h + 1),
           hdi_l = hdi(qbeta, 0.95, shape1 = D_h + 1, shape2 = N_h - D_h + 1)[1],
           hdi_r = hdi(qbeta, 0.95, shape1 = D_h + 1, shape2 = N_h - D_h + 1)[2]
    )  |> 
    left_join(prob_D_h_summary, by = c("group_number", "control", "treat"))
  
  p_H <- ggplot(H_all_summary) +
    scale_fill_brewer(palette = "Dark2") +
    scale_color_brewer(palette = "Dark2") +
    scale_shape_manual(name = 'Insecticide', values = shapes4insecticides) +
    geom_abline(intercept = 0, slope = 1, color = 'grey') +
    geom_smooth(aes(x = MLE, y = pred_median), alpha = 0.6, colour = "black", method="lm", se=FALSE) +
    geom_point(aes(x = MLE, y = pred_median, colour = factor(country, levels = allcountries_sort), shape = factor(insecticide, levels = names(shapes4insecticides_treat_publication)))) +
    geom_linerange(aes(y = pred_median, xmin = q025, xmax = q975, colour = factor(country, levels = allcountries_sort))) +
    geom_linerange(aes(x = MLE, ymin = pred_q025, ymax = pred_q975, colour = factor(country, levels = allcountries_sort))) +
    xlim(0,1) + ylim(0,1) +
    facet_wrap( ~ treat, labeller = label_both) +
    ylab("Predicted EHT mortality [Probability]") + xlab("Actual EHT mortality [Probability]") +
    scale_x_continuous(labels = percent) +
    scale_y_continuous(labels = percent) +
    ggtitle("Experimental hut trial") + 
    labs(colour = "Country", shape = "Insecticide")
  
  p_H_treat <- ggplot(H_all_summary |> filter(treat == 1)) +
    scale_fill_brewer(palette = "Dark2") +
    scale_color_brewer(palette = "Dark2") +
    scale_shape_manual(name = 'Insecticide', values = shapes4insecticides_treat_publication, drop = F) +
    geom_abline(intercept = 0, slope = 1, color = 'grey') +
    geom_smooth(aes(x = MLE, y = pred_median), alpha = 0.6, colour = "black", method="lm", se=FALSE) + 
    geom_point(aes(x = MLE, y = pred_median, colour = factor(country, levels = allcountries_sort),  shape = factor(insecticide, levels = names(shapes4insecticides_treat_publication))), show.legend = TRUE) +
    geom_linerange(aes(y = pred_median, xmin = q025, xmax = q975, colour = factor(country, levels = allcountries_sort))) +
    geom_linerange(aes(x = MLE, ymin = pred_q025, ymax = pred_q975, colour = factor(country, levels = allcountries_sort))) +
    ylab("Predicted EHT mortality [Probability]") + xlab("Actual EHT mortality [Probability]") +
    scale_x_continuous(labels = percent) +
    scale_y_continuous(labels = percent) +
    ggtitle("Experimental hut trial") + 
    labs(colour = "Country", shape = "Insecticide")
}  

# plot actual vs predicted
if (length(prob_D_h) !=0){    
  # combine plots
  plot_comb <- p_B / p_H +
    plot_layout(guides = "collect") & 
    theme(legend.position="bottom", text=element_text(size=9), aspect.ratio=1) 
  
  plot_comb_treat <- p_B_treat + p_H_treat +
    plot_annotation(tag_levels = 'A') +
    plot_layout(guides = "collect") & 
    scale_shape_discrete(drop = F) &
    theme(legend.position="bottom", text=element_text(size=9), aspect.ratio=1) & 
    coord_fixed(ratio = 1, xlim = c(0,1), ylim = c(0,1)) &
    labs(title = NULL) &
    guides(colour = guide_legend(order = 1, nrow = 2), 
           shape = guide_legend(order = 2, nrow = 2)) 
  
  #ggsave(file.path(dir_i, "postpc_comb.png"), plot_comb, width = 10, height = 10)
  ggsave(file.path(dir_i, "postpc_comb_treat_quantile.png"), plot_comb_treat, width = 10, height = 6)
} else {
  ggsave(file.path(dir_i, "postpc_SBonly_treat_quantile.png"), p_B_treat, width = 10, height = 6)
}





