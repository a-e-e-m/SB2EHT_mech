# prior predictive checks for mechanistic models
library("tidyverse")
library("data.table")
library("HDInterval")
library("scales")
library("patchwork")

## mechanistic model without correlation
### only BA 
# nested_binomial_mechanistic_1_separated_p_b_control_onlyBA.stan
 # set sample sizes:
  nT <- 1
  param_ss <- 2000
  
  # setup simulation 
  T_b <- 1:nT
  iter <- 1:param_ss
  
  #function for half normal
  rhalfnorm <- function(size, mu, sigma){
    sample <- numeric()
    while (length(sample) < size){
      draw <- rnorm(1, mu, sigma)
      if (draw > 0){
        sample <- c(sample, draw)
      }
    }
    return(sample)
  }
  
  rhalfcauchy <- function(size, mu, sigma){
    sample <- numeric()
    while (length(sample) < size){
      draw <- rcauchy(1, mu, sigma)
      if (draw > 0){
        sample <- c(sample, draw)
      }
    }
    return(sample)
  }
 
  # sample parameter values from prior
  mu_d <- rnorm(param_ss, 5, 10)
  sigma_d <- rhalfnorm(param_ss, 0, 5)
  sigma_v <- rhalfnorm(param_ss, 0, 5)
  mu_x <- rnorm(param_ss, 0, 5)
  sigma_x <- rhalfnorm(param_ss, 0, 5)
  p_b_control <- rbeta(param_ss, 1,10)
  p_h_control <- rbeta(param_ss, 1,10)
  
  df <- tibble(iter, mu_d, sigma_d, sigma_v, mu_x, sigma_x, p_b_control, p_h_control)
  
  ## plot ID-BA
  ###########################
  c <- seq(0, 30, 0.05)
  dt_sim <- as.data.table(expand_grid(df, c))
  
  # compute probabilities
  dt_sim <- dt_sim[, p := p_b_control + ( 1 - p_b_control) * pnorm( (log(c) - mu_d ) / sqrt(sigma_v^2 + sigma_d^2 ))]
  
  dt_sim_sum <- tibble(dt_sim) |> 
    group_by(c) |>
    dplyr::summarise(p_median = median(p), 
              p_hdi99_l = hdi(p, 0.99)[[1]],
              p_hdi99_r = hdi(p, 0.99)[[2]],
              p_hdi95_l = hdi(p, 0.95)[[1]],
              p_hdi95_r = hdi(p, 0.95)[[2]],
              p_hdi80_l = hdi(p, 0.8)[[1]],
              p_hdi80_r = hdi(p, 0.8)[[2]],
              p_hdi50_l = hdi(p, 0.5)[[1]],
              p_hdi50_r = hdi(p, 0.5)[[2]],
              q_p0005 = quantile(p, 0.005),
              q_p0025 = quantile(p, 0.025),
              q_p0100 = quantile(p, 0.100),
              q_p0250 = quantile(p, 0.250),
              q_p0500 = quantile(p, 0.500),
              q_p0750 = quantile(p, 0.750),
              q_p0900 = quantile(p, 0.900),
              q_p0975 = quantile(p, 0.975),
              q_p0995 = quantile(p, 0.995),
              .groups = "keep")
  
  # generate index sample for plotting sample of dose-response curves only
  iter_sample <- sample(1:max(dt_sim$iter), 500, replace = FALSE)

  # plotting
  # todo: add labels with prob level to quantile curves, see https://community.rstudio.com/t/label-geom-line-with-a-label/71035/2
  # add geom_line and geom_label_repel or geom_label
  # but this needs to reshape dt_sim with variable for prob level
  
  cols = c("sample curves" = 'blue', "median" = "green")
  
  p <- ggplot(dt_sim |> dplyr::filter(iter %in% iter_sample), aes(x = c)) +
    geom_ribbon(data = dt_sim_sum, aes(ymin = q_p0005, ymax = q_p0995, fill = "99%"), colour = "black", alpha = 0.5) +
    geom_ribbon(data = dt_sim_sum, aes(ymin = q_p0025, ymax = q_p0975, fill = "95%"), colour = "black", alpha = 0.5) +
    geom_ribbon(data = dt_sim_sum, aes(ymin = q_p0100, ymax = q_p0900, fill = "80%"), colour = "black", alpha = 0.5) +
    geom_ribbon(data = dt_sim_sum, aes(ymin = q_p0250, ymax = q_p0750, fill = "50%"), colour = "black", alpha = 0.5) +
    scale_colour_manual(values = cols) + 
    geom_line(aes(x = c, y = p, group = iter, colour = "sample curves"), alpha = 0.05 ) +
    geom_line(data = dt_sim_sum, aes(x = c, y = p_median, colour = "median"), alpha = 0.7, linewidth = 0.75) +
    scale_y_continuous(labels = percent) +
    xlab("Insecticide challenge [Dose]") + ylab("Mortality [Probability]") +
    labs(fill = "Area with probability mass", colour = "Dose-response curves") 
  
  p_log <- ggplot(dplyr::filter(dt_sim, iter %in% iter_sample, c > 0.01), aes(x = c)) +
    geom_ribbon(data = dt_sim_sum |> dplyr::filter(c > 0.01), aes(ymin = q_p0005, ymax = q_p0995, fill = "99%"), colour = "black", alpha = 0.5) +
    geom_ribbon(data = dt_sim_sum |> dplyr::filter(c > 0.01), aes(ymin = q_p0025, ymax = q_p0975, fill = "95%"), colour = "black", alpha = 0.5) +
    geom_ribbon(data = dt_sim_sum |> dplyr::filter(c > 0.01), aes(ymin = q_p0100, ymax = q_p0900, fill = "80%"), colour = "black", alpha = 0.5) +
    geom_ribbon(data = dt_sim_sum |> dplyr::filter(c > 0.01), aes(ymin = q_p0250, ymax = q_p0750, fill = "50%"), colour = "black", alpha = 0.5) +
    scale_colour_manual(values = cols) + 
    geom_line(aes(x = c, y = p, group = iter, colour = "sample curves"), alpha = 0.05 ) +
    geom_line(data = dt_sim_sum |> dplyr::filter(c > 0.01), aes(x = c, y = p_median, colour = "median"), alpha = 0.7, linewidth = 0.75) +
    coord_trans(x = "log") +
    scale_x_continuous(breaks = c(0.01, 0.1, 0.25, 0.5, 1, 2, 5, 10, 20)) +
    scale_y_continuous(labels = percent) +
    xlab("Insecticide challenge [Dose]") + ylab("Mortality [Probability]") +
    labs(fill = "Area with probability mass", colour = "Dose-response curves") 

  
  plot_comb <- p_log + p  +
    plot_annotation(tag_levels = list(c("A", "B"))) +
    plot_layout(guides = "collect") &
    theme(legend.position="bottom")
  
  ggsave(file.path("plots_mechmodel_article_SI", "prior_pred_checks_BA.png"), plot_comb, width = 12, height = 6)
  
  
  ## plot EHT
  ###########################
  
  df_sim_EHT <- expand_grid(df, treat = c(FALSE, TRUE))
  
  # compute probabilities 
  df_sim_EHT <- df_sim_EHT |>
    mutate(prob_D_h = ifelse(treat, p_h_control + ( 1 - p_h_control ) * pnorm( (mu_x - mu_d ) / sqrt(sigma_x^2 + sigma_d^2 ) ), p_h_control )    )
  
  
  # plot
  p_EHT <- ggplot(df_sim_EHT) +
    geom_density(aes(x = prob_D_h)) +
    facet_wrap(~ treat, labeller = labeller(treat = c( `FALSE` = "Control EHT", `TRUE` = "Intervention EHT"))) +
    scale_x_continuous(labels = percent) +
    xlab("Mortality [Probability]") 
  
  ggsave(file.path("plots_mechmodel_article_SI", "prior_pred_checks_EHT.png"), p_EHT)
  
  
  
  ## plot DD-BA
  ###########################
  
  df_sim_DDBA <- expand_grid(df, treat = c(FALSE, TRUE))
  
  # compute probabilities 
  df_sim_DDBA <- df_sim_DDBA |>
    mutate(prob_D_b = ifelse(treat, p_b_control + ( 1 - p_b_control ) * pnorm( (0 - mu_d ) / sqrt(sigma_v^2 + sigma_d^2 ) ), p_b_control )    )
  
  # plot
  p_DDBA <- ggplot(df_sim_DDBA) +
    geom_density(aes(x = prob_D_b)) +
    facet_wrap(~ treat, labeller = labeller(treat = c( `FALSE` = "Control SB", `TRUE` = "DD-SB"))) +
    scale_x_continuous(labels = percent) +
    xlab("Mortality [Probability]") 
  
  ggsave(file.path("plots_mechmodel_article_SI", "prior_pred_checks_DD-SB.png"), p_EHT) 
  
  
  ## combine plots
  ###########################
  
  plot_comb <- p_log + p_DDBA + p_EHT  +
    plot_annotation(tag_levels = list(c("A", "B", "C"))) +
    plot_layout(guides = "collect") &
    theme(legend.position="bottom")
  
  ggsave(file.path("plots_mechmodel_article_SI", "prior_pred_checks.png"), plot_comb, width = 12, height = 4.5)
