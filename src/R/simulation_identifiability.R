 # simulation study to see what dose escalation, in terms of LD (relative), is needed for identifying the lethal dose distribution


library(tidyverse)
library(ggplot2)

## fix a few parameter choices, e.g. 6 settings with variable mu_T and sigma_T
mu_T_values <- c(1, 2, 3) # roughly lowest and highest values from Burkina Faso estimates without identifiability issues, plus double the highest
sigma_T_values <- c(0.1, 1, 2) # roughly lowest and highest values from Burkina Faso estimates, plus a middle value

# full factorial combination
param_combinations <- expand_grid(
  mu_T = mu_T_values,
  sigma_T = sigma_T_values
) |>
  mutate(group_number = row_number())

## get the posterior median for sigma_v from the calibrated SB only model
sigma_SB <- 0.1293

## define a control mortality. take average over estimates from Burkina Faso
p_b_control <- 0.016

## produce the dose response curves and plot them

times_disc_dose <- seq(0, 50, 0.05)

# prepare data
# create one data frame per parameter combination
df_list <- param_combinations |>
  pmap(\(mu_T, sigma_T, group_number) {
    tibble(
      group_number = group_number,
      mu_T = mu_T,
      sigma_T = sigma_T,
      c = times_disc_dose
    ) |>
      mutate(
        p = p_b_control +
          (1 - p_b_control) *
          pnorm(
            (log(c) - mu_T) /
              sqrt(sigma_SB^2 + sigma_T^2)
          )
      )
  })

df_all <- bind_rows(df_list)


## compute for each curve LD10, LD20, LD40, LD60, LD80, LD90
LD_probs <- c(0.10, 0.20, 0.40, 0.60, 0.80, 0.90)

df_LD <- param_combinations |>
  tidyr::crossing(LD_prob = LD_probs) |>
  mutate(
    LD = exp(
      mu_T +
        sqrt(sigma_SB^2 + sigma_T^2) *
        qnorm((LD_prob - p_b_control) / (1 - p_b_control))
    ),
    LD_label = paste0("LD", LD_prob * 100)
  )


# two beta parameters for the dispersion setting, one representing binomial and one elevated (look at data)
rho_BB <- 0.1   
rho_values <- c(0, rho_BB)

# total per replicate
N_b <- 50       # set this yourself

# set seed
set.seed(27)

# simulate 5 mosquito death counts for control and all LDs, with total N_b each, corresponding to 5 replicates
rbetabinom <- function(p, N, rho) {
  if (rho == 0) return(rbinom(1, N, p))
  
  phi <- 1 / rho - 1
  rbinom(1, N, rbeta(1, p * phi, (1 - p) * phi))
}

df_sim <- crossing(
  LD_prob = LD_probs,
  rho = rho_values,
  rep = 1:5
) |>
  mutate(
    N_b = N_b,
    D_b = pmap_int(
      list(LD_prob, N_b, rho),
      \(p, N, rho) rbetabinom(p, N, rho)
    )
  )

df_LD_sim <- df_LD |>
  left_join(df_sim, by = "LD_prob")

# control counts
# simulate control outcomes once per rho / replicate
df_control_sim <- crossing(
  rho = rho_values,
  rep = 1:5
) |>
  mutate(
    N_b = N_b,
    D_b = map2_int(rho, N_b, \(rho, N) 
                   rbetabinom(p_b_control, N, rho))
  )

df_control_sim <- param_combinations |>
  crossing(df_control_sim) |>
  mutate(
    LD_prob = p_b_control,
    LD_label = "control",
    LD = 0
  )

df_LD_sim <- bind_rows(df_LD_sim, df_control_sim)

# add labels
df_all <- df_all |>
  mutate(
    mu_label = paste0("mu[phantom(.)*phantom(.)*T] == ", mu_T),
    sigma_label = paste0("sigma[phantom(.)*phantom(.)*T] == ", sigma_T)
  )

df_LD <- df_LD |>
  mutate(
    mu_label = paste0("mu[phantom(.)*phantom(.)*T] == ", mu_T),
    sigma_label = paste0("sigma[phantom(.)*phantom(.)*T] == ", sigma_T)
  )

df_LD_sim <- df_LD_sim |>
  mutate(
    mu_label = paste0("mu[phantom(.)*phantom(.)*T] == ", mu_T),
    sigma_label = paste0("sigma[phantom(.)*phantom(.)*T] == ", sigma_T)
  )


# plot
colours_Dark2 <- RColorBrewer::brewer.pal(8, "Dark2")

LD_colours <- c(
  "0.1" = colours_Dark2[1],
  "0.2" = colours_Dark2[2],
  "0.4" = colours_Dark2[3],
  "0.6" = colours_Dark2[4],
  "0.8" = colours_Dark2[5],
  "0.9" = colours_Dark2[6]
)

plot_doseresponse <- function(rho_i) {
  
  ggplot(df_all) +
    geom_line(aes(x = c, y = p)) +
    geom_linerange(
      data = df_LD,
      aes(
        x = LD,
        ymin = 0,
        ymax = LD_prob,
        colour = factor(LD_prob)
      ),
      show.legend = TRUE
    ) +
    geom_point(
      data = df_LD_sim |> filter(rho == rho_i),
      aes(
        x = LD,
        y = D_b / N_b
      ),
      position = position_jitter(width = 0.15, height = 0),
      size = 0.9,
      alpha = 0.7
    ) +
    scale_colour_manual(
      name = "Dose",
      values = LD_colours,
      breaks = c("0.1", "0.2", "0.4", "0.6", "0.8", "0.9"),
      labels = c("LD10", "LD20", "LD40", "LD60", "LD80", "LD90")
    ) +
    scale_x_continuous(
      name = "Insecticide challenge [Dose]",
      breaks = c(0, 1, 5, 10, 50),
      labels = c("0", "1", "5", "10", "50"),
      limits = c(0, 50)
    ) +
    scale_y_continuous(
      name = "Mortality [Probability]",
      labels = scales::percent,
      limits = c(0, 1)
    ) +
    facet_grid(
      rows = vars(sigma_label),
      cols = vars(mu_label),
      labeller = label_parsed
    ) +
    coord_fixed(expand = TRUE) +
    theme(aspect.ratio = 2 / 3) +
    ggtitle(bquote(rho == .(rho_i)))
}

p_doseresponse_rho0 <- plot_doseresponse(0)

p_doseresponse_rho01 <- plot_doseresponse(0.1)


# For each setting run the SB only model 5 times, for control, LD10 and then adding one LD at a time
library(rstan)
library(tidybayes)
rstan_options(auto_write = TRUE)

LD_max <- c(20, 40, 60, 80, 90)

# compile once
stan_model <- rstan::stan_model("src/stan/nested_binomial_mechanistic_separated_p_b_control_sigmavest_newpriors_onlyBA_simstudy.stan")

# one fit per rho x LD_max
fit_settings <- df_LD_sim |>
  distinct(rho) |>
  tidyr::crossing(LD_max = LD_max)

fits <- vector("list", nrow(fit_settings))

for (j in seq_len(nrow(fit_settings))) {
  
  rho_j    <- fit_settings$rho[j]
  LD_max_j <- fit_settings$LD_max[j]
  
  df_i <- df_LD_sim |>
    filter(
      rho == rho_j,
      LD_label == "control" |
        LD_prob <= LD_max_j / 100
    ) |>
    mutate(
      treat_b = as.integer(LD_label != "control"),
      times_disc_dose = if_else(
        LD_label == "control",
        0,
        LD
      )
    )
  
  stan_data_i <- list(
    nT = n_distinct(df_i$group_number),
    nT_bint = n_distinct(df_i$group_number),
    
    sigma_SB = sigma_SB,
    p_b_control = p_b_control,
    
    S_b = nrow(df_i),
    T_b = df_i$group_number,
    treat_b = df_i$treat_b,
    times_disc_dose = df_i$times_disc_dose,
    N_b = df_i$N_b,
    D_b = df_i$D_b
  )
  
  fits[[j]] <- rstan::sampling(
    object = stan_model,
    data = stan_data_i,
    chains = 4,
    iter = 2000,
    seed = 123
  )
}

fit_settings <- fit_settings |>
  mutate(fit = fits)


# extract from stan, separately for the two rho values
fit_settings_rho0 <- fit_settings |>
  filter(rho == 0)

posterior_rho0 <- map2_dfr(
  fit_settings_rho0$fit,
  fit_settings_rho0$LD_max,
  \(fit_i, LD_max_i) {
    fit_i |>
      tidybayes::spread_draws(
        mu_T[group_number],
        sigma_T[group_number],
        ndraws = 2000,
        seed = 27
      ) |>
      mutate(LD_max = LD_max_i)
  }
)

  
fit_settings_rho01 <- fit_settings |>
  filter(rho == 0.1)

posterior_rho01 <- map2_dfr(
  fit_settings_rho01$fit,
  fit_settings_rho01$LD_max,
  \(fit_i, LD_max_i) {
    fit_i |>
      tidybayes::spread_draws(
        mu_T[group_number],
        sigma_T[group_number],
        ndraws = 2000,
        seed = 27
      ) |>
      mutate(LD_max = LD_max_i)
  }
)


## plotting
# dose grid on original dose scale
dose_grid <- seq(0, sqrt(500), length.out = 1000)^2

# true log-normal distributions
density_true <- param_combinations |>
  tidyr::crossing(Dose = dose_grid) |>
  mutate(
    density = dlnorm(
      Dose,
      meanlog = mu_T,
      sdlog = sigma_T
    )
  )


# posterior mixture densities
get_posterior_density <- function(posterior) {
  
  posterior |>
    select(group_number, LD_max, mu_T, sigma_T) |>
    tidyr::crossing(Dose = dose_grid) |>
    mutate(
      density_draw = dlnorm(
        Dose,
        meanlog = mu_T,
        sdlog = sigma_T
      )
    ) |>
    group_by(group_number, LD_max, Dose) |>
    summarise(
      density = mean(density_draw),
      .groups = "drop"
    )
}

density_rho0   <- get_posterior_density(posterior_rho0)
density_rho01 <- get_posterior_density(posterior_rho01)

# for facet grid
density_rho0 <- density_rho0 |>
  left_join(
    param_combinations,
    by = "group_number"
  )

density_rho01 <- density_rho01 |>
  left_join(
    param_combinations,
    by = "group_number"
  )

# add labels
density_rho0 <- density_rho0 |>
  mutate(
    mu_label = paste0("mu[phantom(.)*phantom(.)*T] == ", mu_T),
    sigma_label = paste0("sigma[phantom(.)*phantom(.)*T] == ", sigma_T)
  )

density_rho01 <- density_rho01 |>
  mutate(
    mu_label = paste0("mu[phantom(.)*phantom(.)*T] == ", mu_T),
    sigma_label = paste0("sigma[phantom(.)*phantom(.)*T] == ", sigma_T)
  )

density_true <- density_true |>
  mutate(
    mu_label = paste0("mu[phantom(.)*phantom(.)*T] == ", mu_T),
    sigma_label = paste0("sigma[phantom(.)*phantom(.)*T] == ", sigma_T)
  )

# plot function
plot_densities <- function(density_posterior, rho_label) {
  
  ggplot() +
    # true log-normal distribution
    geom_line(
      data = density_true,
      aes(
        x = Dose,
        y = density,
        linetype = "Given lethal dose distribution"
      ),
      colour = "black",
      linewidth = 1.5
    ) +
    # posterior mixture distributions
    geom_line(
      data = density_posterior,
      aes(
        x = Dose,
        y = density,
        colour = paste0("LD", LD_max)
      ),
      linewidth = 0.8,
      alpha = 0.6
    ) +
    scale_linetype_manual(
      name = NULL,
      values = c("Given lethal dose distribution" = "solid")
    ) +
    scale_colour_manual(
      name = "Estimated lethal dose distribution\nbased on intensity dose bioassays up to:",
      values = c(
        "LD20" = colours_Dark2[2],
        "LD40" = colours_Dark2[3],
        "LD60" = colours_Dark2[4],
        "LD80" = colours_Dark2[5],
        "LD90" = colours_Dark2[6]
      ),
      breaks = c("LD20", "LD40", "LD60", "LD80", "LD90")
    ) +
    facet_grid(
      rows = vars(sigma_label),
      cols = vars(mu_label),
      scales = "free_y",
      labeller = label_parsed
    ) +
    scale_x_sqrt(
      breaks = c(1, 5, 10, 20, 50, 100)
    ) +
    coord_cartesian(
      xlim = c(0, 50),
      ylim = c(0,2)
    ) +
    labs(
      x = "Lethal dose [Dose]",
      y = "Density",
      title = bquote(rho == .(rho_label))
    ) +
    theme(
      aspect.ratio = 2 / 3
    )
}

p_rho0 <- plot_densities(
  density_rho0,
  rho_label = 0
)

p_rho01 <- plot_densities(
  density_rho01,
  rho_label = 0.1
)



# compose
library(patchwork)

p_comb_rho01 <-
  (p_doseresponse_rho01 + labs(title = NULL)) /
  (p_rho01 + labs(title = NULL)) +
  plot_annotation(tag_levels = list(c("A", "B"))) 

# save plot
ggsave(
  "simulation_rho01.png",
  p_comb_rho01,
  width = 12,
  height = 12
)



