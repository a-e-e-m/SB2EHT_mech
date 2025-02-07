# generating publication plots 
# resistance metric and extrapolation of ITN killing effect for possible dose-response patterns 

# index of model run to fetch sigma_v, mu_x and sigma_x estimates from
i <- 4
# index of corresponding model run with ind dose BA data only
i_BAonly <- 3


# config and prepare
###############################
#libraries
library("rstan")
library("patchwork")
library("ggrepel")
library("scales")
library("viridisLite")
library("latex2exp")
library("tidyverse")
library("tidybayes")


# load runs
models2run <- read.csv(file = './models2run.csv', sep = ';', stringsAsFactors = FALSE)

# folders
dir_i <- paste0("fitting/", models2run$run[i]) # folder
dir_i_BAonly <- paste0("fitting/", models2run$run[i_BAonly]) # folder

## grouping

  # group_number lookup joint model
  B_all <- readRDS(file.path(dir_i, 'B_all.rds'))
  
  group_number_lookup <- B_all |>
    filter(treat == 1) |>
    group_by(group_number, country, site, year, int_dose_available, Trial_code) |>
    dplyr::summarise(insecticide = paste(unique(insecticide), collapse = ', '),
              test_type = paste(unique(test_type), collapse = ', '),
              .groups = "keep"
    )
  
  # group_number lookup BAonly model
  B_all_BAonly <- readRDS(file.path(paste0("fitting/", models2run$run[i_BAonly]) , 'B_all.rds'))
  group_number_conversion <- B_all_BAonly |> 
    filter(treat == 1) |>
    group_by(group_number, country, site, year, int_dose_available, Trial_code) |>
    dplyr::summarise(insecticide = paste(unique(insecticide), collapse = ', '),
              test_type = paste(unique(test_type), collapse = ', '),
              .groups = "keep"
    ) |>
    dplyr::rename(group_number_conversion = group_number)
  
  # group number conversion
  group_number_lookup <- group_number_lookup |>
    dplyr::left_join(group_number_conversion)
  
  # reorder groups
  # first the INT dose groups, then other groups by country, year, insecticide (with this ordering)
  group_number_lookup <- group_number_lookup |> 
    arrange(group_number_conversion, country, site, year, insecticide) |>
    rowid_to_column(var = "group_number_pub") 
  

## load data
# load data joint model (BA and EHT)
fit_i <- readRDS(file.path(dir_i, 'fit.rds'))


# parameter draws and lethal dose computation
# joint model
draws_i_spread <- spread_draws(fit_i, p_b_control, sigma_v, mu_d[group_number], sigma_d[group_number], p_h_control[group_number], mu_x, sigma_x, ndraws = NULL, seed = 27) |>
  dplyr::left_join(group_number_lookup) |> 
  dplyr::rename("Group" = "group_number_pub") |>
  ungroup() |>
  select(!c("group_number", "group_number_conversion")) |>
  group_by(Group) |>
  mutate(Group = factor(Group)
  ) 


# simulated dose-response patterns
mu_s_sim <- seq(-20, 20, 0.1) 
sigma_s_sim <- seq(0,25, 0.1)
p_b_control <- 0
p_h_control <- 0

mu_s_fix <- c(lower = -5, upper = 5)
sigma_s_fix <- c( lower = 0.1, upper = 5)

# exposure parameter draws
df_fit <- draws_i_spread |>
  ungroup() |>
  dplyr::filter(Group == 1) |>
  select(sigma_v, mu_x, sigma_x, .draw)

# fixing mu_s
df_sim_musfix <- expand_grid(df_fit, tibble(mu_s = mu_s_fix, scenario = c("-5", "5")), sigma_s = sigma_s_sim, p_b_control, p_h_control) |>
  dplyr::mutate(
    prob_D_h = p_h_control + ( 1 - p_h_control) * pnorm( (mu_x - mu_s) / sqrt(sigma_x^2 + sigma_s^2 ))
  )

df_sim_musfix_sum <- df_sim_musfix |>
  group_by(scenario, mu_s, sigma_s, p_b_control, p_h_control) |>
  dplyr::summarise(
    prob_D_h_l = quantile(prob_D_h, 0.025),
    prob_D_h_median = median(prob_D_h),
    prob_D_h_mean = mean(prob_D_h),
    prob_D_h_r = quantile(prob_D_h, 0.975),
    .groups = "keep"
            )

# fixing sigma_s
df_sim_sigmasfix <- expand_grid(df_fit, mu_s = mu_s_sim, tibble(sigma_s = sigma_s_fix, scenario = c("0.1", "5")), p_b_control, p_h_control) |>
  mutate(
    prob_D_h = p_h_control + ( 1 - p_h_control) * pnorm( (mu_x - mu_s) / sqrt(sigma_x^2 + sigma_s^2 ))
  )

df_sim_sigmasfix_sum <- df_sim_sigmasfix |>
  group_by(scenario, mu_s, sigma_s, p_b_control, p_h_control) |>
  dplyr::summarise(
    prob_D_h_l = quantile(prob_D_h, 0.025),
    prob_D_h_median = median(prob_D_h),
    prob_D_h_mean = mean(prob_D_h),
    prob_D_h_r = quantile(prob_D_h, 0.975),
    .groups = "keep"
  )


## susceptibility parameters and EHT exposure parameters
params_B <- draws_i_spread |>
  mutate(LD50 = exp(mu_d),
         EHT_killing_effect = 0 + ( 1 - 0) * pnorm( (mu_x - mu_d) / sqrt(sigma_x^2 + sigma_d^2 ))
         ) 

params_B_summary <- params_B |>
  group_by(Group, country, site, year, insecticide, int_dose_available) |>
  dplyr::summarise(mu_d_median = median(mu_d),
            mu_d_mean = mean(mu_d),
            sigma_d_median = median(sigma_d),
            sigma_d_mean = mean(sigma_d),
            EHT_killing_effect_mean = mean(EHT_killing_effect),
            EHT_killing_effect_median = median(EHT_killing_effect),
            EHT_killing_effect_q025 = quantile(EHT_killing_effect, probs = 0.025),
            EHT_killing_effect_q975 = quantile(EHT_killing_effect, probs = 0.975),
            .groups = "keep")

pred <- params_B |>
  median_qi(EHT_killing_effect, LD50, sigma_d, .width = 0.95) |>
  dplyr::left_join(summarise(params_B, corr = cor(mu_d, sigma_d, method = "pearson"))) |>
  select(!c(".width", ".point", ".interval")) 

#function for plotting HDI contours
get_hpd2d_level <- function(x, y, prob, ...) {
  post1 <- MASS::kde2d(x, y)
  dx <- diff(post1$x[1:2])
  dy <- diff(post1$y[1:2])
  sz <- sort(post1$z)
  c1 <- cumsum(sz) * dx * dy
  levels <- sapply(prob, function(x) {
    approx(c1, sz, xout = 1 - x, ties = mean)$y
  })
  return(levels)
}


# heat map wrt either mean or median EHT exposure parameters
df_2d <-  tibble(expand_grid(mu_s = mu_s_sim, sigma_s = sigma_s_sim), p_h_control = 0, mu_x_median = median(df_fit$mu_x), sigma_x_median = median(df_fit$sigma_x), mu_x_mean = mean(df_fit$mu_x), sigma_x_mean = mean(df_fit$sigma_x)) |>
  mutate(
    prob_D_h_median = p_h_control + ( 1 - p_h_control) * pnorm( (mu_x_median - mu_s) / sqrt(sigma_x_median^2 + sigma_s^2 )),
    prob_D_h_mean = p_h_control + ( 1 - p_h_control) * pnorm( (mu_x_mean - mu_s) / sqrt(sigma_x_mean^2 + sigma_s^2 ))
  )

pp <- ggplot(df_2d, aes(x = exp(mu_s), y = sigma_s)) +
  geom_contour_filled(aes(z = prob_D_h_median), show.legend = FALSE) +
  geom_point(data = params_B_summary |> filter(as.numeric(Group) <= 7), aes(x = exp(mu_d_median), y = sigma_d_median, group = Group), colour = "white") +
  geom_vline(aes(xintercept = exp(median(df_fit$mu_x))), colour = "red") +
  geom_hline(aes(yintercept = sigma_s_fix["lower"]), linetype = "dashed", colour = 'black', linewidth = 0.75, show.legend = FALSE) +
  geom_hline(aes(yintercept = sigma_s_fix["upper"]), linetype = "dotted", colour = 'black', linewidth = 0.75, show.legend = FALSE) +
  geom_vline(aes(xintercept = exp(mu_s_fix["lower"])), linetype = "dotdash", colour = 'black', linewidth = 0.75, show.legend = FALSE) +
  geom_vline(aes(xintercept = exp(mu_s_fix["upper"])), linetype = "longdash", colour = 'black', linewidth = 0.75, show.legend = FALSE) +
  coord_fixed(expand = TRUE, ratio = (8 - (-6)) / 6 ) +
  xlab(TeX(r"($LD_{50}$ (exp($\mu_{\ T}$)) \[Dose\] )")) +
  ylab(TeX(r"(\overset{Heterogeneity in tolerance ($\sigma_{\ T}$) }{\[log-Dose\]})")) +
  ylim(0, 6) +
  scale_x_continuous(transform = "log", breaks = c(0.01, 1, 10, 100), labels = c("0.01", "1", "10", "100"), limits = c(exp(-6), exp(8))) +
  labs(fill = "ITN killing effect [probability]")
   

for (kk in 1:7){
  L <- with(params_B |> filter(as.numeric(Group) == kk), get_hpd2d_level(mu_d, sigma_d, prob = 0.95))
  pp <- pp +
    geom_density_2d(data = params_B |> filter(as.numeric(Group) == kk), aes(x = exp(mu_d), y = sigma_d), breaks = L, contour_var = "density", colour = "white", adjust = 5)
}

pp <- pp + geom_label_repel(data = params_B_summary |> filter(as.numeric(Group) <= 7), aes(x = exp(mu_d_median), y = sigma_d_median, group = Group, label = Group), fill = "white", segment.color = "white", size = 2.5)


# EHT killing pred
df_scale <- tibble(killing = rep(0.1,10), prob = seq(0.1, 1, 0.1), x = 4 )

# point_range with colour scale in the back (0-100%)
  p_EHT <- ggplot(df_scale, aes(x = x, y = killing)) +
    scale_fill_viridis_c() +
    geom_col(aes(fill = prob), width = 8, show.legend = FALSE) +
    geom_pointrange(data = params_B_summary |> filter(as.numeric(Group) <= 7), aes(x = as.numeric(Group), y = EHT_killing_effect_median, ymin = EHT_killing_effect_q025, ymax = EHT_killing_effect_q975), colour = "white") +
    geom_label(data = params_B_summary |>  filter(as.numeric(Group) <= 7), aes(x = as.numeric(Group), y = EHT_killing_effect_median, label = as.character(Group)), size = 2.5) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none",
      aspect.ratio = 4
    ) +
    coord_cartesian(expand = TRUE) +
    scale_x_continuous(name = "Assay pair") +
    scale_y_continuous(name = "ITN killing effect [Probability]", labels = percent) +
    theme(legend.position = "none")
  
  ggsave(file.path("plots_mechmodel_article", "ITNest.png"), p_EHT, width = 2, height = 6)
  

## plot over variable sigma_s
  param_snapshots <- c(1, 2, 4, 8)
  param_snapshots_df <- tibble(param_snapshots)
  linetypes_mu_s_fix <- c("-5" = "dotdash", "5" = "longdash")
  p_mu_s_fix <- ggplot(data = df_sim_musfix_sum, aes(x = sigma_s)) +
    scale_linetype_manual(values = linetypes_mu_s_fix, name = TeX(r"($\mu_{\ T}$ fixed to)")) +
    geom_line(aes(y = prob_D_h_median, linetype = scenario), linewidth = 0.75, show.legend = FALSE) +
    geom_point(data = param_snapshots_df, aes(x = param_snapshots, y=-Inf), color="blue") +
    ylab("ITN killing \neffect [Prob]") + xlab(TeX(r"($\sigma_{\ T}$ \[log-Dose\])")) +
    scale_y_continuous(labels = percent) +
    ggtitle("Experimental hut trial") + 
    coord_cartesian(clip = 'off', expand = TRUE) +
    xlim(0,10) +
    theme(axis.title = element_text(size = 9))
  
  # doseresponse snapshot plot
  c <- seq(0, 15, 0.05)
  p_mu_s_fix_dose <- ggplot() +
    scale_linetype_manual(values = linetypes_mu_s_fix, name = TeX(r"($\mu_{\ T}$ fixed to)"))
  
  for (j in seq_along(param_snapshots)){
    df_sim_j <- expand_grid(df_fit, tibble(mu_s = mu_s_fix, scenario = c("-5", "5")), sigma_s = param_snapshots[j], p_b_control, p_h_control, c) |>
      mutate(
        p = p_b_control + ( 1 - p_b_control) * pnorm( (log(c) - mu_s ) / sqrt(sigma_v^2 + sigma_s^2 ))
      )
    
    df_sim_j_sum <- df_sim_j |> 
      group_by(scenario, mu_s, sigma_s, c) |>
      summarise(
        p_l = quantile(p, 0.025),
        p_median = median(p),
        p_r = quantile(p, 0.975),
        .groups = "keep"
        )
    if (j == 1){
      p_mu_s_fix_dose <- p_mu_s_fix_dose +
        geom_line(data = df_sim_j_sum, aes(x = c, y = p_median, linetype = scenario), linewidth = 0.75) 
    }else{
      p_mu_s_fix_dose <- p_mu_s_fix_dose +
        geom_line(data = df_sim_j_sum, aes(x = c, y = p_median, linetype = scenario), linewidth = 0.75, show.legend = FALSE) 
    }
  }
  
  p_mu_s_fix_dose <- p_mu_s_fix_dose +
    facet_grid(cols = vars(sigma_s), labeller = "label_both") +
    theme(
      strip.text.x = element_blank(),
      panel.border = element_rect(colour = "blue", size=1, fill=NA)
    ) +
    ylab("Killing in\nSB [Prob]") + xlab("Insecticide challenge in ID-SB [Dose]") +
    scale_y_continuous(labels = percent) +
    ggtitle("dose-response in bio assay") + 
    labs(fill = bquote(mu[s]), colour = bquote(mu[s])) +
    theme(axis.title = element_text(size = 9))
  
  p_mu_s_comb <- p_mu_s_fix / p_mu_s_fix_dose +
    plot_layout(guides = "collect") &
    theme(legend.position = 'bottom')
  
  # ggsave(file.path("plots_mechmodel_article", "sim_mus_comb.png"), p_mu_s_comb, width = 10, height = 6)
  # ggsave(file.path("plots_mechmodel_article", "sim_mus.png"), p_mu_s_fix, width = 10, height = 6)
  # ggsave(file.path("plots_mechmodel_article", "sim_mus_dose.png"), p_mu_s_fix_dose, width = 10, height = 6)



## plot over variable mu_s
  param_snapshots <- c(-5, 0, log(10), 5)
  param_snapshots_df <- tibble(param_snapshots)
  linetypes_sigma_s_fix <- c("0.1" = "dashed", "5" = "dotted")
  p_sigma_s_fix <- ggplot(data = df_sim_sigmasfix_sum, aes(x = exp(mu_s))) +
    scale_linetype_manual(values = linetypes_sigma_s_fix, name = TeX(r"($\sigma_{\ T}$ fixed to)")) +
    geom_vline(aes(xintercept = exp(median(df_fit$mu_x))), colour = "red") +
    geom_line(aes(y = prob_D_h_median, linetype = scenario), linewidth = 0.75, show.legend = FALSE) +
    geom_point(data = param_snapshots_df, aes(x = exp(param_snapshots), y=-Inf), color="blue") +
    ylab("ITN killing \neffect [Prob]") + 
    xlab(TeX(r"($LD_{50}$ (exp($\mu_{\ T}$)) \[Dose\] )")) +
    scale_y_continuous(labels = percent) +
    ggtitle("Experimental hut trial") + 
    coord_cartesian(clip = 'off', expand = TRUE) +
    scale_x_continuous(transform = "log", breaks = c(0.01, 1, 10, 100), labels = c("0.01", "1", "10", "100"), limits = c(exp(-6), exp(8))) +
    theme(axis.title = element_text(size = 9))
  
  # doseresponse snapshot plot
  c <- seq(0, 15, 0.05)
  p_sigma_s_fix_dose <- ggplot() +
    scale_linetype_manual(values = linetypes_sigma_s_fix, name = TeX(r"($\sigma_{\ T}$ fixed to)"))
  
  for (j in seq_along(param_snapshots)){
    df_sim_j <- expand_grid(df_fit, tibble(sigma_s = sigma_s_fix, scenario = c("0.1", "5")), mu_s = param_snapshots[j], p_b_control, p_h_control, c) |>
      mutate(
        p = p_b_control + ( 1 - p_b_control) * pnorm( (log(c) - mu_s ) / sqrt(sigma_v^2 + sigma_s^2 ))
      )
    
    df_sim_j_sum <- df_sim_j |> 
      group_by(scenario, mu_s, sigma_s, c) |>
      summarise(
        p_l = quantile(p, 0.025),
        p_median = median(p),
        p_r = quantile(p, 0.975),
        .groups = "keep"
      )
    if (j == 1){
      p_sigma_s_fix_dose <- p_sigma_s_fix_dose +
        geom_line(data = df_sim_j_sum, aes(x = c, y = p_median, linetype = scenario), linewidth = 0.75) 
    }else{
      p_sigma_s_fix_dose <- p_sigma_s_fix_dose +
        geom_line(data = df_sim_j_sum, aes(x = c, y = p_median, linetype = scenario), linewidth = 0.75, show.legend = FALSE) 
    }
  }
  
  p_sigma_s_fix_dose <- p_sigma_s_fix_dose +
    facet_grid(cols = vars(mu_s), labeller = "label_both") +
    theme(
      strip.text.x = element_blank(),
      panel.border = element_rect(colour = "blue", size=1, fill=NA)
    ) +
    ylab("Killing in\nSB [Prob]") + xlab("Insecticide challenge in ID-SB [Dose]") +
    scale_y_continuous(labels = percent) +
    ggtitle("dose-response in bio assay") + 
    labs(fill = bquote(sigma[s]), colour = bquote(sigma[s])) +
    theme(axis.title = element_text(size = 9))
  
  p_sigma_s_comb <- p_sigma_s_fix / p_sigma_s_fix_dose +
    plot_annotation(tag_levels = 'A') +
    plot_layout(guides = "collect") &
    theme(legend.position = 'bottom') &
    labs(title = NULL)
  
  # ggsave(file.path("plots_mechmodel_article", "sim_sigmas_comb.png"), p_sigma_s_comb, width = 10, height = 6)
  # ggsave(file.path("plots_mechmodel_article", "sim_sigmas.png"), p_sigma_s_fix, width = 10, height = 6)
  # ggsave(file.path("plots_mechmodel_article", "sim_sigmas_dose.png"), p_sigma_s_fix_dose, width = 10, height = 6)


# best layout so far
  patch_LHS <- pp + labs(tag = 'A') + ( p_EHT + labs(tag = 'B')) & theme(legend.position = "none")
  patch_RHS <-  (p_sigma_s_fix + labs(tag = 'C.1')) / ( p_sigma_s_fix_dose + labs(tag = 'C.2') )  /  ( p_mu_s_fix + labs(tag = 'D.1')) / ( p_mu_s_fix_dose + labs(tag = 'D.2') )
  
  p_allcomb <- patch_LHS + patch_RHS +
    plot_layout(guides = "collect") &
    theme(legend.position = 'bottom', legend.title.position = "top", legend.byrow = TRUE) &
    labs(title = NULL) &
    theme(legend.key.width = unit(2, "cm"))

  ggsave(file.path("plots_mechmodel_article", "IRcharANDITNestANDITNpredANDdoseresp.png"), p_allcomb, width = 12, height = 6)  
  
  
  
  # SI plot
  candidate_groups_small <- pred |> filter(abs(corr) < 0.5) |> pull(Group)
  candidate_groups_small <- c(8, 17, 22, 24, 26, 27, 35)
  candidate_groups_small <- c(23)
  candidate_groups <- seq(1,35)
  
  pp_SI <- ggplot(df_2d, aes(x = exp(mu_s), y = sigma_s)) +
    #scale_fill_viridis_c() +
    #geom_raster(aes(fill=prob_D_h)) + 
    geom_contour_filled(aes(z = prob_D_h_median), show.legend = FALSE) +
    geom_point(data = params_B_summary, aes(x = exp(mu_d_median), y = sigma_d_median, group = Group), colour = "white") +
    geom_vline(aes(xintercept = exp(median(df_fit$mu_x))), colour = "red") +
    geom_label_repel(data = params_B_summary, aes(x = exp(mu_d_median), y = sigma_d_median, label = as.character(Group)), fill = "white", segment.color = "white", size = 2, label.padding = 0.15,  max.overlaps = 15) +
    ylim(0, 11) +
    xlab(TeX(r"($LD_{50}$ (exp($\mu_{\ T}$)) \[Dose\] )")) +
    ylab(TeX(r"(\overset{Heterogeneity in tolerance ($\sigma_{\ T}$) }{\[log-Dose\]})")) +
    scale_x_continuous(transform = "log", breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000), labels = c("0.001", "0.01", "0.1", "1", "10", "100", "1000", "10,000"), limits = c(exp(-10), exp(13))) +
    labs(fill = "ITN killing effect [probability]") +
    coord_fixed(expand = TRUE, ratio = 1.4)
  
  
  
  # EHT killing pred
  df_scale <- tibble(killing = rep(0.1,10), prob = seq(0.1, 1, 0.1), x = 17 )
  
  # point_range with colour scale in the back (0-100%)
  p_EHT <- ggplot(df_scale, aes(x = x, y = killing)) +
    scale_fill_viridis_c() +
    geom_col(aes(fill = prob), width = 35, show.legend = FALSE) +
    geom_pointrange(data = params_B_summary, aes(x = as.numeric(Group)-1, y = EHT_killing_effect_median, ymin = EHT_killing_effect_q025, ymax = EHT_killing_effect_q975), colour = "white") +
    geom_label(data = params_B_summary, aes(x = as.numeric(Group)-1, y = EHT_killing_effect_median, label = as.character(Group)), size = 2.0) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none",
      aspect.ratio = 0.65
    ) +
    xlab("Assay pair") +
    ylab("ITN killing effect [Probability]") +
    coord_cartesian(expand = TRUE) +
    scale_y_continuous(labels = percent) +
    theme(legend.position = "none")
  
  #ggsave(file.path("plots_mechmodel_article_SI", "ITNest_SI.png"), p_EHT, width = 8, height = 8)
  
  
  # combined plot
  p_SI_comb <- wrap_plots( pp_SI, p_EHT)  +
    plot_annotation(tag_levels = list(c("A", "B"))) +
    plot_layout(guides = "collect") &
    theme(legend.position = 'bottom') &
    labs(title = NULL) &
    guides(fill = "none")
  
  ggsave(file.path("plots_mechmodel_article_SI", "IRcharANDITNest_SI.png"), p_SI_comb, width = 12, height = 4.5)
  
  
  
  
 


