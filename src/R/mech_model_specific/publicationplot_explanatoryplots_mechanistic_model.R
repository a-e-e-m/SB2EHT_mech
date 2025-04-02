# generating publication plots 
# illustrative figures

library(tidyverse)
library(ellipse)
library(figpatch)
library(latex2exp)
library(scales)


# made up parameter values
  d <- c(1,5,10)
  sigma_v <- 0.2
  mu_x <- -1
  sigma_x <- 0.5
  mu_s <- seq(-10, 10, 0.05)
  sigma_s <- seq(0, 10, 0.05)
  
  mu_s_fix <- 0.8
  sigma_s_fix <- 1
  
  sample_sigma_s <- c(0.5, 1, 5)
  sample_mu_s <- c(log(0.5), log(2), log(7))

# dose-response curve 
  fun <- function(d, mu_s, sigma_s){
    pnorm( ( log(d) - mu_s ), 0, ( sqrt(sigma_v^2 + sigma_s^2) ))
  }

# data frames for plots
  adi <- as.data.frame(expand.grid(d = d, mu_s = mu_s, sigma_s = sigma_s)) |>
    mutate(p = fun(d, mu_s, sigma_s))

# ITN killing computation
  dfEHT <- tibble(killing = rep(0.1,10), prob = seq(0.1, 1, 0.1), x = 1 )
  EHT_killing <- 0 + ( 1 - 0) * pnorm( (mu_x - mu_s_fix) / sqrt(sigma_x^2 + sigma_s_fix^2 ))


# dose response plot
  times_disc_dose <- seq(0, 15, 0.05)
  p_b_control = 0
  
  # prepare data
  df_d <- tibble(mu_s = mu_s_fix, sigma_s = sigma_s_fix, c = times_disc_dose) |>
    mutate(
      p = p_b_control + ( 1 - p_b_control) * pnorm( (log(c) - mu_s ) / sqrt(sigma_v^2 + sigma_s^2 ))
    )
  
  p_doseresponse <- ggplot(df_d) +
    scale_fill_brewer(palette = "Dark2", name = "Dose") +
    scale_color_brewer(palette = "Dark2", name = "Dose") +
    geom_line(aes(x = c, y = p)) +
    geom_linerange(data = df_d |> filter(c %in% c(1,5,10)), aes(x = c, ymin = 0, ymax = p, colour = factor(c)), show.legend = FALSE) +
    scale_x_continuous(name = "Insecticide challenge [Dose]", breaks = c(0, 1, 5, 10), labels = c("0", "1", "5", "10"), limits = c(0,15)) +
    ylab("Killing [Probability]") + 
    scale_y_continuous(labels = percent, limits = c(0,1)) +
    coord_fixed(expand = TRUE) +
    theme(aspect.ratio=1)


# univariate exposure/tolerance plots  
  #plotting boundaries for 1d and 2d exposure-tolerance plots
  left_bounday = 0.05
  right_boundary = 50
  
  
  # 1d distribution plot
  df_1d <- expand_grid(x = seq(-2.5, 5, 0.01), d = d, colo = c("SB dose 1", "SB dose 5", "SB dose 10"), sigma_v = sigma_v, quantity = "exposure") |>
    mutate(
      density = dnorm(x, mean = log(d), sd = sigma_v)
    )
  
  scale_linetypes <- c("exposure" = "solid", "tolerance" = "dashed")
  scale_size <- c("Exposure" = 0.75, "Log lethal dose" = 0.8)
  
  p_1d <- ggplot() +
    scale_fill_brewer(palette = "Dark2", name = "Exposure in", labels = c("SB dose 1", "SB dose 5", "SB dose 10", "EHT")) +
    scale_color_brewer(palette = "Dark2", name = "Exposure in", labels = c("SB dose 1", "SB dose 5", "SB dose 10", "EHT")) +
    scale_size_manual(values = scale_size) +
    scale_linetype_manual(values = scale_linetypes, name = "Quantity") +
    geom_line(data = df_1d, aes(x = exp(x), y = density, colour = factor(d))) + 
    stat_function(data = data.frame(x = c(left_bounday, right_boundary)), aes(log(x), size = "Log lethal dose"), colour = "darkgrey", fun = rlang::as_function(~ dnorm(x = log(.x), mean = mu_s_fix, sd = sigma_s_fix)), n = 101, args = list(mean = mu_s_fix, sd = sigma_s_fix)) + 
    guides(size = guide_legend(title='')) +
    stat_function(data = data.frame(x = c(0.001, 100)), aes(log(x), colour = "EHT"), fun = rlang::as_function(~ dnorm(x = log(.x), mean = mu_x, sd = sigma_x)), n = 101, args = list(mean = mu_x, sd = sigma_x), alpha = 0.9) + 
    guides(colour = guide_legend(title='Exposure in')) +
    coord_fixed(expand = TRUE) +
    scale_x_continuous(transform = "log", breaks = c(0.1, 1, 5, 10), labels = c("0.1", "1", "5", "10"), limits = c(left_bounday, right_boundary)) +
    theme(aspect.ratio=1) +
    ylab("Density") + xlab("Exposure / Tolerance [Dose]")



# bivariate exposure-susceptibility plot
  df_XS <- expand_grid(d = d, sigma_v = sigma_v, mu_s = mu_s_fix, sigma_s = sigma_s_fix)
  df_XS_exp <- expand_grid(d = d, sigma_v = sigma_v, mu_s = mu_s_fix, sigma_s = sigma_s_fix)
  
  scale_fill_me <- c("Killing area" = "red")
  df_area <- tibble(x = c(exp(seq(log(left_bounday), log(right_boundary), 0.1)), right_boundary), y = c(exp(seq(log(left_bounday), log(right_boundary), 0.1)), left_bounday), g = "Killing area")
  
  p_XS <- ggplot() +
    scale_fill_manual(values = scale_fill_me) +
    scale_color_brewer(palette = "Dark2", name = "Exposure in", labels = c("SB dose 1", "SB dose 5", "SB dose 10", "EHT")) +
    # quick fix for legends
    geom_polygon(data= df_area, aes(x, y, fill = g), alpha = 0.2) +
    geom_abline(color='red') +
    geom_vline(xintercept = 1, color='darkgrey') +
    geom_hline(yintercept = 1, color='darkgrey') +
    coord_fixed(expand = TRUE) +
    scale_x_continuous(name = "Exposure [Dose]", transform = "log", breaks = c(0.1, 1, 5, 10), labels = c("0.1", "1", "5", "10"), limits = c(left_bounday, right_boundary)) +
    scale_y_continuous(name = "Tolerance [Dose]", transform = "log", breaks = c(0.1, 1, 5, 10), labels = c("0.1", "1", "5", "10"), limits = c(left_bounday, right_boundary)) +
    ggtitle("exposure-susceptibility on log-scale") + 
    labs(colour = "Dose", alpha = "Killing probability", fill = "")


  #IDB exposure contours
  library(plyr)  ## not necessary, but convenient
  for (k in 1:nrow(df_XS)){
    m <- c(log(df_XS$d[k]), df_XS$mu_s[k])
    alpha_levels <- seq(0.1,0.9,by=0.2) ## or whatever you want
    names(alpha_levels) <- alpha_levels ## to get id column in result
    contour_data <- ldply(alpha_levels,ellipse,x=0,
                          scale=c(df_XS$sigma_v[k], df_XS$sigma_s[k]),  ## needed for positional matching
                          centre=m) |>
      mutate(d = df_XS$d[k],
             x = exp(x),
             y = exp(y))
    p_XS <- p_XS +
    geom_path(data = contour_data, aes(x = x, y = y, group=.id, colour = factor(d)), show.legend = F)
  }

  #EHT exposure contours
  m <- c(mu_x, df_XS$mu_s[k])
  alpha_levels <- seq(0.1,0.9,by=0.2) ## or whatever you want
  names(alpha_levels) <- alpha_levels ## to get id column in result
  contour_data <- ldply(alpha_levels,ellipse,x=0,
                        scale=c(sigma_x, sigma_s_fix),  ## needed for positional matching
                        centre=m) |>
    mutate(x = exp(x),
           y = exp(y))
  p_XS <- p_XS +
    geom_path(data = contour_data, aes(x = x, y = y, group=.id, colour = "EHT"), show.legend = F)

# SI plots
  #  Parameter identifiability plot (now in SI)
  c1 <- emdbook::curve3d(fun(1,x,y), xlim=c(-10,10), ylim=c(0,10), sys3d="none", n = c(100, 100))
  dimnames(c1$z) <- list(c1$x,c1$y)
  m1 <- reshape2::melt(c1$z) |> mutate(d = 1)
  c5 <- emdbook::curve3d(fun(5,x,y), xlim=c(-10,10), ylim=c(0,10), sys3d="none", n = c(100, 100))
  dimnames(c5$z) <- list(c5$x,c5$y)
  m5 <- reshape2::melt(c5$z) |> mutate(d = 5)
  c10 <- emdbook::curve3d(fun(10,x,y), xlim=c(-10,10), ylim=c(0,10), sys3d="none", n = c(100, 100))
  dimnames(c10$z) <- list(c10$x,c10$y)
  m10 <- reshape2::melt(c10$z) |> mutate(d = 10)
  
  mm <- rbind(m1, m5, m10)
  
  p_contourBA <- ggplot(adi, colour = "black") +
    geom_point(y = sigma_s_fix, x = mu_s_fix, colour = "black", size = 4, shape = 1) +
    scale_color_brewer(palette = "Dark2", name = "Dose in SB") +
    geom_contour(data=m1,
                 aes(x=exp(Var1),y=Var2,z=value, colour = factor(d)),
                 breaks = df_d |> filter(c == 1) |> pull(p), show.legend = FALSE
    ) +
    geom_contour(data=m5,
                 aes(x=exp(Var1),y=Var2,z=value, colour = factor(d)),
                 breaks = df_d |> filter(c == 5) |> pull(p), show.legend = FALSE
    ) +
    geom_contour(data=m10,
                 aes(x=exp(Var1),y=Var2,z=value, colour = factor(d)),
                 breaks = df_d |> filter(c == 10) |> pull(p), show.legend = FALSE
    ) +
    ylim(0, 10) +
    coord_fixed(expand = TRUE) +
    xlab(TeX(r"($LD_{50}$ (exp($\mu_{\ T}$)) \[Dose\] )")) +
    ylab(TeX(r"(\overset{Heterogeneity in tolerance ($\sigma_{\ T}$) }{\[log-Dose\]})")) +
    scale_x_continuous(transform = "log", breaks = c(0.01, 1, 10, 100), labels = c("0.01", "1", "10", "100"), limits = c(exp(-10),exp(10))) +
    ggtitle("resistance metric") +
    labs(colour = "Dose in SB") +
    theme(aspect.ratio=1) 
  
  
  
  
  # ITN killing effect across resistance metrics plot (now in SI)
  df_2d <-  tibble(expand_grid(mu_s = mu_s, sigma_s = sigma_s), p_h_control = 0, mu_x = mu_x, sigma_x = sigma_x) |>
    mutate(
      prob_D_h = p_h_control + ( 1 - p_h_control) * pnorm( (mu_x - mu_s) / sqrt(sigma_x^2 + sigma_s^2 ))
    )
  
  df_2d_exp <-  tibble(expand_grid(mu_s = exp(mu_s), sigma_s = sigma_s), p_h_control = 0, mu_x = mu_x, sigma_x = sigma_x) |>
    mutate(
      prob_D_h = p_h_control + ( 1 - p_h_control) * pnorm( (mu_x - log(mu_s)) / sqrt(sigma_x^2 + sigma_s^2 ))
    )
  
  
  p_contourEHT <- ggplot(df_2d_exp, aes(x = mu_s, y = sigma_s)) +
    geom_contour_filled(aes(z = prob_D_h)) +
    geom_point(y = sigma_s_fix, x = mu_s_fix, colour = "white", size = 4, shape = 16) +
    ylim(0, 10) +
    coord_fixed(expand = TRUE) +
    xlab(TeX(r"($LD_{50}$ (exp($\mu_{\ T}$)) \[Dose\] )")) +
    ylab(TeX(r"(\overset{Heterogeneity in tolerance ($\sigma_{\ T}$) }{\[log-Dose\]})")) +
    scale_x_continuous(transform = "log", breaks = c(0.01, 1, 10, 100), labels = c("0.01", "1", "10", "100"), limits = c(exp(-10), exp(10))) +
    labs(fill = "ITN killing effect [Probability]") +
    theme(aspect.ratio=1)
  
  
  # ITN killing prediction plot (now in SI)
  p_EHT <- ggplot(dfEHT, aes(x = x, y = killing)) +
    scale_fill_viridis_c() +
    geom_col(aes(fill = prob), width = 0.1) +
    geom_label(x = 1, y = EHT_killing, label = paste0(round(EHT_killing*100,1), "%"), size = 3) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none",
      aspect.ratio = 4
    ) +
    ylab("ITN killing effect [Probability]") +
    coord_cartesian(expand = TRUE) +
    scale_y_continuous(labels = percent)
  

  
# composing plots
library(patchwork)

# data source plots
  #load cartoons
  img_1 <- fig(file.path("plots_mechmodel_article", "cartoon_1_new_title.png"))
  img_2 <- fig(file.path("plots_mechmodel_article", "cartoon_2_new_title.png"))
  img_3 <- fig(file.path("plots_mechmodel_article", "cartoon_3_new_title.png"))
  
  # data source only
  p_cartoons <- wrap_plots( img_1, img_2, img_3 )  +
    plot_annotation(tag_levels = list(c("A.1", "A.2", "B"))) +
    plot_layout(guides = "collect") &
    theme(legend.position = 'bottom') &
    labs(title = NULL)
  
  ggsave(file.path("plots_mechmodel_article","cartoons_comb_title.png"), p_cartoons,  height = 3.5, width = 10) 



# model identifiability plot

p_model_expl_identify <- wrap_plots(p_doseresponse + theme(aspect.ratio=1), p_contourBA + theme(aspect.ratio=1), p_contourEHT + theme(aspect.ratio=1), p_EHT,
                                    nrow = 1 )  +
  plot_annotation(tag_levels = list(c("A", "B", "C", "D"))) +
  plot_layout(guides = "collect") &
  theme(legend.position = 'bottom') &
  labs(title = NULL) &
  guides(fill = "none")


ggsave(file.path("plots_mechmodel_article_SI","explanatory_plot_model_identify.png"), p_model_expl_identify,  height = 4.5, width = 12) 
