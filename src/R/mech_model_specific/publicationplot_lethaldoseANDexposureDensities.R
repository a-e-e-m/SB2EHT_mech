# generating publication plots 
# tolerance and exposure 


# config and prepare
###############################
###############################
# packages
library("tidyverse")
library("tidybayes")
library("ggdist")
library("rstan")
library("patchwork")
library("ggnewscale")
library("RColorBrewer")

# load runs
models2run <- read.csv(file = './models2run.csv', sep = ';', stringsAsFactors = FALSE)

# index of model run to fetch sigma_v, mu_x and sigma_x estimates from
i <- 4
dir_i <- paste0("fitting/", models2run$run[i]) # folder

# index of corresponding model run with int dose BA data only
i_BAonly <- 3
dir_i_BAonly <- paste0("fitting/", models2run$run[i_BAonly]) # folder

## grouping

  # group_number lookup joint model
  B_all <- readRDS(file.path(dir_i, 'B_all.rds'))
  
  group_number_lookup <- B_all |>
    filter(treat == 1) |>
    group_by(group_number, country, site, year, int_dose_available, Trial_code) |>
    summarise(insecticide = paste(unique(insecticide), collapse = ', '),
              test_type = paste(unique(test_type), collapse = ', '),
              .groups = "keep"
    )
  
  # group_number lookup BAonly model
  B_all_BAonly <- readRDS(file.path(paste0("fitting/", models2run$run[i_BAonly]) , 'B_all.rds'))
  group_number_conversion <- B_all_BAonly |> 
    filter(treat == 1) |>
    group_by(group_number, country, site, year, int_dose_available, Trial_code) |>
    summarise(insecticide = paste(unique(insecticide), collapse = ', '),
              test_type = paste(unique(test_type), collapse = ', '),
              .groups = "keep"
    ) |>
    dplyr::rename(group_number_conversion = group_number)
  
  # group number conversion
  group_number_lookup <- group_number_lookup |>
    left_join(group_number_conversion)
  
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
draws_i_spread <- spread_draws(fit_i, p_b_control, sigma_v, mu_d[group_number], sigma_d[group_number], p_h_control[group_number], mu_x, sigma_x, ndraws = 2000, seed = 27) |>
  left_join(group_number_lookup) |> 
  rename("Group" = "group_number_pub") |>
  ungroup() |>
  select(!c("group_number", "group_number_conversion")) |>
  group_by(Group) |>
  mutate(Group = factor(Group),
         LD50 = exp(mu_d),
         lethal_dose = rlnorm(n(), mu_d, sigma_d),
         log_lethal_dose = rnorm(n(), mu_d, sigma_d),
         exposure_SB_1 = rlnorm(n(), 0, sigma_v),
         exposure_SB_5 = rlnorm(n(), log(5), sigma_v),
         exposure_SB_10 = rlnorm(n(), log(10), sigma_v),
         exposure_EHT = rlnorm(n(), mu_x, sigma_x)
  ) 

pred_i <- draws_i_spread |>
  median_qi(LD50, sigma_d, .width = 0.95) |>
  select(!c(".width", ".point", ".interval")) 


# exposure to long format with additional column for doses and EHT
draws_i_spread_exp <- pivot_longer(draws_i_spread |> filter(Group ==1) |> ungroup() |> select(!c(mu_d, sigma_d)), cols = starts_with("exposure"), names_to = "exposure_type", names_prefix = "exposure_", values_to = "exposure_value") |>
  mutate(exposure_type = case_when(exposure_type == "SB_1" ~ "SB dose 1",
                                   exposure_type == "SB_5" ~ "SB dose 5",
                                   exposure_type == "SB_10" ~ "SB dose 10",
                                   exposure_type == "SB_15" ~ "SB dose 15",
                                   exposure_type == "EHT" ~ "EHT",
                                   TRUE ~ NA)) 

pred_i_exp <- draws_i_spread_exp |>
  group_by(exposure_type) |>
  median_qi(exposure_value, .width = 0.95) |>
  select(!c(".width", ".point", ".interval")) 


# manual colours from the Dark2 palette
colours_Dark2 <- brewer.pal(8, "Dark2")

# plotting lethal dose, LD50 and variance
# lethal dose plot
p_lethal_dose <- ggplot(data = draws_i_spread |> filter(Group %in% c(1, 2, 3, 4, 5, 6, 7)), aes(x = lethal_dose, y = Group)) +
  scale_fill_identity(guide = "legend", name = "Lethal exposure distribution", labels = c("")) +
  # NOTE: intervals are with respect to 0.5 and 0.95 probability mass
  stat_halfeye(.width = c(0.5, 0.95), aes(fill = "darkgrey"), p_limits = c(NA, 0.999), trim = TRUE, expand = TRUE, show.legend = FALSE) +
  new_scale_fill() +
  scale_color_manual(name = "Doses tested in ID-SB", values = colours_Dark2[c(1, 2, 3, 5)], limits = c("SB dose 1", "SB dose 5", "SB dose 10", "SB dose 15"), labels = c("1", "5", "10", "15")) +
  # add here vlines for all exposure types
  geom_vline(data = pred_i_exp |> dplyr::filter( exposure_type != "EHT"), aes(xintercept = exposure_value, colour = exposure_type)) +
  geom_segment(aes(x = 15, y = 4, xend = 15, yend = 5, colour = "SB dose 15"), show.legend = FALSE) +
  geom_segment(aes(x = 15, y = 6, xend = 15, yend = 7, colour = "SB dose 15"), show.legend = FALSE) +
  # dummy line for improved legend
  geom_vline(aes(xintercept = -10, colour = "SB dose 15")) +
  # the limits remove a few draws but ensure the density plots have more reasonable i.e. more normalised height
  scale_x_sqrt(breaks = c(0, 1, 5, 10, 20, 50, 100, 200), limits = c(0,quantile(draws_i_spread |> filter(Group %in% c(1, 2, 3, 4, 5, 6, 7)) |> pull(lethal_dose), probs = 0.9995))) +
  coord_cartesian(xlim =c(0,250)) +
  xlab("Lethal dose [Dose]") +
  ylab("Assay pair / Density") +
  ggtitle("Lethal exposure [Dose]") +
  scale_y_discrete(limits=rev)


# exposure dose plot  
p_exposure_dose <- ggplot(data = draws_i_spread_exp |> filter(exposure_type %in% c("SB dose 1", "EHT")), aes(x = exposure_value, fill = exposure_type)) +
  scale_color_manual(name = "Exposure distribution in", values = colours_Dark2[c(1, 4)], limits = c("SB dose 1", "EHT"), labels = c("SB", "EHT")) +
  scale_fill_manual(name = "Exposure distribution in", values = colours_Dark2[c(1, 4)], limits = c("SB dose 1", "EHT"), labels = c("SB", "EHT")) +
  tidybayes::stat_slab(normalize = "panels", show.legend = FALSE) +
  # the limits remove a few draws but ensure the density plots have more reasonable i.e. more normalised height
  scale_x_sqrt(breaks = c(0, 0.25, 0.5, 1, 2, 5, 10, 15), limits = c(0, quantile(draws_i_spread_exp$lethal_dose, probs = 0.9995))) +
  coord_cartesian(xlim =c(0,250), ylim = c(0,0.8)) +
  ggtitle("Exposure [Dose]") + 
  xlab("Exposure [Dose]") +
  ylab("Density") +
  theme(axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()
        ) +
  guides(fill = guide_legend(reverse = TRUE))

p_exposure_dose_v <- p_exposure_dose +
  coord_cartesian(xlim =c(0,2), ylim = c(0,0.85))
  
# comb plot
p_comb_v <- wrap_plots(p_exposure_dose_v, p_lethal_dose, nrow = 1, ncol = 2, widths = c(1,3))  +
  plot_annotation(tag_levels = list(c("A", "B"))) +
  plot_layout(guides = "collect") &
  theme(legend.position = 'bottom') &
  labs(title = NULL) 

ggsave(file.path("plots_mechmodel_article", "lethalDoseANDExposure_v.png"), p_comb_v, width = 12, height = 6)

