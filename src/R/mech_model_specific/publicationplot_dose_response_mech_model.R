# generating publication plots 
# BA dose-response data and prediction for both BA and joint model


# config and prepare
###############################
###############################

#libraries
library("HDInterval")
library("tidyr")
library("dplyr")
library("rstan")
library("stringr")
library("scales")
library("ggplot2")
library("RColorBrewer")
library("lemon")

# source script containing needed functions
functionsfolder <- file.path('./src/R/functions')
files.source <- list.files(functionsfolder)
invisible(sapply(files.source, function(x) source(paste0(functionsfolder, "/", x))))
source("./src/R/mech_model_specific/function_load_predictions.R")

models2run <- read.csv(file = './models2run.csv', sep = ';', stringsAsFactors = FALSE)


# loading real data
#############
#BA only model run
i <- 3
# joint model run
i_joint <- 4
dir_i <- paste0("fitting/", models2run$run[i]) # folder
B_all <- readRDS(file.path(dir_i, 'B_all.rds'))


# group_number_conversion with intdose BA only data set
B_all_BAonly <- readRDS(file.path(paste0("fitting/", models2run$run[i]) , 'B_all.rds'))
group_number_conversion <- B_all_BAonly |> 
  filter(!control) |>
  select(country, site, year, insecticide, group_number) |>
  unique() |>
  rename(group_number_conversion = group_number)

# insecticide_conversion for controls
insecticide_conversion <- B_all_BAonly |> 
  filter(!control) |>
  select(country, site, year, insecticide, group_number) |>
  unique() |>
  rename(insecticide_conversion = insecticide)

B_all <- B_all |>
  left_join(insecticide_conversion)

# load model predictions from BA only model
output_BA <-  modelpred_4_pubplot(i)

# load model predictions from joint model 
output_joint <-  modelpred_4_pubplot(i_joint)


# combine predictions and data    
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
  left_join(output_BA[["prob_D_b_summary"]], by = c("group_number", "control", "treat", "times_disc_dose")) |>
  left_join(output_joint[["prob_D_b_summary"]], by = c("group_number", "control", "treat", "times_disc_dose")) |>
  left_join(insecticide_conversion |> mutate(year = as.character(year)))


# add facet_col variable for nice nested facet plot
B_all <- B_all |> 
  filter(int_dose_available == 1, country == "BurkinaFaso") |> 
  mutate(facet = paste0("bold(", "Assay pair ", group_number, ": year ", year),
         facet_col = group_number %% 2,
         facet_col = case_when(facet_col == 1 ~ "first year", facet_col == 0 ~ "second year")
  )

B_all_summary <- B_all_summary |> 
  filter(int_dose_available == 1, country == "BurkinaFaso") |> 
  mutate(facet = paste0("bold(", "Assay pair ", group_number, ": year ", year),
         facet_col = group_number %% 2,
         facet_col = case_when(facet_col == 1 ~ "first year", facet_col == 0 ~ "second year")
  )

# dose response curves per group with facets (only for groups with intensity dose data)
###################
# compute predictions for dose-response
times_disc_dose <- seq(0, 15, 0.05)

# lookup for group_number for int_dose groups
# contains quick hack to exclude data from Tanzania
lookup_group_number <- B_all_summary |>
  filter(int_dose_available == 1, country == "BurkinaFaso") |>
  ungroup() |>
  filter(insecticide != "none") |>
  select(group_number, country, site, year, insecticide) |>
  unique()

# generate index sample for plotting sample of dose-response curves only
iter_sample_BA <- sample(1:max(output_BA[["list_dfs_doseresponse"]][[1]]$iter), 500, replace = FALSE)
iter_sample_joint <- sample(1:max(output_joint[["list_dfs_doseresponse"]][[1]]$iter), 500, replace = FALSE)


# initialise plot
p_B_facet <- ggplot() +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2")

# add dose-response curves for each group (median curve and random sample of curves)
for (ii in seq_along(sort(unique(lookup_group_number$group_number)))){
  p_B_facet <- p_B_facet +
    geom_line(data = output_BA[["list_dfs_doseresponse"]][[ii]] |> 
                filter(iter %in% iter_sample_BA) |> 
                mutate(facet = paste0("bold(", "Assay pair ", group_number, ": year ", year),
                       facet_col = group_number %% 2,
                       facet_col = case_when(facet_col == 1 ~ "first year", facet_col == 0 ~ "second year"),
                       insecticide_conversion = insecticide
                ), 
              aes(x = times_disc_dose, y = p, group = iter, colour = "ID-SB"), alpha = 0.01) +
    geom_line(data = output_BA[["list_dfs_doseresponse_summary"]][[ii]] |> 
                mutate(facet = paste0("bold(", "Assay pair ", group_number, ": year ", year),
                       facet_col = group_number %% 2,
                       facet_col = case_when(facet_col == 1 ~ "first year", facet_col == 0 ~ "second year"),
                       insecticide_conversion = insecticide
                ), 
              aes(x = times_disc_dose, y = p_median, colour = "ID-SB"), linewidth = 0.5) +
    geom_line(data = output_joint[["list_dfs_doseresponse"]][[ii]] |> filter(iter %in% iter_sample_BA) |>
                mutate(year = as.numeric(year)) |>
                left_join(group_number_conversion) |>
                filter(group_number_conversion < 8) |>
                ungroup() |>
                select(!group_number) |>
                rename(group_number = group_number_conversion) |> 
                mutate(facet = paste0("bold(", "Assay pair ", group_number, ": year ", year),
                       facet_col = group_number %% 2,
                       facet_col = case_when(facet_col == 1 ~ "first year", facet_col == 0 ~ "second year"),
                       insecticide_conversion = insecticide
                ), 
              aes(x = times_disc_dose, y = p, group = iter, colour = "DD-SB + ID-SB + EHT"), alpha = 0.01) +
    geom_line(data = output_joint[["list_dfs_doseresponse_summary"]][[ii]] |>
                mutate(year = as.numeric(year)) |>
                left_join(group_number_conversion) |>
                filter(group_number_conversion < 8) |>
                ungroup() |>
                select(!group_number) |>
                rename(group_number = group_number_conversion) |> 
                mutate(facet = paste0("bold(", "Assay pair ", group_number, ": year ", year),
                       facet_col = group_number %% 2,
                       facet_col = case_when(facet_col == 1 ~ "first year", facet_col == 0 ~ "second year"),
                       insecticide_conversion = insecticide
                ), 
              aes(x = times_disc_dose, y = p_median, colour = "DD-SB + ID-SB + EHT"), linewidth = 0.5)
}


# add real data and empirical estimates 
p_B_facet <- p_B_facet +
  geom_pointrange(data = B_all_summary, 
                  mapping = aes(x = times_disc_dose, y = median, ymin = q025, ymax = q975), colour = "black", size = 0.2, alpha = 0.9) +
  geom_jitter(data = B_all, 
              mapping = aes(x = times_disc_dose, y = D_b / N_b), colour = "black", size = 0.5, alpha = 0.9) +
  # geom_text(data = B_all_summary |> filter(treat == 0),
  #                  mapping = aes(x = 4, y = 0.85, label = paste0(year, " (assay pair", group_number, ")")),
  #                  hjust   = -0.1,
  #                  vjust   = -1,
  #                 fontface = "bold") +
  ggh4x::facet_nested(insecticide_conversion ~ site + facet_col) +
  scale_y_continuous(labels = percent) +
  xlab("Insecticide challenge in ID-SB [Dose]") + ylab("Mortality in ID-SB [Probability]") +
  labs(colour = "Fitted to data from", title = NULL) 


# Convert the plot to a gtable after removing the legend
g <- ggplotGrob(p_B_facet + theme(legend.position = "none"))

# Inspect grobs
# gtable::gtable_show_layout(g)
# gtable_show_names(g)
g$layout$name

# Remove empty facet and move axis
g$grobs[[grep("axis-b-4-1", g$layout$name)]] <- g$grobs[[grep("axis-b-4-2", g$layout$name)]]
g$grobs[[grep("axis-b-4-2", g$layout$name)]] <- grid::nullGrob()
g$grobs[[grep("panel-2-4", g$layout$name)]] <- grid::nullGrob()

# put legend back in, in empty facet
g <- reposition_legend(g, legend =  g_legend(p_B_facet), 'bottom', panel='panel-2-4')


# save plot
ggsave(file.path("plots_mechmodel_article", "dose-response-curves_onlypublished_quantile_new.png"), g, width = 12, height = 6)

