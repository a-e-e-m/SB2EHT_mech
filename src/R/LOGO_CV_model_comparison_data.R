## Script to compare LOGO-CV predictions per model

## config
################################
library(tidyverse)
library(tidybayes)
library(ggplot2)
library(extraDistr)
library(patchwork)

# source script containing needed functions
functionsfolder <- file.path('./src/R/functions')
files.source <- list.files(functionsfolder)
invisible(sapply(files.source, function(x) source(paste0(functionsfolder, "/", x))))

# read file with all model runs
models2run <- read.csv(file = './models2run.csv', sep = ';', stringsAsFactors = FALSE)

# model runs to look at
is <- 187

# get data after checking data is the same across all runs
H_all <- readRDS(file = file.path(paste0("fitting/", models2run$run[is]), "H_all.rds"))
B_all <- readRDS(file = file.path(paste0("fitting/", models2run$run[is]), "B_all.rds"))

# summarise data
data_summary <- H_all |>
  filter(!control) |>
  group_by(group_number, secondai_ID, country, site, year, pyrethroid) |>
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
    q975_f = qbeta(0.975, shape1 = AF_h + 1, shape2 = A_h - AF_h + 1),
    run = "actual"
  )

pred_sum <- data_summary

## facet, wrapped, plot

# LLIN types
lookup_vec_secondai <- H_all %>%
  filter(!is.na(secondai)) %>%
  distinct(secondai_ID, secondai) %>%
  { setNames(.$secondai, .$secondai_ID) }
lookup_vec_secondai[1] <- "pyrethroid only"

# prepare data
pred_sum2 <- pred_sum %>%
  ungroup() %>%
  mutate(
    group_number = as.integer(group_number),
    group_number = factor(group_number),
    run = "actual",
    secondai_ID = factor(
      secondai_ID,
      levels = sort(unique(H_all$secondai_ID)),
      labels = lookup_vec_secondai
    ),
    
    facet_label = paste0(
      "Group ", group_number, "\n",
      country, ", ", site, ", ", year, "\n",
      pyrethroid
    )
  ) 

facet_levels <- pred_sum2 %>%
  mutate(group_number = as.integer(group_number)) |>
  distinct(group_number, facet_label) %>%
  arrange(group_number) %>%
  pull(facet_label)

pred_sum2 <- pred_sum2 %>%
  mutate(
    facet_label = factor(facet_label, levels = facet_levels)
  )

pd <- position_dodge(width = 0.6)

p1 <- ggplot(
  filter(pred_sum2, run == "actual"),
  aes(x = secondai_ID, y = MLE, color = run)
)  +
  geom_errorbar(
    aes(ymin = q025, ymax = q975),
    position = pd,
    width = 0.15
  ) +
  geom_point(
    position = pd,
    size = 1.8
  ) +
  facet_wrap(~ facet_label, scales = "free_x") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Treatment",
    y = "Predicted treatment effect",
    color = "Run"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

  
# save plot
ggsave(file = file.path( paste0("LOGO-CV-comp_", is, "_data.png")), plot=p1, width = 18, height = 10)
saveRDS(p1, file = file.path( paste0("LOGO-CV-comp_", is, "_data.rds")))


# produe feeding prob plot if applicable for at least some of the runs
if ( any(models2run$with_feeding[is]) ){
  p1_f <- ggplot(
    pred_sum2, 
    aes(x = secondai_ID, y = MLE_f, color = run)
  ) +
    geom_errorbar(
      aes(ymin = q025_f, ymax = q975_f),
      position = pd,
      width = 0.15
    ) +
    geom_point(
      position = pd,
      size = 1.8
    ) +
    facet_wrap(~ facet_label, scales = "free_x") +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      x = "Intervention",
      y = "Feeding probability",
      color = "Run"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  
  # save plot
  ggsave(file = file.path( paste0("LOGO-CV-comp_", is, "_", "f_data.png")), plot=p1_f, width = 18, height = 10)
  saveRDS(p1_f, file = file.path( paste0("LOGO-CV-comp_", is, "_", "f_data.rds")))
}

