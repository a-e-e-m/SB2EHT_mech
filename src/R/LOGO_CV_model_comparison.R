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
is <- 194


# copying files to local
for (k in is){
  cmd <- sprintf(
    "rsync -avz -e ssh %s %s",
    shQuote(paste0(
      "denadr00@transfer12.scicore.unibas.ch:/scicore/home/chitnis/denadr00/lAIRepmi/BA2EHT/fitting/",
      models2run$run[k]
    )),
    shQuote("./fitting/")
  )
  system(cmd)
}



# collect information on runs
list_dir_i <- vector("list", length(is))
list_trials <- vector("list", length(is))
list_nT <- vector("list", length(is))
list_H_all <- vector("list", length(is))
list_B_all <- vector("list", length(is))

for (k in seq_along(is) ){
  list_dir_i[[k]] <- paste0("fitting/", models2run$run[is[k]]) # folder
  list_trials[[k]] <- readRDS(file.path(list_dir_i[[k]], 'trials.rds'))
  list_nT[[k]] <- list_trials[[k]]$nT 
  list_H_all[[k]] <- readRDS(file = file.path(list_dir_i[[k]], "H_all.rds"))
  list_B_all[[k]] <- readRDS(file = file.path(list_dir_i[[k]], "B_all.rds"))
}

# get data after checking data is the same across all runs
ref_H <- list_H_all[[1]]
ref_B <- list_B_all[[1]]
if ( all(vapply(list_H_all[-1], function(x) identical(x, ref_H), logical(1))) ) {
  H_all <- list_H_all[[1]]
} else {stop("Runs are based on different data!")}

if ( all(vapply(list_B_all[-1], function(x) identical(x, ref_B), logical(1))) ) {
  B_all <- list_B_all[[1]]
} else {stop("Runs are based on different data!")}

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
    q975_f = qbeta(0.975, shape1 = AF_h + 1, shape2 = A_h - AF_h + 1)
  )

# get predictions
# reduce predictions to interventions that are actually in the data
# and add int_dose_available
list_pred_sum <- vector("list", length(is))
for (k in seq_along(is) ){
  list_pred_sum[[k]] <- readRDS(file = file.path(list_dir_i[[k]], "predictions_summary.rds")) |>
    mutate(run = is[k]) |>
    relocate(run) |>
    inner_join(data_summary)
}

# put together, add actual (real data, run == 0) and prepare for plotting
pred_sum <- do.call(rbind, list_pred_sum) |>
  rbind(mutate(data_summary, run = 0L, 
               pred_median = MLE, pred_q025 = q025, pred_q975 = q975,
               pred_f_median = MLE_f, pred_f_q025 = q025_f, pred_f_q975 = q975_f) ) 


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
    
    secondai_ID = factor(
      secondai_ID,
      levels = sort(unique(H_all$secondai_ID)),
      labels = lookup_vec_secondai
    ),
    
    run = factor(
      run,
      levels = c(0, is),
      labels = c("actual", models2run$run[is])
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
  pred_sum2,
  aes(x = secondai_ID, y = pred_median, color = run)
)  +
  geom_errorbar(
    aes(ymin = pred_q025, ymax = pred_q975),
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

# for highlighting groups with IDSB data available if some runs do 
if ( any(str_detect( models2run$bioassay_type[is], "int" ) ) ) {
  run_id_with_IDSB <- match( str_detect( models2run$bioassay_type[is], "int" ), TRUE )
  highlight_df <- pred_sum2 %>%
    distinct(group_number, facet_label) |>
    left_join(unique(select(list_B_all[[run_id_with_IDSB]], group_number, int_dose_available)))
  
  p1 <- p1 +
    geom_rect(
    data = highlight_df,
    aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = NA,
    color = "red",
    linewidth = 1.2
  )
}
  
# save plot
ggsave(file = file.path("plots_LOGO_CV", paste0("LOGO-CV-comp_", paste(is, collapse = "-"), ".png")), plot=p1, width = 18, height = 10)
saveRDS(p1, file = file.path("plots_LOGO_CV", paste0("LOGO-CV-comp_", paste(is, collapse = "-"), "rds")))


# produe feeding prob plot if applicable for at least some of the runs
if ( any(models2run$with_feeding[is]) ){
  p1_f <- ggplot(
    pred_sum2,
    aes(x = secondai_ID, y = pred_f_median, color = run)
  ) +
    geom_errorbar(
      aes(ymin = pred_f_q025, ymax = pred_f_q975),
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
  
  # for highlighting groups with IDSB data available if some runs do 
  if ( any(str_detect( models2run$bioassay_type[is], "int" ) ) ) {
    run_id_with_IDSB <- match( str_detect( models2run$bioassay_type[is], "int" ), TRUE )
    highlight_df <- pred_sum2 %>%
      distinct(group_number, facet_label) |>
      left_join(unique(select(list_B_all[[run_id_with_IDSB]], group_number, int_dose_available)))
    
    p1_f <- p1_f +
      geom_rect(
        data = highlight_df,
        aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
        inherit.aes = FALSE,
        fill = NA,
        color = "red",
        linewidth = 1.2
      )
  }
  
  # save plot
  ggsave(file = file.path("plots_LOGO_CV", paste0("LOGO-CV-comp_", paste(is, collapse = "-"), "_", "f.png")), plot=p1_f, width = 18, height = 10)
  saveRDS(p1_f, file = file.path("plots_LOGO_CV", paste0("LOGO-CV-comp_", paste(is, collapse = "-"), "_", "f.rds")))
}

