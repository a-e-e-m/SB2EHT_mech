## Script to combine plots on LOGO-CV

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
 
 
## NOTE: model order always is Nash_nolocation, Nash_norandom, DDSB_idea2_prob, mechmodel

# Model estimates VS LOGO-CV across models
 is = c(147, 146, 156, 140)
 
 plots_LOGO_CV <- list()
 plots_postpc <- list()
 
 for (j in seq_along(is)) {
   i <- is[j]  # index into your models2run$run
   dir_i <- file.path("fitting", models2run$run[i])
   plots_LOGO_CV[[j]] <- readRDS(file.path(dir_i, "LOGO-CV-actVSpred.rds")) +
     theme_gray(base_size = 12)
   plots_postpc[[j]] <- readRDS(file.path(dir_i, "postpc_H_treat.rds")) +
     theme_gray(base_size = 12)
 }
 
 # Create patchwork rows from your plot lists
 top_row <- wrap_plots(plots_postpc, nrow = 1, guides = "collect")
 
 bottom_row <- wrap_plots(plots_LOGO_CV, nrow = 1, guides = "collect")
 
 # Combine with layout
 combined_plot <- top_row / bottom_row +
   plot_annotation(
     tag_levels = list(c("I.A", "I.B", "I.C", "I.D",
                         "II.A", "II.B", "II.C", "II.D"))
   )
 
 ggsave(file.path("plots_LOGO_CV", "postpcVSlogocv.png"), combined_plot, width = 20, height = 8)
 
 
# LOGO-CV all data vs only West across models
 is_alldata = c(146, 156, 140)
 is_onlywest = c(160, 161, 162) 

 plots_LOGO_CV_alldata <- list()
 plots_LOGO_CV_onlywest <- list()
 
 for (j in seq_along(is_alldata)) {
   i <- is_alldata[j]  # index into your models2run$run
   dir_i <- file.path("fitting", models2run$run[i])
   plots_LOGO_CV_alldata[[j]] <- readRDS(file.path(dir_i, "LOGO-CV-actVSpred.rds")) +
     theme_gray(base_size = 12)
 }
 
 
 for (j in seq_along(is_onlywest)) {
   i <- is_onlywest[j]  # index into your models2run$run
   dir_i <- file.path("fitting", models2run$run[i])
   plots_LOGO_CV_onlywest[[j]] <- readRDS(file.path(dir_i, "LOGO-CV-actVSpred.rds")) +
     theme_gray(base_size = 12)
 }
 
 
 
 # Create patchwork rows from your plot lists
 top_row <- wrap_plots(plots_LOGO_CV_alldata, nrow = 1, guides = "collect")
 
 bottom_row <- wrap_plots(plots_LOGO_CV_onlywest, nrow = 1, guides = "collect")
 
 # Combine with layout
 combined_plot_alldataVSonlywest <- top_row / bottom_row +
   plot_annotation(
     tag_levels = list(c("I.B", "I.C", "I.D",
                         "II.B", "II.C", "II.D"))
   )
 
 ggsave(file.path("plots_LOGO_CV", "LOGO_CV_alldataVSonlywest.png"), combined_plot_alldataVSonlywest, width = 16, height = 8)


# 