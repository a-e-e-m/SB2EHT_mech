# config and prepare
###############################
###############################

#libraries
library("rstan")
rstan_options(auto_write = TRUE) # set FALSE for cluster
library("loo")
library("bayesplot")
#library("Rfast")
library("HDInterval")
library("dplyr")
library("tidyr")
library("tibble")
library("posterior")
library("readr")
library("stringr")

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

  # loading data
  #############
  dir_i <- file.path("fitting", models2run$run[i]) # folder
  fit_i <- readRDS(file.path(dir_i, 'fit.rds'))
  
  
  # reading in parameters relevant for this model
  ###############################################
  parameters_all <- fit_i@model_pars
  
  # select strictly only sampled parameters
  parameters_stansection <- str_match(fit_i@stanmodel@model_code, "(?s)parameters\\{(.*?)\\}")[2]
  parameters_strict <- parameters_all[str_detect(parameters_stansection, paste0(parameters_all, "(;|\\[)"))]  
  
  # select only parameters defined in parameters and transformed parameters block,
  # by excluding all generated quantities.
  # #Warning: This assumes the first 'parameter' defined in the generated quantities block is 'log_lik'.
  # log_lik_index <- grep(x = parameters_all, "log_lik")
  # parameters_all <- parameters_all[1:log_lik_index-1]
  
  parameters_strict_df <- list_to_df(fit_i@par_dims) %>% 
    rename(parameter = name, dim = list.element) %>% 
    #column_to_rownames("parameter") %>%
    relocate(parameter, dim) %>%
    mutate(dim = as.integer(dim)) %>% 
    mutate(dim = case_when(is.na(dim) ~ 1L, TRUE ~ dim)) %>%
    filter(parameter %in% parameters_strict)
  
  parameters_strict_sepdim <- sapply(parameters_strict, FUN = function(x){
    dim_x <- parameters_strict_df |> filter(parameter == x) |> pull(dim)
    if(dim_x == 1){
      x
    } else if (dim_x > 1){
      paste0(x, "[", 1:dim_x, "]")
      }})
  parameters_strict_sepdim <- unlist(parameters_strict_sepdim)
  
  # first level (FL)
  parameters_FL_scalar <- parameters_strict_df %>% filter(dim == 1) %>% pull(parameter)
  
  
  # diagnostics
  #############
  # careful: now sum_df only contains the sampled parameters, not the transformed and no generated quantities
  sum_df <- monitor(rstan::extract(fit_i, pars = parameters_strict, permuted = FALSE))
  sum_df <- as_tibble(sum_df, rownames = "variable")
  saveRDS(sum_df, file = file.path(dir_i, 'sum_df.rds') )
  rhat_max <- max(sum_df$Rhat)
  bulkess_min <- min(sum_df$Bulk_ESS)
  tailess_min <- min(sum_df$Tail_ESS)
  num_divergent_i <- get_num_divergent(fit_i)
  num_maxtreedepth_i <- get_num_max_treedepth(fit_i)

  write_lines( list(paste0("The highest Rhat value across all parameters was ", rhat_max, "."),
                    paste0("The smallest bulk ESS across all parameters was ", bulkess_min, "."),
                    paste0("The smallest tail ESS across all parameters was ", tailess_min, "."),
                    paste0("There were ", num_divergent_i, " divergent transitions."),
                    paste0("The maximal tree depth was hit ", num_maxtreedepth_i, " times.")
                    ), file = file.path(dir_i, 'diagnostics.txt'))

  # if (rhat_max >= 1.05) {write_lines("Rhat is too big.", file = file.path(dir_i, 'diagnostics.txt'), append = TRUE)}
  # if (bulkess_min <=400) {write_lines("Bulk ESS is too small.", file = file.path(dir_i, 'diagnostics.txt'), append = TRUE)}
  # if (tailess_min <=400) {write_lines("Tail ESS is too small.", file = file.path(dir_i, 'diagnostics.txt'), append = TRUE)}
  if (rhat_max >= 1.05) {write_lines(paste0(c("The badly mixing parameters (Rhat > 1.05) were: ", sum_df |> filter(Rhat > 1.05) |> pull(variable))), file = file.path(dir_i, 'diagnostics.txt'), append = TRUE)}


  # chain plots
  ###############################################
  # scalar parameters of first level model
  if (length(parameters_FL_scalar) > 0){
    jpeg(file = file.path(dir_i, paste0('traceplot_','parameters_FL_scalar','.jpeg')),
       width = 800, height = length(parameters_FL_scalar) * 100)

    rstan::traceplot(fit_i, pars = parameters_FL_scalar, nrow = length(parameters_FL_scalar), inc_warmup = TRUE)
    dev.off()
  }

  
  # parameters with bad mixing
  parameters_badmixing <- sum_df |> filter(Rhat > 1.05) |> pull(variable)
  if (length(parameters_badmixing) > 0){
  jpeg(file = file.path(dir_i, paste0('traceplot_','parameters_badmixing','.jpeg')),
       width = 800, height = length(parameters_badmixing) * 100)
  
  rstan::traceplot(fit_i, pars = parameters_badmixing, nrow = length(parameters_badmixing), inc_warmup = TRUE)
  dev.off()
  }

  ## pairs plots
  ######################################################
  # FL_scalar
  if (length(parameters_FL_scalar) > 1){
    jpeg(file = file.path(dir_i, paste0('pairs_','parameters_FL_scalar','.jpeg')),
         width = length(parameters_FL_scalar)*300, height = length(parameters_FL_scalar)*200)
  
    pairs(fit_i, pars = parameters_FL_scalar)
    dev.off()
  }

