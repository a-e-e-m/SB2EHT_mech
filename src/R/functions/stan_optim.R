stan_optim <- function(model, with_feeding, BA_only, add_data = NA, run, iter = 1) {
  
  ## config
  library('tidyverse')
  library('rstan')
  rstan_options(auto_write = TRUE)
  library("loo")
  options(mc.cores = parallel::detectCores())


  # select stan model
  model_file <- paste0("./src/stan/", model, ".stan", collapse = NULL)
  m <- stan_model(file = model_file)
  
  # load data
  B <- readRDS(file.path("fitting", run,  "B.rds"))
  S_b <- nrow(B)
  
  if (!BA_only){
    H <- readRDS(file.path("fitting", run,  "H_all.rds"))
  }else{H <- tibble()}
  S_h <- nrow(H)
  
  if (with_feeding){
    H_f <- readRDS(file.path("fitting", run,  "H_f_all.rds"))
  }else{H_f <- tibble()}
  S_h_f <- nrow(H_f)
  
  # add additional, readily stan readable data, if given
  # if given, copy also to run folder 
  if (!is.na(add_data)){
    add_data_list <- readRDS(file.path("data", "data_add", add_data))
    saveRDS(add_data_list, file = file.path("fitting", run, add_data))
  } else {
    add_data_list <- list()
  }
  
  trials <- readRDS(file.path("fitting", run,  "trials.rds"))
  
  input_stan <- c(as.list(B), as.list(H), as.list(H_f), list(S_b = S_b, S_h = S_h, S_h_f = S_h_f), trials, add_data_list )
  
  # LOO and exactLOO
  input_stan[["LOO"]] <- 0
  input_stan[["exactLOO"]] <- 0
  
  # add empty data for exactLOO
  input_stan[["LOdataset"]] <- numeric(0)
  input_stan[["S_h_LO"]] <- numeric(0)
  input_stan[["N_h_LO"]] <- numeric(0)
  input_stan[["D_h_LO"]] <- numeric(0)
  input_stan[["S_h_f_LO"]] <- numeric(0)
  input_stan[["A_LO"]] <- numeric(0)
  input_stan[["AF_LO"]] <- numeric(0)
  
  saveRDS(input_stan, file=file.path("fitting", run, "input_stan.rds"))
  
  # run stan optimizer  
  if (iter > 1){
    est_iter <- vector(mode = "list", length = iter) 
    loglik_iter <- vector(length = iter)
    for (k in 1:iter){
    est_iter[[k]] <- optimizing(m, data = input_stan, as_vector = FALSE)  
    loglik_iter[k] <- est_iter[[k]][[2]]
    }
    est <- est_iter[[which.max(loglik_iter)]]
  } else {
    est <- optimizing(m, data = input_stan, as_vector = FALSE)  
  }

  ## saving results
  savetopath_est <- file.path("fitting", run, "est.rds")
  saveRDS(est, file=savetopath_est)
}
