validfit <- function(model, init_function = "random", with_feeding, BA_only, add_data = NA, LOO, exactLOO, pareto_k_threshold, LOO_CV, LOO_CV_groups, chains, warmup, iter, run, control_list_name = "") {
  # This function return a list containing the stanfit object for the given model fitted to the whole data, 
  # the loo() output for this with PSIS-approximated elpd (expected log posterior density) 
  # replaced by the exact elpd for the data points with pareto-k value above pareto_k_threshold, 
  # and a vector with the indices of the data points with pareto-k value above pareto-k_threshold.
  # pareto_k_threshold = 0.5 or at most 0.7 is recommended
  # Note: exactLOO means exact LOO only for data points with pareto k parameter above pareto_k_threshold. 
  # Hence, for exactLOO==1, LOO==1 is required and pareto_k_threshold.
  
  ## config
  if (LOO == 1){
    if (chains <= 1) {stop("needs multiple chains")}
  }
  
  library('tidyverse')
  library('rstan')
  rstan_options(auto_write = TRUE)
  library("loo")
  #options(mc.cores = parallel::detectCores())
  options(mc.cores = chains)

  # select stan model
  model_file <- paste0("./src/stan/", model, ".stan", collapse = NULL)
  m <- stan_model(file = model_file)
  
  # load data
  B <- readRDS(file.path("fitting", run,  "B_all.rds"))
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
  
  saveRDS(input_stan, file=file.path("fitting", run, "input_stan.rds"))
  
  # LOO and exactLOO
  input_stan[["LOO"]] <- LOO
  input_stan[["exactLOO"]] <- 0
  
  # add empty data for exactLOO
  input_stan[["LOdataset"]] <- numeric(0)
  input_stan[["S_h_LO"]] <- numeric(0)
  input_stan[["N_h_LO"]] <- numeric(0)
  input_stan[["D_h_LO"]] <- numeric(0)
  input_stan[["S_h_f_LO"]] <- numeric(0)
  input_stan[["A_LO"]] <- numeric(0)
  input_stan[["AF_LO"]] <- numeric(0)
  
  # set initial values for all chains based on passed init_function
  if (init_function == "random" | init_function == ""){
    inits_list <- "random"
  } 
  else if (init_function == "est_all"){
    est <- readRDS(file.path("fitting", run,  "est.rds"))
    inits_list <- list()
    for (k in 1:chains){
      inits_list = append(inits_list, list(as.list(est$par)))
    }
  } 
  else {
    inits_list <- list()
    for (k in 1:chains){
      inits_list = append(inits_list, list(do.call(init_function, list(k))))
    }
  }
  
  # load list to control stan sampler if provided
  if (control_list_name != '' && !is.na(control_list_name)){
    control_list <- readRDS(file = file.path("src", "R", "stan_control_lists", paste0(control_list_name, ".rds")))
    saveRDS(control_list, file = file.path("fitting", run, paste0(control_list_name, ".rds")))
  } else {
    control_list <- list()
  }
  
  # run stan
  fit <- sampling(
    m,
    data = input_stan,
    chains = chains,
    warmup = warmup,
    iter = iter,
    control = control_list,
    show_messages = TRUE,
    init_r = 2,
    init = inits_list,
    refresh = max(iter/50, 1)
  )
  
  # save fit object
  savetopath_fit <- file.path("fitting", run, "fit.rds")
  saveRDS(fit, file=savetopath_fit)

  ## LOO-CV by PSIS with using group association of left-out data point for prediction
  if (LOO==1){
    log_lik <- loo::extract_log_lik(fit, merge_chains = FALSE) # Extract pointwise log-likelihood

    # as of loo v2.0.0 we can optionally provide relative effective sample sizes
    # when calling loo, which allows for better estimates of the PSIS effective
    # sample sizes and Monte Carlo error
    r_eff <- relative_eff(exp(log_lik))
    loo_ <- loo(fit, r_eff = r_eff, save_psis=TRUE) # CHECK ORDERING OF ELEMENTS IN PSIS_OBJECT IN TERMS OF SAMPLE
    # extracting the indices of data points with pareto-k estimate above pareto_k_threshold
    I <- pareto_k_ids(loo_, threshold = pareto_k_threshold)
  }

  ## exact LOO-CV if LOO_CV == 1 (LOO==1 not required) or if exactLOO ==1 for groups with pareto-k estimate above pareto_k_threshold by PSIS (LOO==1 required)
  if (exactLOO==1 | LOO_CV==1){
    # throw error if EHT data excluded by run specification
    if (BA_only){stop("EHT data excluded by run specification. CrossValidation only possible if EHT data included.")}
    
    # if exactLOO then I from pareto k threshold is used, otherwise I is defined here
    if (LOO_CV==1 ){
      if (str_detect(LOO_CV_groups, "all")){
        I <- seq(1, trials[["nT"]])
      } else {
        I <- readRDS(LOO_CV_groups)
      }
    }
    
  # refitting model length(I) times with each time leaving out one i in I and saving corresponding
  # exact expected log posterior density (elpd)
  
  #saving object names before going into the loop
  objects_before_loop <- ls()
  
    for (i in I) {
      rm(list = c("input_stan", "B", "H", "H_f"))
        
      # load data
      B <- readRDS(file.path("fitting", run,  "B.rds"))
      S_b <- nrow(B)
      H <- readRDS(file.path("fitting", run,  "H.rds"))
      
      if (with_feeding){
        H_f <- readRDS(file.path("fitting", run,  "H_f.rds"))
      }else{H_f <- tibble()}
      
      # add additional, readily stan readable data, if given
      if (add_data != ''){
        add_data_list <- readRDS(file.path("data", "data_add", add_data))
      } else {
        add_data_list <- list()
      }
      
      # add grouping info
      trials <- readRDS(file.path("fitting", run,  "trials.rds"))
      
      # put apart EHT treatment data point of data set to leave out as B_LO, H_LO, H_f_LO
      # and rename corresponding variables
      H_LO <- H %>% 
        filter(T_h == i, treat_h == 1) %>% 
        rename(N_h_LO = N_h, D_h_LO = D_h) 
      
      H <- H %>% 
        filter(!(T_h == i & treat_h == 1))     
      
      S_h_LO <- as.array(nrow(H_LO))
      
      H_LO <- as.list(H_LO)
      
      # and convert data to array of dim 1 if only one entry (because of stan input requirements)
      if (S_h_LO == 1){
        H_LO <- lapply(H_LO, as.array)
      }
      
      # correct number of data points in data to fit to
      S_h<- nrow(H)
      
      # same for feeding data if included
      if (with_feeding){
        H_f_LO <- H_f %>% 
          filter(T_h_f == i, treat_h_f == 1) %>% 
          rename(A_LO = A, AF_LO = AF)
        
        H_f <- H_f %>% 
          filter(!(T_h_f == i & treat_h_f == 1)) 
        
        H_f_LO <- as.list(H_f_LO)
        
        if (S_h_f_LO == 1){
          H_f_LO <- lapply(H_LO, as.array)
        }
      }else{
        H_f_LO <- tibble()
      }
      
      S_h_f_LO <- as.array(nrow(H_f_LO))
      S_h_f <- nrow(H_f)

      
      # convert to list
      input_stan <- c(as.list(B), as.list(H), as.list(H_f), 
                      list(S_b = S_b, S_h = S_h, S_h_f = S_h_f), 
                      trials, add_data_list, 
                      H_LO, H_f_LO, 
                      list(S_h_LO = S_h_LO, S_h_f_LO = S_h_f_LO),
                      trials, add_data_list, 
                      list(LOdataset = as.array(i) ))
      
      # add LOO and exactLOO switches
      input_stan[["LOO"]] <- 0
      input_stan[["exactLOO"]] <- 1
      
      
      # save to run folder
      saveRDS(input_stan, file=file.path("fitting", run, paste0("input_stan_LOO_", i,".rds")))
      
      # fit
      fit_lio <- sampling(
        m,
        data = input_stan,
        chains = chains,
        warmup = warmup,
        iter = iter,
        control = list(adapt_delta = 0.8),
        show_messages = TRUE,
        init_r = 2,
        init = inits_list,
        refresh = max(iter/50, 1)
      )
      
      #save fit object
      saveRDS(fit_lio, file=file.path("fitting", run, paste0("fit_LOO_", i, ".rds")))

      # replace pointwise log_lik if exactLOO==1 (and LOO==1)
      if (exactLOO==1){
        # extract log likelihood for data point i
        lio<- loo::extract_log_lik(fit_lio, parameter_name = "lio", merge_chains = FALSE)
        # replace by exact_elpd for left out data point i CHECK IF THIS REALLY CHANGES THE OUTPUT OF loo_
        loo_$pointwise[i,1] <- log(mean(exp(lio)))
        
        # TODO: also replace mcse_elpd_loo and other loo outputs!!!
        
        # compute loo output based on replaced pointwise log expected log prob density
        loo_$estimates["elpd_loo", "Estimate"] <- sum(loo_$pointwise[,"elpd_loo"])
        loo_$elpd_loo <- sum(loo_$pointwise[,"elpd_loo"])
      }
    
      # delete objects created during the loop
      rm(list = setdiff( setdiff(ls(), "objects_before_loop"), objects_before_loop))
    }
  }


  ## saving loo object
  if (LOO==1){
    saveRDS(list(loo_=loo_, I=I), file=file.path("fitting", run, "loo.rds"))
  }

}
