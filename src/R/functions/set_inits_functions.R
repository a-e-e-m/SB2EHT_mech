# Functions for defining different inits for stan sampler


set_inits_C <- function(seed) {
  
  set.seed(seed)
  
  window_small <- 0.1
  window_big <- 1

  sd_small <- 0.1
  sd_big <- 1
  
  list(
    a = rnorm(n=1, mean = 0.5, sd =sd_small),
    a2 = rnorm(n=1, mean = -1.5, sd = sd_small),
    a3 = rnorm(n=1, mean = -0.5, sd = sd_small),
    r = runif(n=1, min = 1 - window_small, max = 1 + window_small),
    kappa = runif(n=1, min = 3 - window_small, max = 3 + window_small),
    sigma_d = runif(n=1, min = 5 - window_big, max = 5 + window_big),
    sigma_c = runif(n=1, min = 5 - window_big, max = 5 + window_big),
    sigma_e = runif(n=1, min = 10 - window_big, max = 10 + window_big),
    sigma_g = runif(n=1, min = 1 - window_small, max = 1 + window_small),
    sigma_h = runif(n=1, min = 0.1 - window_small, max = 0.1 + window_small)
  )
}
