# script to generate lists for controlling stan sampler

list_name <- "adapt_delta95"

control <- list(
  adapt_delta = 0.95,
  max_treedepth = 12
)

saveRDS(control, file = file.path("src", "R", "stan_control_lists", paste0(list_name, ".rds")))


list_name <- "adapt_delta95_plus"

control <- list(
  adapt_delta = 0.95,
  max_treedepth = 13
)

saveRDS(control, file = file.path("src", "R", "stan_control_lists", paste0(list_name, ".rds")))



list_name <- "adapt_delta98_plus"

control <- list(
  adapt_delta = 0.98,
  max_treedepth = 13
)

saveRDS(control, file = file.path("src", "R", "stan_control_lists", paste0(list_name, ".rds")))
