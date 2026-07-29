# extract parameters, cross with left out data (intervention arm of group_number k),
# compute EHT mortality probability and compute log probability mass of left our data points
draws_k_spread <- spread_draws(fit_k, p_b_control, sigma_v, mu_d[group_number], sigma_d[group_number], p_h_control[group_number], mu_x, sigma_x, ndraws = 2000, seed = 27) |>
  filter(group_number == k) |>
  left_join(filter(H_all, group_number == k, !control), by = "group_number", relationship = "many-to-many") |>
  mutate(
    # Note: In contrast to the stan file, the parameters like p_h_control are already the ones for the specific group_number (no '[T_h]' needed)
    prob_D_h = p_h_control + ( 1 - p_h_control ) * pnorm( (mu_x - mu_d ) / sqrt(sigma_x^2 + sigma_d^2 )),
    logprob_D_h = log(prob_D_h),
    log_prob_D_h_2 = log1p((1 - p_h_control) * pnorm((mu_x - mu_d ) / sqrt(sigma_x^2 + sigma_d^2 )) - 1),
    log_1minusprob_D_h_2 = log1p(-p_h_control - (1 - p_h_control) * pnorm((mu_x - mu_d ) / sqrt(sigma_x^2 + sigma_d^2 ))),
    #logit_prob_D_h = p_h_control + ( 1 - p_h_control ) * pnorm( (mu_x - mu_d ) / sqrt(sigma_x^2 + sigma_d^2 )),
    l = dbinom(x = D_h, size = N_h, prob = prob_D_h, log = TRUE),
    l2 = lchoose(N_h, D_h) + D_h * log_prob_D_h_2 + (N_h - D_h) * log_1minusprob_D_h_2
    # for increased precision THIS NEEDS ANALYTIC EXPANSION SINCE prob_D_h is a sum!!!
    # If I manage to find a form with the log of the second sumand, then I can use the log version of pnorm for increased precision.
    # l = dbinom_logit(D_h, N_prob_D_h, logit_prob_D_h, log = TRUE),
  )