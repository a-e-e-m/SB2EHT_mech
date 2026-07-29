# extract parameters, cross with left out data (intervention arm of group_number k),
# compute EHT mortality probability and compute log probability mass of left our data points
draws_k_spread <- spread_draws(fit_k, p_b_control, sigma_v, mu_d[group_number], sigma_d[group_number], p_h_control[group_number], mu_x, sigma_x, ndraws = 2000, seed = 27) |>
  filter(group_number == models2run_LOGO_run$LOGO_ID[k]) |>
  left_join(filter(H_all, group_number == models2run_LOGO_run$LOGO_ID[k], !control), by = "group_number", relationship = "many-to-many") |>
  mutate(
    # Note: In contrast to the stan file, the parameters like p_h_control are already the ones for the specific group_number (no '[T_h]' needed)
    prob_D_h = p_h_control + ( 1 - p_h_control ) * pnorm( (mu_x - mu_d ) / sqrt(sigma_x^2 + sigma_d^2 )),
  )