data{
  // Data set
  int nT; // number of groups
  int nT_bint; // number of groups with intensity dose bioassay data
  
  // fixed parameters
  real<lower=0> sigma_SB;
  real<lower=0, upper=1> p_b_control;

  // BA data
  int S_b; // Number of bio assay data points
  int T_b[S_b]; // group index for bio assay data points
  int treat_b[S_b]; // 1 for treatment and 0 for control
  real times_disc_dose[S_b];
  int N_b[S_b]; // Number total in bioassay
  int D_b[S_b]; // Number dead in bioassay
}

parameters{
  real mu_T[nT];
  real<lower=0> sigma_T[nT];
}

transformed parameters{
  vector[S_b] prob_D_b;

  for (i in 1:S_b){
    if (treat_b[i] == 0){
      prob_D_b[i] = p_b_control;
    } else if (treat_b[i] == 1){
      prob_D_b[i] =
        p_b_control +
        (1 - p_b_control) *
        Phi(
          (log(times_disc_dose[i]) - mu_T[T_b[i]]) /
          sqrt(sigma_SB^2 + sigma_T[T_b[i]]^2)
        );
    }
  }
}

model{
  // likelihood
  D_b ~ binomial(N_b, prob_D_b);

  // priors
  mu_T ~ normal(5, 10);
  sigma_T ~ normal(0, 5);
}
