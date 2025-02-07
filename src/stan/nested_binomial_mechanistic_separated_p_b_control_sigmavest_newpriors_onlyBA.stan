
data{
    // Data set
  int nT; // number of groups
  int nT_bint; // number of groups with intensity dose bioassay data
  int T_bint_indices[nT_bint]; // all occuring group indices with intensity dose bioassay data, sorted
  int T_bdisc_indices[nT - nT_bint]; // all occuring group indices without intensity dose bioassay data, sorted
  
  // BA data
  int S_b; // Number of bio assay data points
  int T_b[S_b]; // group index for bio assay data points
  int treat_b[S_b]; // 1 for treatment and 0 for control
  int times_disc_dose[S_b]; 
  int int_dose_available[S_b]; 
  int N_b[S_b]; // Number total in bioassay
  int D_b[S_b]; // Number dead in bioassay
  
  // switches for LOO
  int<lower=0, upper=1> LOO;
  int<lower=0, upper=1> exactLOO;
  
    // for exactl LOO: left out data point
    int LOdataset[exactLOO ? 1 : 0]; // Group index of left out data set
    // here comes left out EHT treatment data (only treatment part, corresponding control data and BA data (both control and treatment) to be fed in with data above)
      // EHT data all
      int S_h_LO[exactLOO ? 1 : 0]; // Number of all left out hut data points (with full or partial breakdown)
      int N_h_LO[exactLOO ? S_h_LO[1] : 0]; // Number total in hut trial
      int D_h_LO[exactLOO ? S_h_LO[1] : 0]; // Number dead in hut trial
      
      // EHT data full breakdown
      int S_h_f_LO[exactLOO ? 1 : 0]; // Number of left out hut data points with full breakdown
      int A_LO[exactLOO ? S_h_f_LO[1] : 0]; // Number alive in hut trial
      int AF_LO[exactLOO ? S_h_f_LO[1] : 0]; // Number alivefed in hut trial
}

parameters{
  real mu_d[nT];
  real<lower=0> sigma_d[nT];
  real<lower=0> sigma_v;
  real<lower=0, upper=1> p_b_control[nT];
}

transformed parameters{
  // declarations
    vector[S_b] prob_D_b; // logit transformed mortality probability per BA data point

  // statements
    for (i in 1:S_b){
      if (treat_b[i] == 0){
        prob_D_b[i] = p_b_control[T_b[i]];
      } else if (treat_b[i] == 1){
        prob_D_b[i] = p_b_control[T_b[i]] + ( 1 - p_b_control[T_b[i]] ) * Phi( (log(times_disc_dose[i]) - mu_d[T_b[i]] ) / sqrt(sigma_v + sigma_d[T_b[i]]^2 ) ); 
      }
    }
}

model{
  // likelihood
    // likelihood for bioassay mortality
    D_b ~ binomial(N_b, prob_D_b);

  // priors
  mu_d ~ normal(5, 10);
  sigma_d ~ normal(0, 5);
  sigma_v ~ normal(0, 5);
  p_b_control ~ beta(1,10);
}
