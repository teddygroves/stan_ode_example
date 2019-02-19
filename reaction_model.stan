functions {
#include functions.stan
}
data {
  // experimental data
  int<lower=1> N_metabolite;
  int<lower=1> N_controlled;
  int<lower=1> N_experiment;
  real measured_metabolite[N_experiment, N_metabolite];
  real controlled_concentration[N_experiment, N_controlled];
  real measured_flux[N_experiment];
  // hardcoded priors
  vector<lower=0>[N_metabolite] sigma_metabolite;
  real<lower=0> sigma_flux;
  vector[3] enzyme_parameter_mean;
  vector<lower=0>[3] enzyme_parameter_sd;
  real K_I_mean;
  real K_I_sd;
  // model config
  real rel_tol;
  real f_tol;
  int max_steps;
  int<lower=0,upper=1> LIKELIHOOD;
}
parameters {
  vector<lower=0>[4] k_cat;
  vector<lower=0>[4] K_eq;
  vector<lower=0>[4] K_m;
  real<lower=0> K_I;
}
transformed parameters {
  matrix[N_experiment, N_metabolite] metabolite_hat; 
  real flux_hat[N_experiment, 4]; 
  for (e in 1:N_experiment){
    vector[13] theta = append_row(append_row(append_row(k_cat, K_eq), K_m), K_I);
    int x_i[0];
    metabolite_hat[e] = algebra_solver(reaction_steady_state_system, // equation system 
                                       rep_vector(0.1, 3),           // initial guess
                                       theta,                        // parameters
                                       controlled_concentration[e],  // real-valued data
                                       x_i,                          // integer-valued data
                                       rel_tol, f_tol, max_steps)';  // control parameters
    flux_hat[e] = flux_equations(metabolite_hat[e]', theta, controlled_concentration[e]);
  }
}
model {
  // priors
  k_cat ~ lognormal(enzyme_parameter_mean[1], enzyme_parameter_sd[1]);
  K_eq ~ lognormal(enzyme_parameter_mean[2], enzyme_parameter_sd[2]);
  K_m ~ lognormal(enzyme_parameter_mean[3], enzyme_parameter_sd[3]);
  K_I ~ lognormal(K_I_mean, K_I_sd);
  // measurement model
  if (LIKELIHOOD == 1){
    for (e in 1:N_experiment){
      measured_flux[e] ~ normal(flux_hat[e], sigma_flux);
      measured_metabolite[e] ~ normal(metabolite_hat[e], sigma_metabolite);
    }
  }
}
generated quantities {
  real metabolite_pred[N_experiment, N_metabolite];
  real flux_pred[N_experiment, 4];
  for (e in 1:N_experiment){
    for (m in 1:N_metabolite){
      metabolite_pred[e, m] = normal_rng(metabolite_hat[e, m], sigma_metabolite[m]);
    }
    for (r in 1:4){
      flux_pred[e, r] = normal_rng(flux_hat[e, r], sigma_flux);
    }
  }
}
