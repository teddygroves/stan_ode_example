functions {
#include big_k_equations.stan
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
  real KI_mean;
  real KI_sd;
  // model config
  real rel_tol;
  real f_tol;
  int max_steps;
  int<lower=0,upper=1> LIKELIHOOD;
}
parameters {
  vector<lower=0>[4] Kcat;
  vector<lower=0>[4] Keq;
  vector<lower=0>[4] Km;
  real<lower=0> KI;
}
transformed parameters {
  matrix[N_experiment, N_metabolite] metabolite_hat; 
  real flux_hat[N_experiment, 4]; 
  for (e in 1:N_experiment){
    vector[13] theta = append_row(append_row(append_row(Kcat, Keq), Km), KI);
    int x_i[0];
    metabolite_hat[e] = algebra_solver(steady_state_equations_big_k,
                                       rep_vector(0.1, 3),
                                       theta,
                                       controlled_concentration[e],
                                       x_i,
                                       rel_tol, f_tol, max_steps)';
    flux_hat[e] = flux_equations_big_k(metabolite_hat[e]', theta, controlled_concentration[e]);
  }
}
model {
  // priors
  Kcat ~ lognormal(enzyme_parameter_mean[1], enzyme_parameter_sd[1]);
  Keq ~ lognormal(enzyme_parameter_mean[2], enzyme_parameter_sd[2]);
  Km ~ lognormal(enzyme_parameter_mean[3], enzyme_parameter_sd[3]);
  KI ~ lognormal(KI_mean, KI_sd);
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
