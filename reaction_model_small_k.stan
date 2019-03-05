functions {
#include small_k_equations.stan
}
data {
  // experimental data
  int<lower=1> N_metabolite;
  int<lower=1> N_controlled;
  int<lower=1> N_experiment;
  real measured_metabolite[N_experiment, N_metabolite];
  real controlled_concentration[N_experiment, N_controlled];
  vector[N_experiment] measured_flux;
  // hardcoded priors
  real k_loc;
  real<lower=0> k_scale;
  vector<lower=0>[N_metabolite] sigma_metabolite;
  real<lower=0> sigma_flux;
  // model config
  real rel_tol;
  real f_tol;
  int max_steps;
  int<lower=0,upper=1> LIKELIHOOD;
}
parameters {
  vector<lower=0>[26] k;
}
transformed parameters {
  matrix[N_experiment, N_metabolite] metabolite_hat; 
  real flux_hat[N_experiment, 4]; 
  for (e in 1:N_experiment){
    int x_i[0];
    metabolite_hat[e] = algebra_solver(steady_state_equations_small_k,
                                       rep_vector(0.1, 3),
                                       k,
                                       controlled_concentration[e],
                                       x_i,
                                       rel_tol, f_tol, max_steps)';
    flux_hat[e] = flux_equations_small_k(metabolite_hat[e]', k, controlled_concentration[e]);
  }
}
model {
  // priors
  k ~ lognormal(k_loc, k_scale);
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
