functions {
#include ode_function.stan
}
data {
  int<lower=1> N_metabolite;
  int<lower=1> N_controlled;
  int<lower=1> N_experiment;
  real initial_metabolite_concentration[N_metabolite];
  real final_metabolite_concentration[N_experiment, N_metabolite];
  real controlled_concentration[N_experiment, N_controlled];
  int<lower=0,upper=1> LIKELIHOOD;
  real<lower=0> t_0;
  real<lower=0> t_steady;
  vector<lower=0>[N_metabolite] sigma;
  vector[3] enzyme_parameter_mean;
  vector<lower=0>[3] enzyme_parameter_sd;
  real K_I_mean;
  real K_I_sd;
}
parameters {
  real<lower=0> k_cat[4];
  real<lower=0> K_eq[4];
  real<lower=0> K_m[4];
  real<lower=0> K_I;
}
transformed parameters {
  real final_metabolite_concentration_hat[N_experiment, N_metabolite]; 
  for (e in 1:N_experiment){
    real theta[13] = append_array(append_array(append_array(k_cat, K_eq), K_m), {K_I});
    int x_i[0];
    final_metabolite_concentration_hat[e] =
      tail(integrate_ode_rk45(reaction_ode,
                              initial_metabolite_concentration,
                              t_0,
                              {t_steady},
                              theta,
                              controlled_concentration[e],
                              x_i), 1)[1];
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
      final_metabolite_concentration[e] ~ normal(final_metabolite_concentration_hat[e], sigma);
    }
  }
}
generated quantities {
  real final_metabolite_concentration_pred[N_experiment, N_metabolite];
  for (e in 1:N_experiment){
    for (m in 1:N_metabolite){
      final_metabolite_concentration_pred[e, m] = normal_rng(final_metabolite_concentration_hat[e, m], sigma[m]);
    }
  }
}
