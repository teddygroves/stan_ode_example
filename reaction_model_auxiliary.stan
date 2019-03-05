functions {
#include get_small_k_uni_uni.stan
#include get_small_k_uni_uni_inhib.stan
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
  vector<lower=0>[N_metabolite] sigma_metabolite;
  real<lower=0> sigma_flux;
  real dir_enz;
  real dir_rev;
  real gibbs_mean;
  real<lower=0> gibbs_sd;
  real inhibition_modifier_beta_params[2];
  // model config
  real rel_tol;
  real f_tol;
  int max_steps;
  int<lower=0,upper=1> LIKELIHOOD;
}
transformed data {
  int<lower=1> N_enzyme_state_uni_uni = 3;
  int<lower=1> N_enzyme_state_uni_uni_inhib = 4;
  int<lower=1> N_elementary_reaction_uni_uni = 3;
  int<lower=1> N_elementary_reaction_uni_uni_inhib = 4;
}
parameters {
  simplex[N_enzyme_state_uni_uni] enzyme_abundance_uni_uni[3];
  simplex[N_enzyme_state_uni_uni_inhib] enzyme_abundance_uni_uni_inhib;
  simplex[N_elementary_reaction_uni_uni] reversibility[4];
  vector<upper=0>[4] gibbs;
  vector<lower=0>[N_experiment] reaction_flux;
  real<lower=0,upper=1> inhibition_modifier;
}
transformed parameters {
  matrix[N_experiment, N_metabolite] metabolite_hat; 
  for (e in 1:N_experiment){
    int x_i[0];
    vector[26] k;
    k[1:6] = get_small_k_uni_uni(enzyme_abundance_uni_uni[1],
                                 reversibility[1],
                                 gibbs[1],
                                 reaction_flux[e]);
    k[7:14] = get_small_k_uni_uni_inhib(enzyme_abundance_uni_uni_inhib,
                                        reversibility[2],
                                        gibbs[2],
                                        reaction_flux[e],
                                        inhibition_modifier);
    k[15:20] = get_small_k_uni_uni(enzyme_abundance_uni_uni[2],
                                   reversibility[3],
                                   gibbs[3],
                                   reaction_flux[e]);
    k[21:26] = get_small_k_uni_uni(enzyme_abundance_uni_uni[3],
                                   reversibility[4],
                                   gibbs[4],
                                   reaction_flux[e]);
    metabolite_hat[e] = algebra_solver(steady_state_equations_small_k,
                                       rep_vector(0.1, 3),
                                       k,
                                       controlled_concentration[e],
                                       x_i,
                                       rel_tol, f_tol, max_steps)';
  }
}
model {
  // priors
  gibbs ~ normal(gibbs_mean, gibbs_sd);
  inhibition_modifier ~ beta(inhibition_modifier_beta_params[1], inhibition_modifier_beta_params[2]);
  for (r in 1:3){
    enzyme_abundance_uni_uni[r] ~ dirichlet(rep_vector(dir_enz, N_enzyme_state_uni_uni));
    reversibility[r] ~ dirichlet(rep_vector(dir_rev, N_elementary_reaction_uni_uni));
  }
  enzyme_abundance_uni_uni_inhib ~ dirichlet(rep_vector(dir_enz, N_enzyme_state_uni_uni_inhib));
  reaction_flux ~ normal(measured_flux, sigma_flux);
  // measurement model
  if (LIKELIHOOD == 1){
    for (e in 1:N_experiment){
      measured_metabolite[e] ~ normal(metabolite_hat[e], sigma_metabolite);
    }
  }
}
generated quantities {
  real metabolite_pred[N_experiment, N_metabolite];
  for (e in 1:N_experiment){
    for (m in 1:N_metabolite){
      metabolite_pred[e, m] = normal_rng(metabolite_hat[e, m], sigma_metabolite[m]);
    }
  }
}
