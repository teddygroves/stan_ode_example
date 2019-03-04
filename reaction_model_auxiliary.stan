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
  int<lower=1> N_enzyme_state_uninhib;
  int<lower=1> N_enzyme_state_inhib;
  int<lower=1> N_elementary_reaction_uninhib;
  int<lower=1> N_elementary_reaction_inhib;
  real measured_metabolite[N_experiment, N_metabolite];
  real controlled_concentration[N_experiment, N_controlled];
  real measured_flux;
  // hardcoded priors
  vector<lower=0>[N_metabolite] sigma_metabolite;
  real<lower=0> sigma_flux;
  // model config
  real rel_tol;
  real f_tol;
  int max_steps;
  int<lower=0,upper=1> LIKELIHOOD;
}
parameters {
  simplex[N_enzyme_state_uninhib] enzyme_abundance_uni_uni[3];
  simplex[N_enzyme_state_inhib] enzyme_abundance_uni_uni_inhib;
  simplex[N_elementary_reaction_uninhib] reversibility_uni_uni[3];
  simplex[N_elementary_reaction_inhib] reversibility_uni_uni_inhib;
  vector<upper=0>[4] gibbs;
  real reaction_flux;
}
transformed parameters {
  vector[4] Kcat;
  vector[4] Keq;
  vector[4] Km;
  real Ki;
  matrix[N_experiment, N_metabolite] metabolite_hat; 
  real flux_hat[N_experiment, 4]; 
  for (e in 1:N_experiment){
    vector[26] small_k;
    small_k[1:6] = get_small_k_uni_uni(enzyme_abundance_uni_uni[1]
                                       reversibility_uni_uni[1],
                                       gibbs[1],
                                       reaction_flux);
    small_k[7:14] = get_small_k_uni_uni_inhib(enzyme_state_proportion_uni_uni_inhib,
                                              reversibility_uni_uni_inhib,
                                              elementary_flux_uni_uni_inhib);
    small_k[15:20] = get_small_k_uni_uni(enzyme_abundance_uni_uni[2]
                                         reversibility_uni_uni[2],
                                         gibbs[2],
                                         reaction_flux);
    small_k[21:26] = get_small_k_uni_uni(enzyme_abundance_uni_uni[3]
                                         reversibility_uni_uni[3],
                                         gibbs[3],
                                         reaction_flux);
  int x_i[0];
  metabolite_hat[e] = algebra_solver(steady_state_equations_small_k, // equation system 
                                       rep_vector(0.1, 3),             // initial guess
                                       small_k,
                                       controlled_concentration[e],    // real-valued data
                                       x_i,                            // integer-valued data
                                       rel_tol, f_tol, max_steps)';    // control parameters
    flux_hat[e] = flux_equations(metabolite_hat[e]', theta, controlled_concentration[e]);
  }
}
model {
  // priors
  for (enzyme in 1:3){
    enzyme_abundance_uni_uni[enzyme] ~ dirichlet(2);
    reversibility_uni_uni[enzyme] ~ dirichlet(2);
  }
  enzyme_abundance_uni_uni_inhib ~ dirichlet(2);
  reversibility_uni_uni_inhib ~ dirichlet(2);
  reaction_flux ~ normal(measured_flux, sigma_flux);
  gibbs ~ normal(-12, 2);
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
