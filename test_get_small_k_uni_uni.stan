functions{
#include get_small_k_uni_uni.stan
}
data {
  vector[3] enzyme_abundance_raw;
  vector[3] reversibility_proportional;
  real gibbs;
  real reaction_flux;
}
generated quantities{
  vector[6] small_k = get_small_k_uni_uni(enzyme_abundance_raw,
                                          reversibility_proportional,
                                          gibbs,
                                          reaction_flux);
}
     
