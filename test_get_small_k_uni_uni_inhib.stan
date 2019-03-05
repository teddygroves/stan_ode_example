functions{
#include get_small_k_uni_uni_inhib.stan
}
data {
  vector[4] enzyme_abundance_raw;
  vector[4] reversibility_proportional;
  real gibbs;
  real reaction_flux;
  real elementary_flux_modifier;
}
generated quantities{
  vector[8] small_k = get_small_k_uni_uni_inhib(enzyme_abundance_raw,
                                                reversibility_proportional,
                                                gibbs,
                                                reaction_flux,
                                                elementary_flux_modifier);
}
     
