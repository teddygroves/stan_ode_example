vector get_small_k_uni_uni(vector enzyme_abundance_raw,
                           vector reversibility_proportional,
                           real gibbs,
                           real reaction_flux){
  /* 
     Get elementary rate constants for a reaction with the following mechanism:
     1. E + A <-> EA
     2. EA <-> EP
     3. EP <-> E + P
     
     Method is mostly copied from the GRASP file `calculateKineticParams.m`,
     however the transformation from `reversibility_proportional` to
     `reversibility` is copied from `sampleGeneralReversibilities.m`.
  */

  real RT = 2.4788;
  int enzyme_index[6] = {1, 2, 2, 3, 3, 1};
  int irrev_ix[3] = {1, 3, 5};
  int rev_ix[3] = {2, 4, 6};

  vector[3] reversibility = exp(reversibility_proportional * (gibbs / RT));
  vector[6] enzyme_abundance = enzyme_abundance_raw[enzyme_index];

  vector[6] elementary_flux;
  elementary_flux[irrev_ix] = inv(1-reversibility);
  elementary_flux[rev_ix] = reversibility .* inv(1-reversibility);

  return reaction_flux * elementary_flux .* inv(enzyme_abundance);
}
