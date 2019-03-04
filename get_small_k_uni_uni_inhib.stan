vector get_small_k_uni_uni_inhib(vector enzyme_abundance_raw,
                                 vector reversibility_proportional,
                                 real gibbs,
                                 real reaction_flux,
                                 real elementary_flux_modifier){
  /* 
     Get elementary rate constants for a reaction with the following mechanism:
     1. E + A <-> EA
     2. E + I <-> EI
     3. EA <-> EP
     4. EP <-> E + P

     Method is mostly copied from the GRASP file `calculateKineticParams.m`,
     however the transformation from `reversibility_proportional` to
     `reversibility` is copied from `sampleGeneralReversibilities.m`.
  */

  real RT = 2.4788;
  int enzyme_index[8] = {1, 2, 2, 3, 2, 4, 4, 1};
  int irrev_ix[4] = {1, 3, 5, 7};
  int rev_ix[4] = {2, 4, 6, 8};
  int modified_ix[2] = {3, 4};

  vector[4] reversibility = exp(reversibility_proportional * (gibbs / RT));
  vector[8] enzyme_abundance = enzyme_abundance_raw[enzyme_index];

  vector[8] elementary_flux;
  elementary_flux[irrev_ix] = inv(1-reversibility);
  elementary_flux[rev_ix] = reversibility .* inv(1-reversibility);
  elementary_flux[modified_ix] = rep_vector(elementary_flux_modifier, 2);

  return reaction_flux * elementary_flux .* inv(enzyme_abundance);
}
