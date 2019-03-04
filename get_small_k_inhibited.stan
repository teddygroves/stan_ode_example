vector get_small_k_inhibited(vector[4] enzyme_state_fraction,
                             vector[4] reversibility,
                             vector[4] reaction_flux){
  /* 
     Get elementary rate constants for a reaction with the following mechanism:
     
     1. E + A <-> EA
     2. E + I <-> EI
     3. EA <-> EP
     4. EP <-> E + P
     
     The order of the 4 enzyme states is E, EI, EA, EP. 
     
     The order of the 8 elementary constants is k1+, k2+, k3+, k4+, k1-, k2-, k3-, k4-.

     Since there are two paths through the reaction, a branching vector is
     needed. This is found by multiplying a path (reversibility?) matrix N,
     i.e.

     [[1, 0],
     [0, 1],
     [1, 0],
     [1, 0]]
      
     by a weight vector [w1, w2].
  */

  int n_arrow = 8;
  int n_enzyme_state = 4;
  vector[n_arrow] small_k;
  // which enzyme state is produced by each elementary reaction?
  int[n_arrow] enzyme_index = [3, 2, 4, 1, 1, 1, 3, 4];
  // which reaction parameter goes with which elementary reaction?
  int[n_arrow] reaction_index = [1, 2, 3, 4, 1, 2, 3, 4];

  vector[n_arrow] P = enzyme_fraction[enzyme_index];
  vector[n_arrow] gamma;
  for (i in 1:n_arrow){
    int r = reaction_index[i];
    real R = reaction_flux[r] > 0 ? reversibility[r] : inv(reversibility[r]);
    gamma[i] = i < n_reaction / 2 ? inv(1-R) : R / (1 - R);
  }
  return  gamma / P;
}
