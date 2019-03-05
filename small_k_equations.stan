real[] flux_equations_small_k(vector X, vector small_k, real[] controlled_concentrations){
  real S = controlled_concentrations[5];
  real P = controlled_concentrations[6];
  real E[4] = controlled_concentrations[1:4];
  vector[6] k1 = segment(small_k, 1, 6);
  vector[8] k2 = segment(small_k, 7, 8);
  vector[6] k3 = segment(small_k, 16, 6);
  vector[6] k4 = segment(small_k, 21, 6);
  real fluxes[4];
  fluxes[1] = E[1] * (k1[1]*k1[3]*k1[5]*S - k1[2]*k1[4]*k1[6]*X[1]) /
    (k1[3]*k1[5]
     + k1[1]*(k1[4] + k1[3] + k1[5])*S
     + k1[4]*k1[6]*X[1]
     + k1[3]*k1[6]*X[1]
     + k1[2]*(k1[4] + k1[5] + k1[6]*X[1]));
  fluxes[2] = E[2] * (k2[8]*(k2[1]*k2[3]*k2[5]*X[1] - k2[2]*k2[4]*k2[6]*X[2])) /
    (k2[3]*k2[5]*k2[8]
     + k2[1]*(k2[4] + k2[3] + k2[5])*k2[8]*X[1]
     + k2[4]*k2[6]*k2[8]*X[2]
     + k2[3]*k2[6]*k2[8]*X[2]
     + k2[3]*k2[5]*k2[7]*X[3]
     + k2[2]*(k2[6]*k2[8]*X[2] + k2[4]*(k2[8] + k2[7]*X[3]) + k2[5]*(k2[8] + k2[7]*X[3])));
  fluxes[3] = E[3] * (k3[1]*k3[3]*k3[5]*X[2] - k3[2]*k3[4]*k3[6]*X[3]) /
    (k3[3]*k3[5]
     + k3[1]*(k3[4] + k3[3] + k3[5])*X[2]
     + k3[4]*k3[6]*X[3]
     + k3[3]*k3[6]*X[3]
     + k3[2]*(k3[4] + k3[5] + k3[6]*X[3]));
  fluxes[4] = E[4] * (-(k4[2]*k4[4]*k4[6]*P) + k4[1]*k4[3]*k4[5]*X[3]) /
    (k4[3]*k4[5]
     + k4[4]*k4[6]*P
     + k4[3]*k4[6]*P
     + k4[2]*(k4[4] + k4[5] + k4[6]*P) + k4[1]*(k4[4] + k4[3] + k4[5])*X[3]);

  return fluxes;
}
vector steady_state_equations_small_k(vector y, vector k, real[] x_r, int[] x_i){
  real flux[4] = flux_equations_small_k(y, k, x_r);
  vector[3] roots;

  roots[1] = flux[1] - flux[2];
  roots[2] = flux[2] - flux[3];
  roots[3] = flux[3] - flux[4];
  return roots;
}
