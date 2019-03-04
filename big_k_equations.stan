real[] flux_equations_big_k(vector X, vector theta, real[] x_r){
  real E[4] = x_r[1:4];
  real S = x_r[5];
  real P = x_r[6];
  real k_cat[4] = to_array_1d(theta[1:4]);
  real K_eq[4] = to_array_1d(theta[5:8]);
  real K_m[4] = to_array_1d(theta[9:12]);
  real K_I = theta[13];
  real fluxes[4];

  fluxes[1] = E[1] * k_cat[1] * (S - X[1] / K_eq[1]) * inv(S + K_m[1] + X[1] / K_eq[1]);
  fluxes[2] = E[2] * k_cat[2] * (X[1] - X[2] / K_eq[2]) * inv(X[1] + K_m[2] + X[2] / K_eq[2]) * inv(1 / (1 + X[3] / K_I));
  fluxes[3] = E[3] * k_cat[3] * (X[2] - X[3] / K_eq[3]) * inv(X[2] + K_m[3] + X[3] / K_eq[3]);
  fluxes[4] = E[4] * k_cat[4] * (X[3] - P / K_eq[4]) * inv(X[3] + K_m[4] + P    / K_eq[4]);
  return fluxes;
}
vector steady_state_equations_big_k(vector y, vector theta, real[] x_r, int[] x_i){
  real flux[4] = flux_equations_big_k(y, theta, x_r);
  vector[3] roots;

  roots[1] = flux[1] - flux[2];
  roots[2] = flux[2] - flux[3];
  roots[3] = flux[3] - flux[4];
  return roots;
}
