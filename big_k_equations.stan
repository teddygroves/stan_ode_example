real[] flux_equations_big_k(vector X, vector theta, real[] x_r){
  real E[4] = x_r[1:4];
  real S = x_r[5];
  real P = x_r[6];
  real Kcat[4] = to_array_1d(theta[1:4]);
  real Keq[4] = to_array_1d(theta[5:8]);
  real Km[4] = to_array_1d(theta[9:12]);
  real KI = theta[13];
  real fluxes[4];
  fluxes[1] = E[1] * Kcat[1] * (S - X[1] / Keq[1]) / (S + Km[1] + X[1] / Keq[1]);
  fluxes[2] = E[2] * Kcat[2] * (X[1] - X[2] / Keq[2]) / (X[1] + Km[2] + X[2] / Keq[2]) / (1 + X[3] / KI);
  fluxes[3] = E[3] * Kcat[3] * (X[2] - X[3] / Keq[3]) / (X[2] + Km[3] + X[3] / Keq[3]);
  fluxes[4] = E[4] * Kcat[4] * (X[3] - P / Keq[4]) / (X[3] + Km[4] + P    / Keq[4]);
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
