real[] flux_equations(real[] E,
                      real[] k_cat,
                      real[] K_eq,
                      real[] K_m,
                      real KI,
                      real S,
                      real P,
                      real[] X){
  real fluxes[4];
  fluxes[1] = E[1] * (k_cat[1] * (S-X[1] / K_eq[1]) / (S + K_m[1] + X[1] / K_eq[1]));
  fluxes[2] = E[2] * (k_cat[2] * ((X[1]-X[2] / K_eq[2]) / (X[1] + K_m[2] + X[2] / K_eq[2])) * (1 / (1 + X[3] / KI)));
  fluxes[3] = E[3] * (k_cat[3] * (X[2]-X[3] / K_eq[3]) / (X[2] + K_m[3] + X[3] / K_eq[3]));
  fluxes[4] = E[4] * (k_cat[4] * (X[3]-P / K_eq[4]) / (X[3] + K_m[4] + P / K_eq[4]));
  return fluxes;
}
real[] reaction_ode(real t, real[] X, real[] theta, real[] x_r, int[] x_i){

  real E[4] = x_r[1:4];
  real S = x_r[5];
  real P = x_r[6];
  real k_cat[4] = theta[1:4];
  real K_eq[4] = theta[5:8];
  real K_m[4] = theta[9:12];
  real KI = theta[13];

  real v[4] = flux_equations(E, k_cat, K_eq, K_m, KI, S, P, X);

  real dXdt[3];
  dXdt[1] = v[1] - v[2];
  dXdt[2] = v[2] - v[3];
  dXdt[3] = v[3] - v[4];
  return dXdt;
}

