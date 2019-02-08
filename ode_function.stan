real[] reaction_ode(real t, real[] y, real[] theta, real[] x_r, int[] x_i){

  real X_1 = y[1];
  real X_2 = y[2];
  real X_3 = y[3];
  
  real E_1 = x_r[1];
  real E_2 = x_r[2];
  real E_3 = x_r[3];
  real E_4 = x_r[4];
  real S = x_r[4];
  real P = x_r[5];

  real kcat_1 = theta[1];
  real kcat_2 = theta[2];
  real kcat_3 = theta[3];
  real kcat_4 = theta[4];

  real Keq_1 = theta[5];
  real Keq_2 = theta[6];
  real Keq_3 = theta[7];
  real Keq_4 = theta[8];

  real Km_1 = theta[9];
  real Km_2 = theta[10];
  real Km_3 = theta[11];
  real Km_4 = theta[12];

  real KI_2 = theta[13];

  real v1 = E_1 * (kcat_1 * (S-X_1 / Keq_1) / (S + Km_1 + X_1 / Keq_1));
  real v2 = E_2 * (kcat_2 * ((X_1-X_2 / Keq_2) / (X_1 + Km_2 + X_2 / Keq_2)) * (1 / (1 + X_3 / KI_2)));
  real v3 = E_3 * (kcat_3 * (X_2-X_3 / Keq_3) / (X_2 + Km_3 + X_3 / Keq_3));
  real v4 = E_4 * (kcat_4 * (X_3-P / Keq_4) / (X_3 + Km_4 + P / Keq_4));
  
  real dmetdt[3];
  dmetdt[1] = v1 - v2;   // X1
  dmetdt[2] = v2 - v3;   // X2
  dmetdt[3] = v3 - v4;   // X3
  return dmetdt;
}
