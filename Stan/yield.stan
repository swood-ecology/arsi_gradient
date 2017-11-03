data {
  // Dimensions
  int<lower=0> N;
  int<lower=0> K;
  
  // Variables
  vector[N] fert;
  vector[N] maom;
  vector[N] pom;
  vector[N] trd;
  vector[N] y;
}
parameters  {
  real inter;
  real alpha;
  vector[K] beta;
  real<lower=0> sigma;
}
model {
  inter ~ normal(0,10);
  alpha ~ normal(0,10);
  beta ~ normal(0,10);
  sigma ~ cauchy(0,5);
  

  y ~ skew_normal(fert*beta[1] + maom*beta[2] + pom*beta[3] + trd*beta[4] + inter, sigma, alpha);
}
generated quantities {
  vector[N] y_tilde;
  for (n in 1:N)
    y_tilde[n] = skew_normal_rng(fert[n]*beta[1] + maom[n]*beta[2] + pom[n]*beta[3] + trd[n]*beta[4] + inter, sigma, alpha);
}

