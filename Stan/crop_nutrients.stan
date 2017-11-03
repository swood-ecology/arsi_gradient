data {
  // Dimensions
  int<lower=0> N;
  int<lower=0> K;
  
  // Variables
  matrix[N,K] x;
  vector[N] y;
}
parameters {
  real alpha;
  vector[K] beta;
  real<lower=0> sigma;
}
model {
  alpha ~ normal(0,10);
  beta ~ normal(0,10);
  sigma ~ cauchy(0,5);

  y ~ normal(x*beta + alpha, sigma);
}
generated quantities {
  vector[N] y_tilde;
  for (n in 1:N)
    y_tilde[n] = normal_rng(x[n]*beta + alpha, sigma);
}
