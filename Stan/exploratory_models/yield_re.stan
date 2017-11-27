data {
  // Dimensions
  int<lower=0> N;
  int<lower=1> J;
  // Variables
  vector[N] trd;
  vector[N] y;
  int location[N];
}
parameters  {
  // Random intercept and slope parameters
  vector[J] a;
  vector[J] b_trd;
  real mu_a;
  real mu_b_trd;
  // Variance parameters
  real<lower=0> sigma;
  real<lower=0> sigma_a;
  real<lower=0> sigma_b_trd;
}
transformed parameters {
  vector[N] y_hat;
  
  for(n in 1:N)
    y_hat[n] = a[location[n]] + b_trd[location[n]]*trd[n];
}
model {
  mu_a ~ normal(0,100);
  mu_b_trd ~ normal(0,100);
  
  a ~ normal(mu_a, sigma_a);
  b_trd ~ normal(mu_b_trd, sigma_b_trd);
  y ~ lognormal(y_hat,sigma);
}
generated quantities {
  vector[N] y_tilde;
  
  for (n in 1:N)
    y_tilde[n] = lognormal_rng(a[location[n]] + b_trd[location[n]]*trd[n],sigma);
}
