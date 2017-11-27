data {
  // Dimensions
  int<lower=0> N;
  int<lower=1> J;
  // Variables
  vector[N] pom;
  vector[N] maom;
  vector[N] fert;
  vector[N] pH;
  vector[N] y;
  int location[N];
}
parameters  {
  // Random intercept and slope parameters
  vector[J] b_pom;
  real a;
  real b_maom;
  real b_fert;
  real b_pH;
  real mu_b_pom;
  // Variance parameters
  real<lower=0> sigma;
  real<lower=0> sigma_b_pom;
}
transformed parameters {
  vector[N] y_hat;
  
  for(n in 1:N)
    y_hat[n] = a + b_maom*maom[n] + b_pom[location[n]]*pom[n] + b_fert*fert[n] + b_pH*pH[n];
}
model {
  a ~ normal(0,10);
  mu_b_pom ~ normal(0,10);
  b_pH ~ normal(0,10);
  b_fert ~ normal(0,10);
  b_maom ~ normal(0,10);
  
  b_pom ~ normal(mu_b_pom, sigma_b_pom);
  y ~ lognormal(y_hat,sigma);
}
