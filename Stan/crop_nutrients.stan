data {
  // Dimensions
  int<lower=0> N;
  int<lower=0> K;
  
  // Variables
  matrix[N,K] x;
  vector[N] y;
}
transformed data {
  // Define transformed matrix
  matrix[N,K] x_std;
  
  // Transform data matrix to standardized
  for (k in 1:K)
    x_std[,k] = (x[,k] - mean(x[,k])) / (2*sd(x[,k]));
}
parameters {
  // Define standardized parameters
  real alpha_std;
  vector[K] beta_std;
  real<lower=0> sigma_std;
}
model {
  alpha_std ~ normal(0,10);
  beta_std ~ normal(0,10);
  sigma_std ~ cauchy(0,5);

  y ~ lognormal(x_std*beta_std + alpha_std, sigma_std);
}
generated quantities {
  vector[N] y_tilde;
  vector[K] beta;
  
  for(k in 1:K)
    beta[k] = beta_std[k] * sd(y) / sd(x[,k]);
  
  for (n in 1:N)
    y_tilde[n] = lognormal_rng(x_std[n]*beta_std + alpha_std, sigma_std);
}
