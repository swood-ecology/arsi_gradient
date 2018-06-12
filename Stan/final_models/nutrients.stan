data {
  // Dimensions
  int<lower=0> N;
  int<lower=0> K;
  
  // Variables
  vector[N] fert;
  vector[N] maom;
  vector[N] pom;
  vector[N] pH;
  vector[N] nf;
  vector[N] mid;
  vector[N] y;
}
parameters  {
  // Define parameters
  real alpha;
  vector[K] beta;
  real<lower=0> sigma;
}
model {
  // Uninformative priors
  alpha ~ normal(0,10);
  beta ~ normal(0,10);
  sigma ~ cauchy(0,5);

  // Run model
  y ~ lognormal(fert*beta[1] + maom*beta[2] + pom*beta[3] + pH*beta[4] + nf*beta[5] + mid*beta[6] + alpha,sigma);
}
generated quantities {
  // Estimate nutrition impacts
  
  // Convert parameter from g per kg to %
  real norm_maom = beta[2] / 10;
  
  // Calculated people potentionally nourished per hectare
  vector[N] pro_nourished;
  vector[N] zn_nourished;
  vector[N] fe_nourished;
  
  for (n in 1:N){
    pro_nourished[n] = ((y[n] * 1000 * norm_maom)/365)/28.622;
    zn_nourished[n] = ((y[n] * 1000 * norm_maom)/365)/6.232;
    fe_nourished[n] = ((y[n] * 1000 * norm_maom)/365)/32.502;
  }
}
