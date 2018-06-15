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
  vector[N] yield;
  vector[N] y;
}
transformed data {
  vector[N] fert_std;
  vector[N] maom_std;
  vector[N] pom_std;
  vector[N] pH_std;
  vector[N] nf_std;
  vector[N] mid_std;
  
  // Standardize variables
  fert_std = (fert - mean(fert)) / (2*sd(fert));
  maom_std = (maom - mean(maom)) / (2*sd(maom));
  pom_std = (pom - mean(pom)) / (2*sd(pom));
  pH_std = (pH - mean(pH)) / (2*sd(pH));
  nf_std = nf - mean(nf);
  mid_std = mid - mean(mid);
}
parameters  {
  // Define parameters
  real alpha_std;
  vector[K] beta;
  real<lower=0> sigma_std;
}
model {
  // Uninformative priors
  alpha_std ~ normal(0,10);
  beta ~ normal(0,10);
  sigma_std ~ cauchy(0,5);

  // Run model
  y ~ lognormal(fert_std*beta[1] + maom_std*beta[2] + pom_std*beta[3] + pH_std*beta[4] + nf_std*beta[5] + mid_std*beta[6] + alpha_std,sigma_std);
}
generated quantities {
  // Estimate nutrition impacts
  // Convert parameter from g per kg to %
  real maom_coeff = beta[2] * sd(y) / sd(maom);
  real pom_coeff = beta[3] * sd(y) / sd(pom);
  real norm_maom = maom_coeff / 10;
  real norm_pom = pom_coeff / 10;
  
  // Calculated people potentionally nourished per hectare
  real pro_rda = 28.622;  // g / d
  real zn_rda = 6.232;    // mg / d
  real fe_rda = 32.502;   // mg / d
  
  // Multiply by 10 not 1000 b/c protein is already in XX format, not 0.XX
  vector[N] pro_nourished = ((yield * 10 * (norm_maom + norm_pom)) / 365) / pro_rda;
  // No multiplication needed b/c micronutrients in same units as RDA
  vector[N] zn_nourished = ((yield * (norm_maom + norm_pom)) / 365) / zn_rda;
  vector[N] fe_nourished = ((yield * (norm_maom + norm_pom)) / 365) / fe_rda;

  // Generating estimates for SOM differences between HG and WF
  // Multiply by 10 not 1000 b/c protein is already in XX format, not 0.XX
  vector[N] pro_nour_hg = ((yield * 1.65* 10 * (norm_maom + norm_pom)) / 365) / pro_rda;
  // No multiplication needed b/c micronutrients in same units as RDA
  vector[N] zn_nour_hg = ((yield * 1.65 * (norm_maom + norm_pom)) / 365) / zn_rda;
  vector[N] fe_nour_hg = ((yield * 1.65 * (norm_maom + norm_pom)) / 365) / fe_rda;
}
