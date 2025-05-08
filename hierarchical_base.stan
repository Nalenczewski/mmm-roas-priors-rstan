functions {
  matrix Hill(matrix x, real kappa, real slope) {
    return (1/(1+(x/kappa)^(-slope)));
    }
}

data {
  int<lower=1> N; // Rows of data
  int<lower=1> G; // Columns of data
  matrix<lower=0>[N, G] X; // Raw spend
  matrix[N, G] C; // Control
  matrix[N, G] Y; // Revenue
  int<lower=0> L; // Max carryover duration
} 

parameters {
  vector[G] intercept;
  vector[G] gamma; // Control parameter
  vector<lower=0>[G] beta; // Set lower bound of 0
  real<lower=0> sigma;
  real<lower=0, upper=1> alpha; // Adstock decay
  real<lower=0> kappa; // Half-saturation point
  real<lower=0> slope; // Curvature
  
  // population level parameters for hierarchical model
  real intercept_pop;
  real<lower=0> intercept_pop_sd;
  real gamma_pop;
  real<lower=0> gamma_pop_sd;
  real<lower=0> beta_pop; // Set lower bound of 0
  real<lower=0> beta_pop_sd;
}

transformed parameters {
  matrix<lower=0>[N, G] x_adstocked;
  matrix<lower=0>[N, G] x_adstocked_saturated;
  real<lower=0> adstock_numerator;
  real<lower=0> adstock_denominator;
  matrix[N, G] mu; // Vector of the mean response
    
  # Adstock calculation
  for (g in 1:G) {  
    for (n in 1:N) {
      adstock_numerator = 0;
      adstock_denominator = 0;
      if (n <= L) {
        for (i in 0:(n-1)) {
          adstock_numerator += pow(alpha,i)*(X[n-i, g]);
          adstock_denominator += pow(alpha,i);
          }
      } else {
        for (i in 0:L) {
          adstock_numerator += pow(alpha,i)*(X[n-i, g]);
          adstock_denominator += pow(alpha,i);
          }
      }
      x_adstocked[n, g] = adstock_numerator / adstock_denominator;
    }
  }
    
  # Apply hill transformation
  x_adstocked_saturated = Hill(x_adstocked, kappa, slope);

  # Linear model
  for (g in 1:G) {
    for (n in 1:N) {
      mu[n, g] = intercept[g] + beta[g]*x_adstocked_saturated[n, g] + gamma[g]*C[n, g];
    }
  }
  
}

model {
  // population-level priors; uninformative
  intercept_pop ~ normal(0, 10);
  intercept_pop_sd ~ normal(0, 1);
  gamma_pop ~ normal(0, 1);
  gamma_pop_sd ~ normal(0, 1);
  beta_pop ~ normal(0, 2);
  beta_pop_sd ~ normal(0, 1);
  
  
  // Adstock and Hill parameters are assumed to be the same across geos
  kappa ~ gamma(3, 1);
  slope ~ gamma(3, 1);
  alpha ~ beta(3, 3); // restricted between 0 and 1
  
  // geo-level priors
  // geo-level intercepts
  for (g in 1:G) {
    intercept[g] ~ normal(intercept_pop, intercept_pop_sd); // Provides hierarchical structure
  }
  
  // geo-level control variable
  for (g in 1:G) {
    gamma[g] ~ normal(gamma_pop, gamma_pop_sd); // Provides hierarchical structure
  }
  
  // geo-level beta variable
  for (g in 1:G) {
    beta[g] ~ normal(beta_pop, beta_pop_sd); // Provides hierarchical structure
  }

  // final model
  for (g in 1:G) {
    for (n in 1:N) {
      Y[n, g] ~ normal(mu[n, g], sigma);
    }
  }

}