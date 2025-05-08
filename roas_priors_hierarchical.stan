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
  real<lower=0> sigma;
  real<lower=0, upper=1> alpha; // Adstock decay
  real<lower=0> kappa; // Half-saturation point
  real<lower=0> slope; // Curvature
  real<lower=0> roas; // Full window average iROAS
  real<lower=0> eta; // Variable needed for reparameterization
  vector[G] Z; // Variable needed for reparameterization
  real intercept_pop;
  real<lower=0> intercept_pop_sd;
  real gamma_pop;
  real<lower=0> gamma_pop_sd;
}

transformed parameters {
  matrix<lower=0>[N, G] x_adstocked;
  matrix<lower=0>[N, G] x_adstocked_saturated;
  real<lower=0> adstock_numerator;
  real<lower=0> adstock_denominator;
  real<lower=0> beta = 0;
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

  beta = (sum(X .* roas) - eta*sum(x_adstocked_saturated * Z)) / sum(x_adstocked_saturated); 
  
  # Linear model
  for (g in 1:G) {
    for (n in 1:N) {
      mu[n, g] = intercept[g] + beta*x_adstocked_saturated[n, g] + gamma[g]*C[n, g];
    }
  }
  
}

model {
  // population-level priors
  intercept_pop ~ normal(0, 10); // Uninformative prior from google MMM paper 2017
  intercept_pop_sd ~ normal(0, 1); // Uninformative prior; constrained to positive in parameters block
  gamma_pop ~ normal(0, 1); // Uninformative prior from google MMM paper 2017
  gamma_pop_sd ~ normal(0, 1); // Uninformative prior; constrained to positive in parameters block
  roas ~ lognormal(0.4408287, 0.041068); // Informative prior from 100 simulations
  eta ~ normal(0, 1);
  sigma ~ normal(0, 5);
  
  
  // Adstock and Hill parameters are assumed to be the same across geos; from google MMM paper 2017
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
  
  // geo-level variation
  for (g in 1:G) {
    Z[g] ~ normal(0,1); // Necessary as a result of reparameterizing beta in terms of ROAS
  }

  // final model
  for (g in 1:G) {
    for (n in 1:N) {
      Y[n, g] ~ normal(mu[n, g], sigma);
    }
  }

}