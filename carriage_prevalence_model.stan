data {
  int<lower=1> N_countries;  // Number of countries
  int<lower=1> N_serotypes;  // Number of serotypes
  int<lower=0> N[N_countries];  // Total samples per country
  int<lower=0> x[N_countries];  // Positive samples per country
  int<lower=0> y[N_countries, N_serotypes];  // Serotype counts for positive samples
}

parameters {
  real<lower=0, upper=1> theta[N_countries];  // Prevalence for each country
  simplex[N_serotypes] phi[N_countries];  // Serotype probabilities for each country
}

model {
  // Prior for prevalence
  theta ~ beta(1, 1);  // Uniform prior
  
  // Prior for serotype probabilities
  for (i in 1:N_countries) {
    phi[i] ~ dirichlet(rep_vector(1, N_serotypes));  // Uniform Dirichlet prior
  }
  
  // Likelihood for prevalence
  for (i in 1:N_countries) {
    x[i] ~ binomial(N[i], theta[i]);
  }
  
  // Likelihood for serotype distribution
  for (i in 1:N_countries) {
    y[i] ~ multinomial(phi[i]);
  }
}

generated quantities {
  int x_pred[N_countries];
  int y_pred[N_countries, N_serotypes];
  
  for (i in 1:N_countries) {
    x_pred[i] = binomial_rng(N[i], theta[i]);
    y_pred[i] = multinomial_rng(phi[i], x_pred[i]);
  }
}