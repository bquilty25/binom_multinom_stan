# Load required libraries
library(rstan)
library(tidyverse)
library(gtools)  # For rdirichlet function

# Set seed for reproducibility
set.seed(123)

# Simulation parameters
N_countries <- 10
N_serotypes <- 5
min_samples <- 100
max_samples <- 1000

# Simulate data
simulate_data <- function() {
  N <- sample(min_samples:max_samples, N_countries, replace = TRUE)
  theta <- rbeta(N_countries, 2, 5)  # True prevalence for each country
  x <- rbinom(N_countries, N, theta)  # Positive samples
  
  phi <- matrix(0, nrow = N_countries, ncol = N_serotypes)
  y <- matrix(0, nrow = N_countries, ncol = N_serotypes)
  
  for (i in 1:N_countries) {
    phi[i,] <- rdirichlet(1, alpha = rep(1, N_serotypes))
    y[i,] <- rmultinom(1, x[i], phi[i,])
  }
  
  list(
    N_countries = N_countries,
    N_serotypes = N_serotypes,
    N = N,
    x = x,
    y = y,
    true_theta = theta,
    true_phi = phi
  )
}

# Simulate the data
sim_data <- simulate_data()

# Prepare data for Stan
stan_data <- list(
  N_countries = sim_data$N_countries,
  N_serotypes = sim_data$N_serotypes,
  N = sim_data$N,
  x = sim_data$x,
  y = sim_data$y
)

# Compile the Stan model
stan_model <- stan_model(file = "carriage_prevalence_model.stan")

# Run the Stan model
fit <- sampling(stan_model, data = stan_data, 
                chains = 4, iter = 2000, warmup = 1000, 
                cores = parallel::detectCores())

# Extract summary statistics
summary_stats <- summary(fit)$summary

# Prepare prevalence data with uncertainty
prevalence_data <- tibble(
  Country = 1:N_countries,
  True_Prevalence = sim_data$true_theta,
  Estimated_Prevalence = summary_stats[paste0("theta[", 1:N_countries, "]"), "mean"],
  Lower_CI = summary_stats[paste0("theta[", 1:N_countries, "]"), "2.5%"],
  Upper_CI = summary_stats[paste0("theta[", 1:N_countries, "]"), "97.5%"]
)

# Plot prevalence with uncertainty
ggplot(prevalence_data, aes(x = True_Prevalence, y = Estimated_Prevalence)) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  geom_pointrange(aes(ymin = Lower_CI, ymax = Upper_CI)) +
  labs(title = "True vs Estimated Prevalence with 95% CI",
       x = "True Prevalence",
       y = "Estimated Prevalence") +
  theme_minimal() +
  coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1))

ggsave("prevalence_with_uncertainty.png", width = 10, height = 8)

# Prepare serotype data with uncertainty
serotype_data <- tibble(
  Country = rep(1:N_countries, each = N_serotypes),
  Serotype = rep(paste0("Serotype_", 1:N_serotypes), N_countries),
  True_Proportion = as.vector(t(sim_data$true_phi)),
  Estimated_Proportion = summary_stats[paste0("phi[", rep(1:N_countries, each = N_serotypes), ",", rep(1:N_serotypes, N_countries), "]"), "mean"],
  Lower_CI = summary_stats[paste0("phi[", rep(1:N_countries, each = N_serotypes), ",", rep(1:N_serotypes, N_countries), "]"), "2.5%"],
  Upper_CI = summary_stats[paste0("phi[", rep(1:N_countries, each = N_serotypes), ",", rep(1:N_serotypes, N_countries), "]"), "97.5%"]
)

# Plot serotype distribution with uncertainty
ggplot(serotype_data, aes(x = Serotype, y = Estimated_Proportion, color = Serotype)) +
  geom_pointrange(aes(ymin = Lower_CI, ymax = Upper_CI), position = position_dodge(width = 0.5)) +
  geom_point(aes(y = True_Proportion), color = "red", size = 2, position = position_dodge(width = 0.5)) +
  facet_wrap(~Country, scales = "free_y", ncol = 2) +
  labs(title = "Serotype Distribution by Country with 95% CI",
       subtitle = "Points: Estimated Proportion, Red Dots: True Proportion",
       x = "Serotype",
       y = "Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  coord_cartesian(ylim = c(0, 1))

ggsave("serotype_distribution_with_uncertainty.png", width = 12, height = 16)

print(prevalence_data)
print(serotype_data)