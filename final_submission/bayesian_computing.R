# Load required libraries
library(MASS) # for multivariate normal functions
library(coda) # for mcmc objects
library(mcmcplots) # for mcmcplot
library(truncnorm) # for truncated normal functions


# get the obs count per ID
obs <- data %>%
    group_by(id) %>%
    summarise(n = n()) %>%
    arrange(desc(n)) 

obs

cols_to_remove <- c("id","date","sleep_duration_hr","sleep_binary")
model_data <- data[,!(names(data) %in% cols_to_remove)]
dim(model_data)
head(model_data)

# use our data
set.seed(123)
IDs <- unique(data$id)
N <- nrow(model_data)
nGroups <- length(IDs)
p <- ncol(model_data)
X <- cbind(1, as.matrix(model_data))

# design matrix D that maps the U random intercepts to the observations for each individual
# D is a block diagonal matrix with nGroups blocks where each block has a vector of 1s of length n_i
# for each individual i generate a matrix where the column i is a vector of 1s of length n_i
# first creatre a matrix of 0s of size N x nGroups
D <- matrix(0, N, nGroups)
for (i in 1:nGroups) {
  # column i is a vector of 1s of length n_i
  # n_i is the number of observations for individual i
  # determined using obs
  n_i <- sum(data$id == IDs[i])
  D[which(data$id == IDs[i]), i] <- rep(1, n_i)
}

y1 <- data$TotalMinutesAsleep # continuous outcome
y <- data$sleep_binary # binary outcome

# Number of iterations and burn-in
niter <- 20000
burnin <- 10000

# Initialize arrays to store samples
beta_samples <- matrix(NA, nrow = niter, ncol = p+1)
u_samples <- matrix(NA, nrow = niter, ncol = nGroups)
tau_samples <- rep(NA, niter)
z_samples <- matrix(NA, nrow = niter, ncol = N)

# Initialize starting values, which are matrices of 0s
beta <- matrix(0, p, 1)
u <- matrix(0, nGroups, 1)
z <- matrix(0, N, 1)
tau <- 1

# loop over the number of iterations
for (i in 1:niter) {
  # Sample beta using MVN
  beta_mu <- solve(t(X) %*% X) %*% t(X) %*% (z - D %*% u)
  beta_sigma <- solve(t(X) %*% X)
  beta <- mvrnorm(n = 1, mu = beta_mu, Sigma = beta_sigma)

  # Sample U using MVN
  I <- diag(nGroups)
  u_mu <- solve(t(D) %*% D + tau * I) %*% t(D) %*% (z - X %*% beta)
  u_sigma <- solve(t(D) %*% D + tau * I)
  u <- mvrnorm(n = 1, mu = u_mu, Sigma = u_sigma)

  # Sample tau using gamma
  tau_shape <- nGroups/2
  tau_rate <- 1/2 * sqrt(sum(u^2))
  tau <- rgamma(n = 1, shape = tau_shape, rate = tau_rate)

  # sample z using truncated normal
  z_sd = 1
  z_mean = X %*% beta + D %*% u 
  # TODO - how to condition on y? 
  z <- if_else(y == 1, 
    rtruncnorm(n = 1, a=0, b=Inf, mean = z_mean, sd = 1),
    rtruncnorm(n = 1, a=-Inf, b=0, mean = z_mean, sd = 1))

  # Store samples
  beta_samples[i, ] <- beta
  u_samples[i, ] <- u
  z_samples[i, ] <- z
  tau_samples[i] <- tau
}

# Combine samples into an mcmc object
tau_samples <- matrix(tau_samples)
# first name the columns of each
colnames(beta_samples) <- paste0("beta_", 1:(p+1))
colnames(u_samples) <- paste0("u_", 1:nGroups)
colnames(tau_samples) <- "tau"
colnames(z_samples) <- paste0("z_", 1:N)
samples <- cbind(beta_samples, u_samples, tau_samples, z_samples)
samples_mcmc <- as.mcmc(samples[-(1:burnin), ]) # Remove burn-in samples

# Visualize convergence using mcmcplot for the betas and the u's
mcmcplot(samples_mcmc[, c(colnames(beta_samples),colnames(u_samples), colnames(tau_samples))])

# Visualize posterior distributions
# change the bins to 100
# save plot to PNG
# remove the z vars
vars <- colnames(samples_mcmc)
vars <- vars[!grepl("z", vars)]
for (var in vars) {
  print(var)
  png(paste0("plots/",var,"_posterior_distributions.png"), width = 800, height = 800, res = 300)
  hist(
    samples_mcmc[, var], 
    main = var,
    xlab = var,
    breaks = 100)

  # add vertical line at 0, the mean, and the 2.5% and 97.5% quantiles
  abline(v = 0, col = "red")
  abline(v = mean(samples_mcmc[, var]), col = "blue")
  # add box for 2.5 to 97.5 quantiles
  box <- quantile(samples_mcmc[, var], probs = c(0.025, 0.975))
  rect(box[1], 0, box[2], 1000, col = "blue", density = 20, angle = 45, border = NA)
  dev.off()
}

# get probability that beta_2 is greater than 0
mean(samples_mcmc[, var] > 0)

beta_df <- as.data.frame(beta_samples)
# get the 2.5%, 50%, and 97.5% quantiles, and the mean
beta_quantiles <- apply(beta_df, 2, quantile, probs = c(0.025, 0.5, 0.975))
beta_quantiles <- rbind(beta_quantiles, colMeans(beta_df))
# name the rows
rownames(beta_quantiles) <- c("2.5%", "50%", "97.5%", "mean")
beta_quantiles

u_df <- as.data.frame(u_samples)
# get the 2.5%, 50%, and 97.5% quantiles, and the mean
u_quantiles <- apply(u_df, 2, quantile, probs = c(0.025, 0.5, 0.975))
u_quantiles <- rbind(u_quantiles, colMeans(u_df))
# name the rows
rownames(u_quantiles) <- c("2.5%", "50%", "97.5%", "mean")
u_quantiles

# Convergence of the algorithm aside from inspection. 

# Geweke diagnostic for the betas and the u's
geweke_dig_data <- geweke.diag(samples_mcmc[, c(colnames(beta_samples),colnames(u_samples), colnames(tau_samples))])
# get counts of cases where absolute value is larger than 2.58
sum(abs(geweke_dig_data) > 2.58)
