source("~/qPCR_2/functions.R", echo=TRUE)

# simulated data -------

numCycles <- 45

# sites
n <- 100

# sampling occasions
M <- 1 + rpois(n, lambda = 3)
idx_site <- rep(1:n, M)

# technical replicates
K <- 1 + rep(3, sum(M))
idx_site_K <- rep(idx_site, K)
idx_sample_K <- rep(1:sum(M), K)

# PCR plates
n_P <- n
P <- idx_site_K

# hyperprior
{
  tau <- 1
  sigma <- 2
  sigma_y <- 1
  beta_b0 <- 0
  beta_w0 <- -1
  theta0 <- .02
  p0 <- .05
  nu_0 <- 20
  nu <- -1
  sigma_nu <- 1
  
  alpha1_0 <- 45
  alpha2_0 <- -1.7
  
  lambda <- 20
  sigma_lambda <- 10
  
  phi_0 <- -26.643
  phi_1 <- 8.3
}

# standards number and concentration
{
  K_standards <- 3
  qty <- c(3e+07, 3e+06, 3e+05, 3e+04, 3e+03, 3e+02, 3e+01)
  data_standard <- data.frame(P = rep(1:n, each = K_standards * length(qty)),
                              Qty = rep(rep(qty, each = K_standards), times = n),
                              Ct = rep(NA, n * K_standards * length(qty)))
}

# covariates
ncov_z <- 2
X_b <- matrix(rnorm(ncov_z * n), n, ncov_z)
ncov_w <- 2
X_w <- matrix(rnorm(ncov_w * sum(M)), sum(M), ncov_w)

# covariates coefficients
beta_b <- sample(c(-1,1), size = ncov_z, replace = T)# rnorm(ncov_z)
beta_w <- rnorm(ncov_w)

# amount of biomass at the site
l <- beta_b0 + X_b %*% beta_b + rnorm(n, sd = tau)
# b <- ifelse(logb > 1, logb, 0)
# z <- as.numeric(l > 0)
# b <- z * (exp(l) - 1)

# z_samples <- z[idx_site]
l_samples <- l[idx_site]
# b_samples <- l[idx_site]

# STAGE 1 -----
# w <- rep(NA, sum(M))
v <- rep(NA, sum(M)) # collected biomass
# delta <- rep(NA, sum(M))
v_tilde <- rep(NA, sum(M)) # biomass from contamination
# gamma <- rep(NA, sum(M))

# true positives 
v <- l_samples + beta_w0 + X_w %*% beta_w + rnorm(sum(M), sd = sigma)
w_hat <- (v > 0) * (exp(v) - 1)
# v[z == 1 & v < 0] <- 0

# assign to 0 sites with log-collected biomass less than 1
v_tilde[v > 0] <- 0
v_tilde[v <= 0] <- rnorm(sum(v <= 0), nu, sd = sigma_nu)
w_tilde <- (v_tilde > 0) * (exp(v_tilde) - 1)
# v[z==0] <- rnorm(sum(z==0), nu, sd = sigma_nu)
# delta[z==0] <- v[z==0] > 0
# v[z == 0 & v < 0] <- 0
w <- w_hat + w_tilde

w_pcr <- w[idx_sample_K]
# w <- ifelse(v[idx_sample_K] > 0, exp(v[idx_sample_K] - 1), 0)
# w <- ifelse(delta == 1, w_tilde, 0)

# STAGE 2 --------

alpha1 <- rnorm(n_P, mean = alpha1_0, sd = 1)
alpha2 <- rnorm(n_P, mean = alpha2_0, sd = .1)

# standards
p_imk_star <- logistic(phi_0 + phi_1 * log(data_standard$Qty))
delta_star <- sapply(p_imk_star, function(p){
  rbinom(1, 1, p)
})

delta_star <- sapply(delta_star, function(x){
  if(x == 0){
    2 * rbinom(1, 1, p0)
  } else {
    x
  }
})

data_standard$Ct <- ifelse(delta_star == 1, 
                           alpha1[data_standard$P] + alpha2[data_standard$P] * 
                             log(data_standard$Qty) + 
                             rnorm(nrow(data_standard), sd = sigma_y), 
                           ifelse(delta_star == 1, 
                                  rnorm(length(data_standard$P), lambda, sigma_lambda), 
                                  0))

# real
data_real <- data.frame(Site = idx_site_K,
                        Sample = idx_sample_K,
                        Plate = idx_site_K)

p_imk <- logistic(phi_0 + phi_1 * log(w_pcr))
delta <- sapply(p_imk, function(p){
  rbinom(1, 1, p)
})

delta <- sapply(delta, function(x){
  if(x == 0){
    2 * rbinom(1, 1, p0)
  } else {
    x
  }
})

Ct <- ifelse(delta == 1, 
             alpha1[data_real$Plate] + alpha2[data_real$Plate] * log(w_pcr) +
               rnorm(nrow(data_real), sd = sigma_y), 
             ifelse(delta == 2, 
                    rnorm(length(P), lambda, sigma_lambda), 
                    0))

# false positives
# idx_0 <- which(delta[idx_sample_K] == 0 & gamma[idx_sample_K] == 0)
# zeta <- rbinom(length(idx_0), 1, p0)
# Ct[idx_0] <- zeta * rnorm(length(idx_0), nu_0, sd = sigma_y) 

data_real$Ct <- Ct

# save true values -----------

{
  l_true <- l
  v_true <- v
  
  tau_true <- tau
  vtilde_true <- v_tilde
  
  alpha1_true <- alpha1
  alpha2_true <- alpha2
  
  beta_0b_true <- beta_b0
  beta_b_true <- beta_b
  
  beta_w0_true <- beta_w0
  beta_w_true <- beta_w
  
  delta_true <- delta
  delta_star_true <- delta_star
  
  phi_0_true <- phi_0
  phi_1_true <- phi_1
  
  sigma_y_true <- rep(sigma_y, n_P)
  sigma_true <- sigma
}
