source("~/qPCR_2/functions.R", echo=TRUE)

# simulated data -------

numCycles <- 45

# sites
n <- 100

# sampling occasions
M <- 1 + rep(2, n)#rpois(n, lambda = 3)
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
  beta_b0 <- 2
  beta_w0 <- 0
  theta0 <- .02
  p0 <- .05
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
ncov_b <- 2
X_b <- matrix(rnorm(ncov_b * n), n, ncov_b)
ncov_w <- 2
X_w <- matrix(rnorm(ncov_w * sum(M)), sum(M), ncov_w)
ncov_wt <- 1
X_wt <- matrix(rnorm(ncov_wt * sum(M)), sum(M), ncov_wt)

# covariates coefficients
beta_b <- sample(c(-1,1), size = ncov_b, replace = T)# rnorm(ncov_b)
beta_w <- sample(c(-1,1), size = ncov_w, replace = T)# rnorm(ncov_b)rnorm(ncov_w)
beta_wt <- sample(c(-1,1), size = ncov_wt, replace = T)# rnorm(ncov_b)rnorm(ncov_w)

# amount of biomass at the site
l <- beta_b0 + X_b %*% beta_b + rnorm(n, sd = tau)

# STAGE 1 -----

v1_tilde <- rep(NA, sum(M)) # collected log-biomass
v2_tilde <- rep(NA, sum(M)) # log-biomass from contamination

# true positives 
l_samples <- l[idx_site]
v1 <- l_samples + beta_w0 + X_w %*% beta_w + rnorm(sum(M), sd = sigma)
w1 <- convertPositive(v1)

# assign to 0 sites with log-collected biomass less than 1
v2 <- nu + X_wt %*% beta_wt + rnorm(sum(M), 0, sd = sigma_nu)
w2 <- convertPositive(v2)

w <- w1 + w2
w_pcr <- w[idx_sample_K]

# v <- ifelse(!is.na(v1_tilde), v1, v2)

# v_pcr <- v[idx_sample_K]
# v_pcr <- ifelse(!is.na(v_pcr), v_pcr, -Inf)

# STAGE 2 --------

alpha1 <- rnorm(n_P, mean = alpha1_0, sd = 1)
alpha2 <- rnorm(n_P, mean = alpha2_0, sd = .1)

# standards
p_imk_star <- logistic(phi_0 + phi_1 * log(data_standard$Qty))
delta_star <- sapply(p_imk_star, function(p){
  rbinom(1, 1, p)
})

gamma_star <- sapply(delta_star, function(x){
  if(x == 0){
    rbinom(1, 1, p0)
  } else {
    0
  }
})

data_standard$Ct <- ifelse(delta_star == 1, 
                           alpha1[data_standard$P] + alpha2[data_standard$P] * 
                             log(data_standard$Qty) + 
                             rnorm(nrow(data_standard), sd = sigma_y), 
                           ifelse(gamma_star == 1, 
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

gamma <- sapply(delta, function(x){
  if(x == 0){
    rbinom(1, 1, p0)
  } else {
    0
  }
})

Ct <- ifelse(delta == 1, 
             alpha1[data_real$Plate] + 
               alpha2[data_real$Plate] * log(w_pcr) +
               rnorm(nrow(data_real), sd = sigma_y), 
             ifelse(gamma == 1, 
                    rnorm(length(P), lambda, sigma_lambda), 
                    0))

data_real$Ct <- Ct

# save true values -----------

{
  l_true <- l
  v1_true <- v1
  
  tau_true <- tau
  v2_true <- v2
  
  alpha1_true <- alpha1
  alpha2_true <- alpha2
  
  beta_b0_true <- beta_b0
  beta_b_true <- beta_b
  
  beta_w0_true <- beta_w0
  beta_w_true <- beta_w
  
  nu_true <- nu
  beta_wt_true <- beta_wt
  sigma_nu_true <- sigma_nu
  
  delta_true <- delta
  delta_star_true <- delta_star
  gamma_true <- gamma
  gamma_star_true <- gamma_star
  
  phi_0_true <- phi_0
  phi_1_true <- phi_1
  
  p0_true <- p0
  
  sigma_y_true <- rep(sigma_y, n_P)
  sigma_true <- sigma
}
