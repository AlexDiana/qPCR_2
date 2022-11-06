# simulate data 

simulateCfromDelta <- function(delta_star, phi_0, phi_1, w_star, 
                               alpha1, alpha2, lambda, sigma_lambda){

  C_star <- ifelse(delta_star == 1, 
                             alpha1[P_star] + alpha2[P_star] * 
                               log(w_star) + 
                               rnorm(nrow(data_standard), sd = sigma_y), 
                             ifelse(delta_star == 1, 
                                    rnorm(length(P_star), lambda, sigma_lambda), 
                                    0))
  
  C_star
  
}

# STARTING VALUES ----


phi_0 <- -26
phi_1 <- 8
n_star <- 10
n_P <- 3
P_star <- sample(1:n_P, n_star, replace = T)
w_star <- rgamma(10, 30, 1)
round(logistic(phi_0 + phi_1 * log(w_star)), 3)

lambda <- 10
sigma_lambda <- 10

# MCMC ------

niter <- 1000

delta_output <- matrix(NA, length(P_star), niter)
gamma_output <- matrix(NA, length(P_star), niter)

# starting values
{
  p_imk_star <- logistic(phi_0 + phi_1 * log(w_star))
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
}

for (iter in 1:niter) {
  
  print(iter)
  
  C_star <- simulateCfromDelta(delta_star, phi_0, phi_1, w_star, 
                               alpha1, alpha2, lambda, sigma_lambda)
  
  
  list_deltagamma <- updateDeltaGammaCpp(Ct, C_star, w_pcr, w_star,
                                         phi_0, phi_1, p0, alpha1,
                                         alpha2, sigma_y, P, P_star, lambda,
                                         sigma_lambda)
  delta <- list_deltagamma$delta
  gamma <- list_deltagamma$gamma
  delta_star <- list_deltagamma$delta_star
  gamma_star <- list_deltagamma$gamma_star
  
  delta_output[,iter] <- delta_star
  gamma_output[,iter] <- gamma_star
    
}


apply(delta_output, 1, mean)
logistic(phi_0 + phi_1 * log(w_star))
