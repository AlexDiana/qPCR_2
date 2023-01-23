logistic <- function(x){
  1 / (1 + exp(-x))
}

simulateData <- function(dataParams,
                         hyperParams){
  
  list2env(dataParams,environment())
  list2env(hyperParams,environment())
    
  idx_site <- rep(1:n, M)
  idx_site_K <- rep(idx_site, K)
  idx_sample_K <- rep(1:sum(M), K)

  # standards number and concentration
  {
    data_standard <- data.frame(P = rep(1:n, each = K_standards * length(qty)),
                                Qty = rep(rep(qty, each = K_standards), times = n),
                                Ct = rep(NA, n * K_standards * length(qty)))
  }
  
  # covariates
  # ncov_b <- 2
  X_b <- matrix(rnorm(ncov_b * n), n, ncov_b)
  X_w <- matrix(rnorm(ncov_w * sum(M)), sum(M), ncov_w)
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
  
  trueParams <- list(
    l_true = l,
    v1_true = v1,
    tau_true = tau,
    v2_true = v2,
    alpha1_true = alpha1,
    alpha2_true = alpha2,
    beta_b0_true = beta_b0,
    beta_b_true = beta_b,
    beta_w0_true = beta_w0,
    beta_w_true = beta_w,
    nu_true = nu,
    beta_wt_true = beta_wt,
    sigma_nu_true = sigma_nu,
    delta_true = delta,
    delta_star_true = delta_star,
    gamma_true = gamma,
    gamma_star_true = gamma_star,
    phi_0_true = phi_0,
    phi_1_true = phi_1,
    p0_true = p0,
    sigma_y_true = rep(sigma_y, n_P),
    sigma_true = sigma  
  )
  
  dataParams <- list("data_real" = data_real,
                     "data_standard" = data_standard,
                     "X_b" = X_b,
                     "X_w" = X_w,
                     "X_wt" = X_wt
  )
  
  list("trueParams" = trueParams, 
       "dataParams" = dataParams)
}

convertPositive <- function(x){
  (x > 0) * (exp(x) - 1)
}

# erf <- function(x){
#   ifelse(x > 0,
#      pnorm(x, sd = 1 / sqrt(2)) - pnorm(-x, sd = 1 / sqrt(2)),
#      - (pnorm(-x, sd = 1 / sqrt(2)) - pnorm(x, sd = 1 / sqrt(2))))
#   
# }
# 
# int_1 <- function(b, tau, mu){
#   sqrt(pi)*tau*erf(sqrt(2)*b/(2*tau) - mu*sqrt(2)/(2*tau))/(2*sqrt(pi*tau^2))
# }
# 
# int_2 <- function(x, mu, tau, a_b, b_b, c_b,M_i){
#   sqrt(2)*sqrt(pi) * 
#     exp(-mu^2/(2*tau^2) - c_b/(2*sigma^2) - (mu/tau^2 + b_b/sigma^2)^2/(4*(-1/(2*tau^2) - a_b/(2*sigma^2)))) * 
#     erf(sqrt(2/tau^2 + 2*a_b/sigma^2)*x/2 - (mu/tau^2 + b_b/sigma^2)/sqrt(2/tau^2 + 2*a_b/sigma^2)) / 
#     (2*sqrt(pi*tau^2)*(sqrt(2)*sqrt(pi*sigma^2))^M_i*sqrt(2/tau^2 + 2*a_b/sigma^2))
# }
# 
# tnorm <- function(mu, sigma){
#   
#   x <- rnorm(1, mu, sigma)
#   
#   while(x > 0){
#     x <- rnorm(1, mu, sigma)
#   }
#   
#   x
# }
# 
# updateBtilde <- function(wtilde, sigma, tau, X_b, beta_b, beta_b0, 
#                          X_w, beta_w, beta_w0, nu, sigma_nu){
#   
#   n <- nrow(X_b)
#   btilde <- rep(NA, n)
#   Xbeta_b <- X_b %*% beta_b
#   Xbeta_w <- X_w %*% beta_w
#   
#   for (i in 1:n) {
#     
#     wtilde_i <- wtilde[1:M[i] + sum(M[seq_len(i-1)])] - beta_w0 - (Xbeta_w)[1:M[i] + sum(M[seq_len(i-1)])]
#     
#     priorbtilde <- beta_b0 + Xbeta_b[i]
#     p_0 <- prod(dnorm(wtilde_i, nu, sigma_nu)) * (int_1(0, tau, priorbtilde) - int_1(-1000, tau, priorbtilde))
#     
#     a_b <- M[i]
#     b_b <- sum(wtilde_i)
#     c_b <- sum(wtilde_i^2)
#     p_1 <- int_2(1000, priorbtilde, tau, a_b, b_b, c_b, M[i]) - int_2(0, priorbtilde, tau, a_b, b_b, c_b, M[i])
#     
#     p_b0 <- p_1 / (p_1 + p_0)
#     
#     b0 <- rbinom(1, 1, p_b0) 
#     
#     if(b0 == 1){
#       
#       # b greater than 0
#       
#       prior_var <- 1 / tau^2 + M[i] / sigma^2
#       
#       sum_lik <- 0
#       
#       for(m in 1:M[i]){
#         
#         idx_m <- m + sum(M[seq_len(i-1)])
#         sum_lik <- sum_lik + (wtilde[idx_m] - beta_0w - Xbeta_w[idx_m]) 
#         
#       }
#       
#       mu_btilde <- priorbtilde / tau^2 + sum_lik / sigma^2
#       
#       mean_btilde <- mu_btilde / prior_var
#       
#       btilde[i] <- rnorm(1, mean_btilde, sqrt(1 / prior_var))
#       
#         
#     } else {
#       
#       # b less than 0
#       btilde[i] <- tnorm(priorbtilde, tau)
#       
#     }
#     
#   }
#   
#   btilde
#   
# }

updateOmega <- function(phi_0, phi_1, Ct, v_tilde_pcr){
  
  omega <- sapply(1:length(Ct), function(k){
    rpg(1, phi_0 + phi_1 * v_tilde_pcr[k])
  })
  
  omega
}

updateBetaB <- function(l, X_b, tau, 
                        prior_betab0, 
                        sigma_beta0, sigma_beta){
  
  X_b1 <- cbind(1, X_b) 
  
  Lambda_beta <- t(X_b1) %*% X_b1 / tau^2 + 
    diag(1 / c(sigma_beta0^2, rep(sigma_beta^2, ncov_b)), nrow = ncov_b + 1)
  
  mu_beta <- sapply(seq_along(l), function(i){
    X_b1[i,] * l[i] / tau^2
  })  
  
  if(ncov_b > 0){
    mu_beta <- apply(mu_beta, 1, sum) 
  } else {
    mu_beta <- sum(mu_beta)
  }
  
  mu_beta <- mu_beta + 
    c(prior_betab0, rep(0, ncov_b)) / c(sigma_beta0^2, rep(sigma_beta^2, ncov_b))
  
  mean_beta <- solve(Lambda_beta) %*% mu_beta 
    
  Cov_beta <- solve(Lambda_beta)
  
  beta0beta <- mvrnorm(1, mean_beta, Cov_beta)
  
  list("beta_b0" = beta0beta[1],
    "beta_b" = beta0beta[-1])
}

updateBetaW <- function(l, v, X_w, sigma, sigma_beta, idx_site,
                        isBetaW0){
  
  ncov_w <- ncol(X_w)
  
  if(isBetaW0 != T | ncov_w != 0){
    
    l_samples <- l[idx_site]
    
    if(isBetaW0){
      X_w1 <- X_w
    } else {
      X_w1 <- cbind(1, X_w)
    }
    
    X_w1 <- X_w1#[l_samples > 0,]
    
    v_res <- v - l_samples
    # v_res <- v_res[l_samples > 0]
    
    Lambda_beta <- t(X_w1) %*% X_w1 / sigma^2 + diag(1 / sigma_beta^2, nrow = ncov_w + (1 - isBetaW0))
    
    mu_beta <- sapply(seq_along(v_res), function(i){
      X_w1[i,] * v_res[i] / sigma^2
    })
    
    mu_beta <- apply(mu_beta, 1, sum)
    
    mean_beta <- solve(Lambda_beta) %*% mu_beta
    
    Cov_beta <- solve(Lambda_beta)
    
    beta0beta <- mvrnorm(1, mean_beta, Cov_beta)
    
    if(isBetaW0){
      beta_w0 <- 0
      beta_w <- beta0beta
    } else {
      beta_w0 <- beta0beta[1]
      beta_w = beta0beta[-1]
    }
    
    
    
  }
  
  list("beta_w0" = beta_w0,
       "beta_w" = beta_w)
}

updateBetaWt <- function(v_tilde, v, X_wt, 
                         nu0, sigma_nu, sigma_beta, 
                        ncov_wt){
  
  if(any(v < 0)){
    X_wt_v0 <- X_wt[v < 0,,drop = F]
    
    X_wt_v0 <- cbind(1, X_wt_v0)
    
    v_tilde0 <- v_tilde[v < 0]
    
    Lambda_beta <- t(X_wt_v0) %*% X_wt_v0 / sigma_nu^2 + 
      diag(1 / sigma_beta^2, nrow = ncov_wt + 1)
    
    mu_beta <- sapply(seq_along(v_tilde0), function(i){
      X_wt_v0[i,] * v_tilde0[i] / sigma_nu^2
    })
    
    mu_beta <- apply(mu_beta, 1, sum) + 
      c(nu0, rep(0, ncov_wt)) / rep(1 / sigma_beta^2, ncov_wt + 1)
    
    mean_beta <- solve(Lambda_beta) %*% mu_beta
    
    Cov_beta <- solve(Lambda_beta)
    
    beta0beta <- mvrnorm(1, mean_beta, Cov_beta)
    
    nu <- beta0beta[1]
    beta_wt <-  beta0beta[-1]  
  
  } else {
    
    nu <- rnorm(1, nu0, sigma_nu)
    beta_wt <- rnorm(ncov_wt, 0, sigma_beta)
    
  }
  
  
  list("nu" = nu,
       "beta_wt" = beta_wt)
}

updateSigmaNu <- function(v_tilde, v, X_wt, nu, 
                          a_sigma_nu, b_sigma_nu){
  
  Xwtbetat <- nu + X_wt %*% beta_wt
  
  v_tilde0 <- v_tilde[v < 0]
  
  X_wt_v0 <- Xwtbetat[v < 0]
  
  wtilde_res <- (v_tilde0 - X_wt_v0)^2
  
  sumsq <- sum(wtilde_res)
  n_samples <- length(wtilde_res)
  
  sigma_nu <- sqrt(rinvgamma(a_sigma_nu + n_samples / 2, b_sigma_nu + sumsq / 2))
  
  sigma_nu
}

rinvgamma <- function(a, b){
  1 / rgamma(1, a, b)
}

updateTau <- function(btilde, X_b, beta_b0, beta_b, a_tau, b_tau){
  
  X_b1 <- cbind(1, X_b) 
  beta0beta <- c(beta_b0, beta_b)
  Xbeta <- X_b1 %*% beta0beta
  
  btilde_res <- (btilde - Xbeta)^2
  
  sumsq <- sum(btilde_res)
  n_samples <- length(btilde_res)
  
  tau <- sqrt(rinvgamma(a_tau + n_samples / 2, b_tau + sumsq / 2))
  
  tau
}

updateCtilde <- function(delta, w_tilde, w_tilde_star, alpha){
  
  w_tilde_PCR <- w_tilde[idx_sample_K]
  delta_PCR <- delta[idx_sample_K]
  
  for (i in 1:length(Ctilde)) {
    if(Ctilde[i] == numCycles){
      if(delta_PCR[i] == 1){
        
        print(i)
        
      }  
    }
    
  }
  
}

logf_v <- function(v_i, Ct_i, v_mean, sigma, 
                   alpha1_i, alpha2_i, w_tilde_i, 
                   sigma_y, p0, phi_0, phi_1,
                   lambda, sigma_lambda){
  
  logprior <- dnorm(v_i, v_mean, sd = sigma, log = T)
  
  # pgivendelta0 <- sapply(1:length(Ct_i), function(i){
  #   if(Ct_i[i] == 0){
  #     
  #     1 - p0
  #     
  #   } else {
  #     
  #     p0 * dnorm(Ct_i[i], lambda, sigma_lambda)
  #     
  #   }
  # })
  
  if(v_i > 0){
    
    w_i <- exp(v_i) - 1
    
    loglik <- sum(sapply(1:length(Ct_i), function(i){
      if(Ct_i[i] > 0){
        
        # log(dbinom(1, 1, 
        #        logistic(phi_0 + phi_1 * log(w_i)), 
        #        log = F) *
        #   dnorm(Ct_i[i], alpha1_i[i] + alpha2_i[i] * log(w_i), 
        #         sd = sigma_y[i]) + 
        #     dbinom(0, 1, 
        #            logistic(phi_0 + phi_1 * log(w_i)), 
        #            log = F) *
        #     p0 * dnorm(Ct_i[i], lambda, sigma_lambda))
        log(dnorm(Ct_i[i], alpha1_i[i] + alpha2_i[i] * log(w_i), 
                  sd = sigma_y[i]) + 
              exp(-(phi_0 + phi_1 * log(w_i))) *
              p0 * dnorm(Ct_i[i], lambda, sigma_lambda)) - 
          log(1 + exp(-(phi_0 + phi_1 * log(w_i))))
          
        
      } else {
        
        dbinom(0, 1, 
                   logistic(phi_0 + phi_1 * log(w_i)), 
                   log = T) + log(1 - p0)
           
      }
      
    }))
    
  } else {
    
    loglik <- sum(sapply(1:length(Ct_i), function(i){
      if(Ct_i[i] > 0){
        
        # log(dbinom(1, 1, 
        #            logistic(phi_0 + phi_1 * log(w_tilde_i)), 
        #            log = F) *
        #       dnorm(Ct_i[i], alpha1_i[i] + alpha2_i[i] * log(w_tilde_i), 
        #             sd = sigma_y[i], log = F) + 
        #       dbinom(0, 1, 
        #              logistic(phi_0 + phi_1 * log(w_tilde_i)), 
        #              log = F) *
        #       p0 * dnorm(Ct_i[i], lambda, sigma_lambda))
        
        log(dnorm(Ct_i[i], alpha1_i[i] + alpha2_i[i] * log(w_tilde_i), 
                  sd = sigma_y[i]) + 
              exp(-(phi_0 + phi_1 * log(w_tilde_i))) *
              p0 * dnorm(Ct_i[i], lambda, sigma_lambda)) - 
          log(1 + exp(-(phi_0 + phi_1 * log(w_tilde_i))))
        
        
      } else {
        
        dbinom(0, 1, 
               logistic(phi_0 + phi_1 * log(w_tilde_i)), 
               log = T) + log(1 - p0)
        
      }
      
    }))
    
  }
  
  # loglik <- 0
  # logprior <- 0
  logprior + loglik
  
  
}

logf_v_ratio <- function(v_i, v_i_star,
                         v_tilde_i, v_tilde_i_star,
                         Ct_i, v_mean, sigma, 
                         alpha1_i, alpha2_i,  
                         sigma_y, p0, phi_0, phi_1,
                         lambda, sigma_lambda){
  
  # v_i <- v_current
  # v_tilde_i <- v_tilde_current
  # v_i_star <- v_star
  # v_tilde_i_star <- v_tilde_star
  
  logprior_ratio <- 
    dnorm(v_i_star, v_mean, sd = sigma, log = T) -
  dnorm(v_i, v_mean, sd = sigma, log = T)
  
  w_i <- ifelse(v_i > 0,  
                exp(v_i) - 1, 
                exp(v_tilde_i) - 1)
  w_star_i <- ifelse(v_i_star > 0,  
                exp(v_i_star) - 1, 
                exp(v_tilde_i_star) - 1)
  
  loglik_ratio <- loglik_vv_tilde_ratio(v_i, v_tilde_i, 
                                        v_i_star, v_tilde_i_star,
                                        Ct_i, 
                                        alpha1_i, alpha2_i, 
                                        sigma_y, p0, phi_0, phi_1,
                                        lambda, sigma_lambda)
  
  # loglik_ratio <- sum(sapply(1:length(Ct_i), function(i){
  #   if(Ct_i[i] > 0){
  #     
  #     if(w_i > 0 & w_star_i > 0){
  #       
  #       return(
  #         log(
  #           (dnorm(Ct_i[i], alpha1_i[i] + alpha2_i[i] * log(w_star_i), 
  #                  sd = sigma_y[i]) + 
  #              exp(-(phi_0 + phi_1 * log(w_star_i))) *
  #              p0 * dnorm(Ct_i[i], lambda, sigma_lambda)) /
  #             (dnorm(Ct_i[i], alpha1_i[i] + alpha2_i[i] * log(w_i), 
  #                    sd = sigma_y[i]) + 
  #                exp(-(phi_0 + phi_1 * log(w_i))) *
  #                p0 * dnorm(Ct_i[i], lambda, sigma_lambda))
  #         ) - 
  #           (log(1 + exp(-(phi_0 + phi_1 * log(w_star_i)))) -
  #              log(1 + exp(-(phi_0 + phi_1 * log(w_i)))) )
  #         
  #         )
  #       
  #     } else if(w_i > 0 & w_star_i < 0){
  #       
  #       return(-Inf)
  #       
  #     } else if(w_i < 0 & w_star_i > 0){
  #       
  #       return(Inf)
  #       
  #     } else {
  #       
  #       print("Unlikely scenario")
  #       
  #     }
  #     
  #   } else {
  #     
  #     # dbinom(0, 1, 
  #     #            logistic(phi_0 + phi_1 * log(w_i)), 
  #     #            log = T) + log(1 - p0)
  #     if(w_i > 0 & w_star_i > 0){
  #       
  #       linPred_current <- phi_0 + phi_1 * log(w_i)
  #       linPred_star <- phi_0 + phi_1 * log(w_star_i)
  #       - (linPred_star) - log(1 + exp(-linPred_star)) - 
  #         (- linPred_current - log(1 + exp(-linPred_current)))
  #       
  #       
  #     } else if(w_i > 0 & w_star_i < 0){
  #       
  #       - dbinom(0, 1, 
  #              logistic(phi_0 + phi_1 * log(w_i)),
  #              log = T)
  #       
  #     } else if(w_i < 0 & w_star_i > 0){
  #       
  #       dbinom(0, 1, 
  #              logistic(phi_0 + phi_1 * log(w_star_i)),
  #              log = T)
  #       
  #     } else {
  #       
  #       0
  #       
  #     }
  #     
  #     
  #   }
  #   
  # }))
    
  logprior_ratio + loglik_ratio
  
}

logf_v_int <- function(v_i, Ct_i, v_mean, sigma, 
                   alpha1_i, alpha2_i, w_tilde_i, 
                   sigma_y, p0, phi_0, phi_1,
                   lambda, sigma_lambda, Q_v,
                   nu, sigma_nu){
  
  logprior <- dnorm(v_i, v_mean, sd = sigma, log = T)
  
  # pgivendelta0 <- sapply(1:length(Ct_i), function(i){
  #   if(Ct_i[i] == 0){
  #     
  #     1 - p0
  #     
  #   } else {
  #     
  #     p0 * dnorm(Ct_i[i], lambda, sigma_lambda)
  #     
  #   }
  # })
  
  if(v_i > 0){
    
    w_i <- exp(v_i) - 1
    
    loglik <- sum(sapply(1:length(Ct_i), function(i){
      if(Ct_i[i] > 0){
        
        log(dbinom(1, 1, 
               logistic(phi_0 + phi_1 * log(w_i)), 
               log = F) *
          dnorm(Ct_i[i], alpha1_i[i] + alpha2_i[i] * log(w_i), 
                sd = sigma_y[i]) + 
            dbinom(0, 1, 
                   logistic(phi_0 + phi_1 * log(w_i)), 
                   log = F) *
            p0 * dnorm(Ct_i[i], lambda, sigma_lambda))
          
        
      } else {
        
        dbinom(0, 1, 
                   logistic(phi_0 + phi_1 * log(w_i)), 
                   log = T) + log(1 - p0)
           
      }
      
    }))
    
  } else {
    
    loglik_int <- mean(sapply(1:Q_v, function(i_q){
      
      v_tilde_i <- rnorm(1, nu, sigma_nu)
      w_tilde_i <- (v_tilde_i > 0) * (exp(v_tilde_i) - 1)
      
      sum(sapply(1:length(Ct_i), function(i){
        if(Ct_i[i] > 0){
          
          if(w_tilde_i > 0){
            return(log(dbinom(1, 1, 
                       logistic(phi_0 + phi_1 * log(w_tilde_i)), 
                       log = F) *
                  dnorm(Ct_i[i], alpha1_i[i] + alpha2_i[i] * log(w_tilde_i), 
                        sd = sigma_y[i], log = F) + 
                  dbinom(0, 1, 
                         logistic(phi_0 + phi_1 * log(w_tilde_i)), 
                         log = F) *
                  p0 * dnorm(Ct_i[i], lambda, sigma_lambda)))
          } else {
            return(-Inf)
          }
         
          
          
        } else {
          
          if(w_tilde_i > 0){
            return(dbinom(0, 1, 
                   logistic(phi_0 + phi_1 * log(w_tilde_i)), 
                   log = T) + log(1 - p0))
          } else {
            return(log(1 - p0))
          }
          
        }
      
    }))
  
    }))
    
  }
  
  # loglik <- 0
  # logprior <- 0
  logprior + loglik
  
  
}

logf_vtilde <- function(x, Ct_i, vtilde_mean, sigma_vtilde, 
                   alpha1_i, alpha2_i, 
                   sigma_y, p0, phi_0, phi_1,
                   lambda, sigma_lambda){
  
  logprior <- dnorm(x, vtilde_mean, sd = sigma_vtilde, log = T)

  if(x > 0){
    
    w_i <- exp(x) - 1
    
    loglik <- sum(sapply(1:length(Ct_i), function(i){
      if(Ct_i[i] > 0){
        
        log(dbinom(1, 1, 
               logistic(phi_0 + phi_1 * log(w_i)), 
               log = F) *
          dnorm(Ct_i[i], alpha1_i[i] + alpha2_i[i] * log(w_i), 
                sd = sigma_y[i]) + 
            dbinom(0, 1, 
                   logistic(phi_0 + phi_1 * log(w_i)), 
                   log = F) *
            p0 * dnorm(Ct_i[i], lambda, sigma_lambda))
          
        
      } else {
        
        dbinom(0, 1, 
                   logistic(phi_0 + phi_1 * log(w_i)), 
                   log = T) + log(1 - p0)
           
      }
      
    }))
    
  } else {
    
    loglik <- sum(sapply(1:length(Ct_i), function(i){
      if(Ct_i[i] > 0){
        
        log(p0 * dnorm(Ct_i[i], lambda, sigma_lambda))
        
        
      } else {
        
        log(1 - p0)
        
      }
      
    }))
    
  }
  
  logprior + loglik
  
}

findMustar_vtilde <- function(Ct_i, 
                              vtilde_mean, sigma_vtilde,
                              alpha1_i, alpha2_i, w_tilde_i, 
                              sigma_y, p0, phi_0, phi_1,
                              lambda, sigma_lambda){

  fn_optim <- function(x){
    
    - logf_vtilde_cpp(x, Ct_i, v_tilde_mean, sigma_nu, 
                alpha1_i, alpha2_i, 
                sigma_y, p0, phi_0, phi_1,
                lambda, sigma_lambda)
    
  }
  
  vtilde_star1 <- optimize(fn_optim, lower = -20, upper = 0)
  vtilde_star2 <- optimize(fn_optim, lower = 0, upper = 20)
  
  v_tilde_star <- ifelse(vtilde_star1$objective < vtilde_star2$objective,
                         vtilde_star1$minimum, vtilde_star2$minimum)
  
  v_grid <- seq(-5, 5, length.out = 100)
  #
  # v_grid <- seq(v_true[idx_current] - 10, v_true[idx_current] + 3, length.out = 100)
  y_grid <- sapply(v_grid, function(x){
    logf_vtilde(x, Ct_i, v_tilde_mean, sigma_nu,
                alpha1_i, alpha2_i,
                sigma_y, p0, phi_0, phi_1,
                lambda, sigma_lambda)
  })
  # 
  qplot(v_grid, y_grid)
  v_tilde_star
  
}

updateV <- function(v, v_tilde, Ct, l,
                    beta_w0, X_w, beta_w,
                    omega,
                    alpha1, alpha2, P, 
                    nu, X_wt, beta_wt, sigma_nu, sigmas_v,
                    phi_0, phi_1, lambda, sigma_lambda){
  
  Xwbetaw <- beta_w0 + X_w %*% beta_w
  Xwtbetawt <- nu + X_wt %*% beta_wt
  
  for (i in 1:n) {
    
    l_i <- l[i]
    
    for (m in 1:M[i]) {
      
      idx_current <- m + sum(M[seq_len(i-1)])
      v_mean <- l_i + Xwbetaw[idx_current]
      
      v_current <- v[idx_current]
      
      Ct_i <- Ct[idx_sample_K == idx_current]
      
      alpha1_i <- alpha1[P[idx_sample_K == idx_current]]
      alpha2_i <- alpha2[P[idx_sample_K == idx_current]]
      
      sigma_y_i <- sigma_y[P[idx_sample_K == idx_current]]
      omega_i <- omega[idx_sample_K == idx_current]
      delta_i <- delta[idx_sample_K == idx_current]
      k_i <- delta_i - .5
      
      # sample v1_tilde > or < 0
      {
        m_k_lik <- 
          sum(((Ct_i - alpha1_i) / alpha2_i) / sigma_y_i^2) +
          sum(phi_1 * (k_i - phi_0 * omega))
        prec_lik <- sum(1 / sigma_y_i^2) + sum(phi_1^2 * omega)
        
        m_k_lik <- 
          sum(phi_1 * (k_i - phi_0 * omega))
        prec_lik <- sum(phi_1^2 * omega)
        
        var_lik <- 1 / prec_lik
        
        mean_k_lik <-  m_k_lik * var_lik
        
        rnorm(1, mean_k_lik, sqrt(var_lik))
      }

      v_tilde_mean <- Xwtbetawt[idx_current]
      
      v_star <- rnorm(1, v_current, sd = sigmas_v[idx_current])
      
      if(v_current > 0 & v_star > 0){
        
        logproposal_v <- 0
        
        logprior_v <- 0
        
        v_tilde_current <- NA
        
        v_tilde_star <- NA

      } else if(v_current > 0 & v_star < 0){

        v_tilde_star <- findMustar_vtilde_cpp(Ct_i, 
                                              v_tilde_mean, sigma_nu,
                                              alpha1_i, alpha2_i, 
                                              sigma_y_i, p0, phi_0, phi_1,
                                              lambda, sigma_lambda)
        
        # v_tilde_star <- findMustar_vtilde(Ct_i, 
        #                                   vtilde_mean, sigma_nu,
        #                                   alpha1_i, alpha2_i, w_tilde_i, 
        #                                   sigma_y, p0, phi_0, phi_1,
        #                                   lambda, sigma_lambda)
        
        v_tilde_current <- rnorm(1, v_tilde_star, 1)
        
        logproposal_v <- - dnorm(v_tilde_current, v_tilde_star, 1, log = T)
        
        logprior_v <- dnorm(v_tilde_current, v_tilde_mean, sd = sigma_nu,
                            log = T)

      } else if(v_current < 0 & v_star > 0){

        v_tilde_current <- v_tilde[idx_current]
        
        v_tilde_star <- findMustar_vtilde_cpp(Ct_i, 
                                          v_tilde_mean, sigma_nu,
                                          alpha1_i, alpha2_i,  
                                          sigma_y_i, p0, phi_0, phi_1,
                                          lambda, sigma_lambda)
        
        logproposal_v <- dnorm(v_tilde_current, v_tilde_star, 1, log = T)
        
        logprior_v <- - dnorm(v_tilde_current, v_tilde_mean, sd = sigma_nu,
                              log = T)

      } else {

        logproposal_v <- 0
        v_tilde_current <- v_tilde[idx_current]
        v_tilde_star <- v_tilde[idx_current]
        
        logprior_v <- 0

      }
      
      # w_tilde_i <- (exp(v_tilde_current) - 1) * (v_tilde_current > 0)
      
      # w_tilde_i <- w_tilde[idx_current]
      
      #
      
      # f_to_integrate <- function(x){
      #   exp(logf_v(x, Ct_i, v_mean, sigma, 
      #          alpha1_i, alpha2_i, w_tilde_i, 
      #          sigma_y, p0, phi_0, phi_1,
      #          lambda, sigma_lambda))
      # }
      # 
      # p1_prop <- integrate(Vectorize(f_to_integrate), upper = 15, lower = .001, subdivisions = 10000)$value
      # p0_prop <- integrate(Vectorize(f_to_integrate), upper = -.001, lower = -15, subdivisions = 1000)$value
      # 
      # p1_prop / (p1_prop + p0_prop)
      
      #
      # v_grid <- seq(v_current - 1, 12, length.out = 100)
      # 
      # # v_grid <- seq(v_true[idx_current] - 3, v_true[idx_current] + 3, length.out = 100)
      # y_grid <- sapply(v_grid, function(x){
      #   logf_v_cpp(x, Ct_i, v_mean, sigma,
      #   # logf_v(x, Ct_i, v_mean, sigma,
      #              alpha1_i, alpha2_i, w_tilde_i,
      #              sigma_y_i, p0, phi_0, phi_1,
      #              lambda, sigma_lambda)
      # })
      # # 
      # qplot(v_grid, y_grid)
      # #+
      # #   geom_vline(aes(xintercept = v_true[idx_current]))
      # qplot(v_grid, exp(y_grid) / sum(exp(y_grid)))# +
      #   geom_vline(aes(xintercept = v_current))
      
      # logposterior_current <- logf_v_cpp(v_current, Ct_i, v_mean, sigma, 
      # # logposterior_current <- logf_v(v_current, Ct_i, v_mean, sigma,
      #                                alpha1_i, alpha2_i, w_tilde_i,
      #                                sigma_y, p0, phi_0, phi_1,
      #                                lambda, sigma_lambda)
      # 
      # logposterior_star <- logf_v_cpp(v_star, Ct_i, v_mean, sigma, 
      # # logposterior_star <- logf_v(v_star, Ct_i, v_mean, sigma, 
      #                                alpha1_i, alpha2_i, w_tilde_i, 
      #                                sigma_y, p0, phi_0, phi_1,
      #                                lambda, sigma_lambda)
      
      logposterior_ratio <- logf_v_ratio(v_current, v_star,
                                         v_tilde_current, v_tilde_star,
                                         Ct_i, v_mean, sigma, 
                                         alpha1_i, alpha2_i,  
                                         sigma_y, p0, phi_0, phi_1,
                                         lambda, `sigma_lambda`) 
      
      logratio <- 
        logposterior_ratio  +
        # logposterior_star - logposterior_current + 
        logprior_v + logproposal_v
      
      # exp(logposterior_star - logposterior_current)
      if(runif(1) < exp(logratio)){
        
        v[idx_current] <- v_star
        
      }
      
      if(v[idx_current] < 0){
        
        v_tilde[idx_current] <- v_tilde_current
        
      } else {
        
        v_tilde[idx_current] <- NA
        
      }
        
    }
    
  }
  
  list("v" = v,
       "v_tilde" = v_tilde)
  
}

logf_v_tilde <- function(v_tilde_i, Ct_i, nu, sigma_nu, 
                   alpha1_i, alpha2_i, 
                   sigma_y, p0, phi_0, phi_1,
                   lambda, sigma_lambda){
  
  logprior <- dnorm(v_tilde_i, nu, sd = sigma_nu, log = T)
  
  if(v_tilde_i > 0){
    
    w_tilde_i <- exp(v_tilde_i) - 1
    
    loglik <- sum(sapply(1:length(Ct_i), function(i){
      if(Ct_i[i] > 0){
        
        log(
          dbinom(1, 1, 
                 logistic(phi_0 + phi_1 * log(w_tilde_i)), 
                 log = F) *
            dnorm(Ct_i[i], alpha1_i[i] + alpha2_i[i] * log(w_tilde_i), 
                  sd = sigma_y) + 
            dbinom(0, 1, 
                   logistic(phi_0 + phi_1 * log(w_tilde_i)), 
                   log = F) *
            p0 * dnorm(Ct_i[i], lambda, sigma_lambda)
        )
        
      } else {
        
        dbinom(0, 1, 
               logistic(phi_0 + phi_1 * log(w_tilde_i)), 
               log = T) + log(1 - p0)
        
      }
      
    }))
    
  } else {
    
    loglik <- sum(sapply(1:length(Ct_i), function(i){
      if(Ct_i[i] > 0){
        
        log(p0) + dnorm(Ct_i[i], lambda, sigma_lambda, log = T)
        
        
      } else {
        
        log(1 - p0)
        
      }
      
    }))
    
  }
  
  logprior + loglik
  
  
}


updateVtilde <- function(v_tilde, v, Ct, l, beta_w0, X_w, beta_w,
                         alpha1, alpha2, P, nu, sigma_nu,
                         phi_0, phi_1, lambda, sigma_lambda){
  
  Xwtbetawt <- nu + X_wt %*% beta_wt
  
  for (i in 1:n) {
    
    l_i <- l[i]
    
    for (m in 1:M[i]) {
      
      idx_current <- m + sum(M[seq_len(i-1)])
      
      if(v[idx_current] < 0){
        
        v_tilde_current <- v_tilde[idx_current]
        
        Ct_i <- Ct[idx_sample_K == idx_current]
        
        alpha1_i <- alpha1[P[idx_sample_K == idx_current]]
        alpha2_i <- alpha2[P[idx_sample_K == idx_current]]
        sigmay_i <- sigma_y[P[idx_sample_K == idx_current]]
        
        v_grid <- seq(v_tilde_current - 6, v_tilde_current + 5, length.out = 100)
        y_grid <- sapply(v_grid, function(x){
          logf_v_tilde_cpp(x, Ct_i, nu, sigma_nu,
                 alpha1_i, alpha2_i, 
                 sigmay_i, p0, phi_0, phi_1,
                 lambda, sigma_lambda)
        })

        qplot(v_grid, y_grid) #+
        # geom_vline(aes(xintercept = v_current))
        # qplot(v_grid, exp(y_grid) / sum(exp(y_grid))) +
        # geom_vline(aes(xintercept = v_current))
        
        v_tilde_star <- rnorm(1, v_tilde_current, sd = .5)
        
        logposterior_current <- logf_v_tilde(v_tilde_current, Ct_i, nu, sigma_nu, 
                                       alpha1_i, alpha2_i, 
                                       sigma_y, p0, phi_0, phi_1,
                                       lambda, sigma_lambda)
        
        logposterior_star <- logf_v_tilde(v_tilde_star, Ct_i, nu, sigma_nu, 
                                    alpha1_i, alpha2_i, 
                                    sigma_y, p0, phi_0, phi_1,
                                    lambda, sigma_lambda)
        
        if(runif(1) < exp(logposterior_star - logposterior_current)){
          
          v_tilde[idx_current] <- v_tilde_star
          
          
        }  
        
      }
        
    }
    
  }
  
  v_tilde
  
}

loglik_vv_tilde_ratio <- function(v_i, v_tilde_i, 
                                v_i_star, v_tilde_i_star,
                                Ct_i, 
                                alpha1_i, alpha2_i, 
                                sigma_y, p0, phi_0, phi_1,
                                lambda, sigma_lambda){
  
  w_i <- ifelse(v_i > 0,  
                exp(v_i) - 1, 
                exp(v_tilde_i) - 1)
  w_star_i <- ifelse(v_i_star > 0,  
                     exp(v_i_star) - 1, 
                     exp(v_tilde_i_star) - 1)
  
  loglik_ratio <- sum(sapply(1:length(Ct_i), function(i){
    if(Ct_i[i] > 0){
      
      if(w_i > 0 & w_star_i > 0){
        
        # log(
        #   (dnorm(Ct_i[i], alpha1_i[i] + alpha2_i[i] * log(w_star_i), 
        #          sd = sigma_y[i]) + 
        #      exp(-(phi_0 + phi_1 * log(w_star_i))) *
        #      p0 * dnorm(Ct_i[i], lambda, sigma_lambda))/ 
        #     (dnorm(Ct_i[i], alpha1_i[i] + alpha2_i[i] * log(w_i), 
        #            sd = sigma_y[i]) + 
        #        exp(-(phi_0 + phi_1 * log(w_i))) *
        #        p0 * dnorm(Ct_i[i], lambda, sigma_lambda))
        #   )  - 
        #   log( (1 + exp(-(phi_0 + phi_1 * log(w_star_i)))) / 
        #           (1 + exp(-(phi_0 + phi_1 * log(w_i)))) ) 
        
        log(
          (dnorm(Ct_i[i], alpha1_i[i] + alpha2_i[i] * log(w_star_i),
                 sd = sigma_y[i]) +
             exp(-(phi_0 + phi_1 * log(w_star_i))) *
             p0 * dnorm(Ct_i[i], lambda, sigma_lambda)) /
            (dnorm(Ct_i[i], alpha1_i[i] + alpha2_i[i] * log(w_i),
                   sd = sigma_y[i]) +
               exp(-(phi_0 + phi_1 * log(w_i))) *
               p0 * dnorm(Ct_i[i], lambda, sigma_lambda))
        ) -
          (log(1 + exp(-(phi_0 + phi_1 * log(w_star_i)))) -
             log(1 + exp(-(phi_0 + phi_1 * log(w_i)))) )
        
      } else if(w_star_i > 0 & w_i < 0){
        
        log(
          (dbinom(1, 1, prob = 
                    logistic(phi_0 + phi_1 * log(w_star_i)), log = F) * 
             dnorm(Ct_i[i], alpha1_i[i] + alpha2_i[i] * log(w_star_i), 
                 sd = sigma_y[i], log = F) + 
             dbinom(0, 1, prob = 
                      logistic(phi_0 + phi_1 * log(w_star_i)), log = F) *
             p0 * dnorm(Ct_i[i], lambda, sigma_lambda))) -
            (log(p0) + dnorm(Ct_i[i], lambda, sigma_lambda, log = T))
        
      } else if(w_i > 0 & w_star_i < 0){
        
        (log(p0) + dnorm(Ct_i[i], lambda, sigma_lambda, log = T)) - 
        log(
          (dbinom(1, 1, prob = 
                    logistic(phi_0 + phi_1 * log(w_i)), log = F) * 
             dnorm(Ct_i[i], alpha1_i[i] + alpha2_i[i] * log(w_i), 
                   sd = sigma_y[i], log = F) + 
             dbinom(0, 1, prob = 
                      logistic(phi_0 + phi_1 * log(w_i)), log = F) *
             p0 * dnorm(Ct_i[i], lambda, sigma_lambda)))
          
        
      } else {
        
        0
        
      }
      
    } else {
      
      if(w_i > 0 & w_star_i > 0){
        
        linPred_current <- phi_0 + phi_1 * log(w_i)
        linPred_star <- phi_0 + phi_1 * log(w_star_i)
        - (linPred_star) - log(1 + exp(-linPred_star)) - 
          (- linPred_current - log(1 + exp(-linPred_current)))
        
      } else if(w_star_i > 0 & w_i < 0){ 
        
        dbinom(1, 1, prob = 
                 logistic(phi_0 + phi_1 * log(w_star_i)), log = T)
        
        
      } else if(w_i > 0 & w_star_i < 0){
        
        - dbinom(1, 1, prob = 
                 logistic(phi_0 + phi_1 * log(w_i)), log = T)
        
        
      } else {
        
        0
        
      }
      
    }
    
  }))
  
  loglik_ratio
  
}

logf_v_vtilde_ratio <- function(v_i, v_tilde_i, 
                                v_i_star, v_tilde_i_star,
                                Ct_i, mean_v, sigma,
                                mean_vtilde, sigma_nu, 
                                alpha1_i, alpha2_i, 
                                sigma_y, p0, phi_0, phi_1,
                                lambda, sigma_lambda){
  
  # v_i <- v_current
  # v_tilde_i <- v_tilde_current
  # v_i_star <- v_star
  # v_tilde_i_star <- v_tilde_star
  
  
  logprior_ratio <- 
    dnorm(v_i_star, mean_v, sd = sigma, log = T) -
    dnorm(v_i, mean_v, sd = sigma, log = T)
  
  # logprior <- dnorm(v_i, mean_v, sd = sigma, log = T)
  
  if(v_i_star < 0){
  
    logprior_ratio <- logprior_ratio + 
      dnorm(v_tilde_i_star, mean_vtilde, sd = sigma_nu, log = T)  
    
  }
  
  if(v_i < 0){
  
    logprior_ratio <- logprior_ratio - 
      dnorm(v_tilde_i, mean_vtilde, sd = sigma_nu, log = T)  
    
  }
  
  loglik_ratio <- loglik_vv_tilde_ratio(v_i, v_tilde_i, 
                                        v_i_star, v_tilde_i_star,
                                        Ct_i, 
                                        alpha1_i, alpha2_i, 
                                        sigma_y, p0, phi_0, phi_1,
                                        lambda, sigma_lambda)
  
  # w_i <- ifelse(v_i > 0,  
  #               exp(v_i) - 1, 
  #               exp(v_tilde_i) - 1)
  # w_star_i <- ifelse(v_i_star > 0,  
  #                    exp(v_i_star) - 1, 
  #                    exp(v_tilde_i_star) - 1)
  # 
  # loglik_ratio <- sum(sapply(1:length(Ct_i), function(i){
  #   if(Ct_i[i] > 0){
  #     
  #     if(w_i > 0 & w_star_i > 0){
  #      
  #       log(
  #         (dnorm(Ct_i[i], alpha1_i[i] + alpha2_i[i] * log(w_star_i), 
  #                sd = sigma_y[i]) + 
  #            exp(-(phi_0 + phi_1 * log(w_star_i))) *
  #            p0 * dnorm(Ct_i[i], lambda, sigma_lambda)) /
  #           (dnorm(Ct_i[i], alpha1_i[i] + alpha2_i[i] * log(w_i), 
  #                  sd = sigma_y[i]) + 
  #              exp(-(phi_0 + phi_1 * log(w_i))) *
  #              p0 * dnorm(Ct_i[i], lambda, sigma_lambda))
  #       ) - 
  #         (log(1 + exp(-(phi_0 + phi_1 * log(w_star_i)))) -
  #            log(1 + exp(-(phi_0 + phi_1 * log(w_i)))) ) 
  #       
  #     } else if(w_i > 0 & w_star_i < 0){
  #       
  #       
  #       
  #       
  #     } else if(w_i > 0 & w_star_i < 0){
  #       
  #       
  #       
  #     } else {
  #       
  #       
  #       
  #     }
  #     
  #     
  #     
  #   } else {
  #     
  #     # dbinom(0, 1, 
  #     #            logistic(phi_0 + phi_1 * log(w_i)), 
  #     #            log = T) + log(1 - p0)
  #     linPred_current <- phi_0 + phi_1 * log(w_i)
  #     linPred_star <- phi_0 + phi_1 * log(w_star_i)
  #     - (linPred_star) - log(1 + exp(-linPred_star)) - 
  #       (- linPred_current - log(1 + exp(-linPred_current)))
  #     
  #     
  #   }
  #   
  # }))
  
  logprior_ratio + loglik_ratio
  
  
}

dtnorm <- function(x, mu, sigma, pos){
  
  if(pos){
    dnorm(x, mu, sigma) / (1 - pnorm(0, mu, sigma))
  } else {
    dnorm(x, mu, sigma) / pnorm(0, mu, sigma)  
  }
  
}


rtnorm_old_neg <- function(mu, sigma){
  
  x <- rnorm(1, mu, sigma)
  while(x > 0){
    x <- rnorm(1, mu, sigma)
  }
  
  x
  
}

rtnorm <- function(mu, sigma, pos){
  if(pos){
    tnorm(mu, sigma, 0)
  } else {
    - tnorm(-mu, sigma, 0)  
  }
  
}

updateV_Vtilde <- function(v, v_tilde, Ct, l, beta_w0, X_w, beta_w,
                           alpha1, alpha2, P, nu, sigma_nu, sigmas_vvtilde,
                           phi_0, phi_1, lambda, sigma_lambda){
  
  Xwbetaw <- beta_w0 + X_w %*% beta_w
  Xwtbetawt <- nu + X_wt %*% beta_wt
  
  for (i in 1:n) {
    
    l_i <- l[i]
    
    for (m in 1:M[i]) {
      
      idx_current <- m + sum(M[seq_len(i-1)])
      
      v_current <- v[idx_current]
      v_tilde_current <- v_tilde[idx_current]
      
      if(v_current > 0 | v_tilde_current > 0){
        
        Ct_i <- Ct[idx_sample_K == idx_current]
        alpha1_i <- alpha1[P[idx_sample_K == idx_current]]
        alpha2_i <- alpha2[P[idx_sample_K == idx_current]]
        sigma_y_i <- sigma_y[P[idx_sample_K == idx_current]]
        
        mean_v <- l_i + Xwbetaw[idx_current]
        mean_vtilde <- Xwtbetawt[idx_current]
        
        if(v_current > 0){ #  1 to 2
          
          v_star <- rtnorm(mean_v, sigma, 0)
          v_tilde_star <- rnorm(1, v_current, sd = sigmas_vvtilde[i])
          
          logproposal <- 
            log(dtnorm(v_current, v_tilde_star, sigmas_vvtilde[i], 1)) - 
            (log(dtnorm(v_star, mean_v, sigma, 0)) + 
               dnorm(v_tilde_star, v_current, sigmas_vvtilde[i], log = T))
          
        } else { # 2 to 1
          
          v_star <- rtnorm(v_tilde_current, sigmas_vvtilde[i], 1)
          v_tilde_star <- NA
          
          logproposal <- 
            (log(dtnorm(v_current, mean_v, sigma, 0)) + 
               dnorm(v_tilde_current, v_star, sigmas_vvtilde[i], log = T)) - 
            log(dtnorm(v_star, v_tilde_current, sigmas_vvtilde[i], 1))
          
        } 
        
        logposterior_ratio <- logf_v_vtilde_ratio(v_current, v_tilde_current, 
                                                  v_star, v_tilde_star,
                                                  Ct_i, mean_v, sigma,
                                                  mean_vtilde, sigma_nu, 
                                                  alpha1_i, alpha2_i, 
                                                  sigma_y, p0, phi_0, phi_1,
                                                  lambda, sigma_lambda) 
        
        log_mh_ratio <- 
          logposterior_ratio + 
          logproposal
        
        if(runif(1) < exp(log_mh_ratio)){
          
          v[idx_current] <- v_star
          v_tilde[idx_current] <- v_tilde_star
          
        }    
        
      }
      
    }
    
  }
  
  list("v" = v,
       "v_tilde" = v_tilde)
  
}

erf_l <- function(x, sigma, w){
  sqrt(pi)*sigma*erf(sqrt(2)*x/(2*sigma) - w*sqrt(2)/(2*sigma))/(2*sqrt(pi*sigma^2))
}

logf_l <- function(l_i, v_i, 
                   mean_l, tau, Xwbeta_i, sigma){
  
  logprior <- dnorm(l_i, mean_l, tau, log = T)
  
  loglik <- sum(sapply(1:length(v_i), function(i){
    dnorm(v_i[i], l_i + Xwbeta_i[i], sd = sigma, log = T)
  }))
  
  logprior + loglik
}

logfp_l <- function(l_i, v_i, 
                   w_i, 
                   mean_l, tau, Xwbeta_i, sigma){
  
  logprior_p <- - (1 / (2 * tau^2)) * (l_i - mean_l) 
  # logprior <- dnorm(l_i, mean_l, tau, log = T)
  
  b_i <- (l_i > 0) * (exp(l_i) - 1)
  
  if(l_i > 0){
    
    loglik <- sum(sapply(1:length(w_i), function(i){
      if(v_i[i] == -Inf){
        log(erf_l(0, sigma, log(b_i) + Xwbeta_i[i]) - erf_l(-Inf, sigma, log(b_i) + Xwbeta_i[i]))
      } else {
        dnorm(v_i[i], log(b_i) + Xwbeta_i[i], sd = sigma, log = T)
      }
      
    }))
    
  } else {
    
    loglik <- sum(sapply(1:length(w_i), function(i){
      if(w_i[i] == 0){
        return(0)
      } else {
        return(-Inf)
      }
    }))
    
  }
  
  logprior + loglik
}

updateL <- function(v, X_b, beta_b0, beta_b,
                    beta_w0, X_w, beta_w, tau){
  
  Xbeta <- beta_b0 + X_b %*% beta_b
  Xwbeta <- beta_w0 + X_w %*% beta_w
  
  for (i in 1:n) {
    
    l_current <- l[i]
    
    priormean_l <- Xbeta[i]
    
    v_i <- v[1:M[i] + sum(M[seq_len(i-1)])]
    Xwbeta_i <- Xwbeta[1:M[i] + sum(M[seq_len(i-1)])]
    
    prior_var <- 1 / tau^2 + M[i] / sigma^2
    
    sum_lik <- sum(v_i - Xwbeta_i)
    
    # sum_lik <- 0
    
    # for(m in 1:M[i]){
    #   
    #   idx_m <- m + sum(M[seq_len(i-1)])
    #   sum_lik <- sum_lik + (wtilde[idx_m] - beta_0w - Xbeta_w[idx_m]) 
    #   
    # }
    
    mu_l <- priormean_l / tau^2 + sum_lik / sigma^2
    
    mean_l <- mu_l / prior_var
    
    l[i] <- rnorm(1, mean_l, sqrt(1 / prior_var))
    
    # l_grid <- seq(l_current - 3, l_current + 3, length.out = 1000)
    # y_grid <- sapply(l_grid, function(x){
    #   logf_l(x, v_i, mean_l, tau, Xwbeta_i, sigma)
    # })
    # 
    # qplot(l_grid, exp(y_grid) / sum(exp(y_grid)))
    
    
      
  }
  
  l
}

updateDeltaGamma <- function(Ct, C_star, w_pcr, w_star,
                             phi_0, phi_1, p0){
  
  delta <- rep(NA, length(Ct))
  gamma <- rep(NA, length(Ct))
  
  for (i in seq_along(w_pcr)) {
    
    if(Ct[i] > 0){
      
      w_current <- w_pcr[i]
      
      if(w_current > 0){
        
        Xb_current <- phi_0 + phi_1 * log(w_current)
        
        logprior_delta1 <- - log(1 + 1 / exp(Xb_current))
        
        logprior_delta0 <- - Xb_current - log(1 + 1 / exp(Xb_current))
        
        logprior_gamma_1 <- log(p0)
        
        loglik_delta1 <- dnorm(Ct[i], alpha1[P[i]] + alpha2[P[i]] * log(w_current), 
                               sd = sigma_y[P[i]], log = T)
        
        loglik_gamma1 <- dnorm(Ct[i], lambda, sd = sigma_lambda, log = T)
        
        log_delta1 <- logprior_delta1 + loglik_delta1
        
        log_gamma1 <- logprior_delta0 +  logprior_gamma_1 + loglik_gamma1
        
        p_delta1 <- exp(log_delta1 - log_gamma1) / (1 + exp(log_delta1 - log_gamma1))
        
        delta[i] <- rbinom(1, 1, p_delta1)
        gamma[i] <- 1 - delta[i]
        
      } else {
        
        delta[i] <- 0
        gamma[i] <- 1
        
      }
      
    } else {
      
      delta[i] <- 0
      gamma[i] <- 0
      
    }
    
  }
  
  delta_star <- rep(NA, length(C_star))
  gamma_star <- rep(NA, length(C_star))
  
  for (i in seq_along(w_star)) {
    
    if(C_star[i] > 0){
      
      w_current <- w_star[i]
      
      Xb_current <- phi_0 + phi_1 * log(w_current)
      
      logprior_delta1 <- - log(1 + 1 / exp(Xb_current))
      
      logprior_delta0 <- - Xb_current - log(1 + 1 / exp(Xb_current))
      
      logprior_gamma_1 <- log(p0)
      
      loglik_delta1 <- dnorm(C_star[i], alpha1[P_star[i]] + alpha2[P_star[i]] * log(w_current), 
                             sd = sigma_y[P_star[i]], log = T)
      
      loglik_gamma1 <- dnorm(C_star[i], lambda, sd = sigma_lambda, log = T)
      
      log_delta1 <- logprior_delta1 + loglik_delta1
      
      log_gamma1 <- logprior_delta0 +  logprior_gamma_1 + loglik_gamma1
      
      # p_delta1 <- exp(log_delta1 - log_gamma1) / (1 + exp(log_delta1 - log_gamma1))
      
      p_delta1 <- exp(log_delta1) / (exp(log_gamma1) + exp(log_delta1))
      
      delta_star[i] <- rbinom(1, 1, p_delta1)
      gamma_star[i] <- 1 - delta_star[i]
      
    } else {
      
      delta_star[i] <- 0
      gamma_star[i] <- 0
      
    }
    
  }
  
  list("delta" = delta,
       "gamma" = gamma,
       "delta_star" = delta_star,
       "gamma_star" = gamma_star)
}

loglik_phi <- function(phi_01,
                       X_phi, 
                       delta_w){
  
  Xphi <- X_phi %*% phi_01
  
  sum(sapply(1:nrow(X_phi), function(i){
    
    if(delta_w[i] == 1){
      
      return(- log(1 + exp(- Xphi[i])))
      
    } else {
      
      return(- Xphi[i] - log(1 + exp(- Xphi[i])))
      
    }
    
  }))
  
}

der2_loglik_phi <- function(phi_01,
                            X_phi){
  
  Xphi <- X_phi %*% phi_01
  
  Hess <- matrix(NA, 2, 2)
  
  Hess[1,1] <- sum(sapply(1:nrow(X_phi), function(i){
    
    # if(delta_w[i] == 1){
      
      return(- exp(-Xphi[i])/((1 + exp(-Xphi[i]))^2))
      
    # } else {
    #   
    #   return(- Xphi[i] - log(1 + exp(- Xphi[i])))
    #   
    # }
    
  }))
  
  Hess[1,2] <- sum(sapply(1:nrow(X_phi), function(i){
    
    # if(delta_w[i] == 1){
      
      return(-X_phi[i,2] * exp(-Xphi[i])/(1 + exp(-Xphi[i]))^2)
      
    # } else {
    #   
    #   return(- Xphi[i] - log(1 + exp(- Xphi[i])))
    #   
    # }
    
  }))
  
  Hess[2,2] <- sum(sapply(1:nrow(X_phi), function(i){
    
    # if(delta_w[i] == 1){
      
      return(-X_phi[i,2]^2*exp(-Xphi[i])/(1 + exp(-Xphi[i]))^2)
      
    # } else {
    #   
    #   return(- Xphi[i] - log(1 + exp(- Xphi[i])))
    #   
    # }
    
  }))
  
  Hess[2,1] <- Hess[1,2]
  
  Hess
}

mvrnorm_chol <- function(mu, chol_Sigma){
  
  Y <- rnorm(ncol(chol_Sigma_star))
  
  mu + t(Y %*% chol_Sigma)
  # arma::mat Y = arma::randn(1, ncols);
  # return mu + arma::trans(Y * arma::chol(Sigma))
  
}

updatePhi <- function(Ct, C_star, w, w_star, 
                      delta, delta_star,
                      phi_0, phi_1, 
                      mean_phi, sigma_phi,
                      Sigma_prop){
  
  
  w_all <- c(w_pcr, w_star)
  delta_all <- c(delta, delta_star)
  
  phi_current <- c(phi_0, phi_1)
  
  idx_w <- which(w_all > 0)
  
  delta_w <- as.numeric(delta_all[idx_w] == 1)
  
  X_phi <- cbind(1, log(w_all[idx_w]))
  
  # phi_01_star <- mvrnorm(1, c(phi_0, phi_1), Sigma = Sigma_prop)
  
  phi_star <- optim(phi_current, fn = function(phi){
    - loglik_phi_cpp(phi, X_phi[,2], delta_w, 
                     mean_phi, diag(sigma_phi, nrow = 2))
  }, method = "BFGS", control = list(reltol = .0001))$par
  
  Hess_phi_star <- hes_loglik_phi_cpp(phi_star, X_phi[,2], delta_w)
  M <- - Hess_phi_star
  Sigma_star <- solve(M)
  
  phi_new <- mvrnorm_inv(phi_star, M)
  # phi_new <- mvrnorm(mu = phi_star, Sigma = Sigma_star)
  
  logproposal_current <- dmvnorm_cpp_inv(phi_current, phi_star, M, 1)
  logproposal_star <- dmvnorm_cpp_inv(phi_new, phi_star, M, 1)
  # logproposal_current <- dmvnorm_cpp(phi_current, phi_star, Sigma_star, 1)
  # logproposal_star <- dmvnorm_cpp(phi_new, phi_star, Sigma_star, 1)
  
  loglik_star <- loglik_phi(phi_new, X_phi, delta_w)
  loglik_current <- loglik_phi(phi_current, X_phi, delta_w)
  
  logprior_star <- dmvnorm_cpp(phi_star, mean_phi, diag(sigma_phi, nrow = 2), 1)
  logprior_current <- dmvnorm_cpp(phi_current, mean_phi, diag(sigma_phi, nrow = 2), 1)
  
  logratio_star <- loglik_star + logprior_star - logproposal_star
  logratio_current <- loglik_current + logprior_current - logproposal_current
  
  if(runif(1) < exp(logratio_star - logratio_current)){
    
    list_phi <- phi_new
    
  } else {
    
    list_phi <- phi_current
    
  }
  
  # phi_01_star <- optim(phi_01, 
  #                      fn = function(x){
  #                        - loglik_phi(x, X_phi, delta_w)
  #                      })$par
  # 
  # Hess_star <- der2_loglik_phi(phi_01_star, X_phi)
  # 
  # solve(- Hess_star)
  
  # phi_0 <- phi_0_true
  # phi_1 <- phi_1_true
  # phi_0_grid <- seq(phi_0_true - 1, phi_0_true + 1, length.out = 100)
  # y_grid <- sapply(phi_0_grid, function(x){
  #   loglik_phi(x, phi_1,
  #              X_phi,
  #              delta_w)
  # })
  # 
  # qplot(phi_0_grid, exp(y_grid) / sum(exp(y_grid))) + geom_vline(aes(xintercept = phi_0_true))
  # 
  # glm(delta_w ~ X_phi - 1, family = binomial(link = "logit"))
  # 
  # phi_0 <- phi_0_true
  # phi_1 <- phi_1_true
  # phi_1_grid <- seq(phi_1_true - 1, phi_1_true + 1, length.out = 100)
  # y_grid <- sapply(phi_1_grid, function(x){
  #   loglik_phi(phi_0, x,
  #              X_phi, 
  #              delta_w)
  # })
  # 
  # qplot(phi_1_grid, y_grid) + geom_vline(aes(xintercept = phi_1_true))
  
  # list_phi <- sample_betaPG(phi_01, X_phi, 
  #                                  mean_phi, diag(sigma_phi, nrow = 2),
  #                                  rep(1, nrow(X_phi)), k = delta_w - .5)
  
  # phi_0_star <- 
  
  
  
  list("phi_0" = list_phi[1],
       "phi_1" = list_phi[2])
  
}

updatePhi_PG <- function(Ct, C_star, w, w_star, 
                      delta, delta_star,
                      phi_0, phi_1, 
                      mean_phi, sigma_phi,
                      Sigma_prop){
  
  
  w_all <- c(w_pcr, w_star)
  delta_all <- c(delta, delta_star)
  
  phi_current <- c(phi_0, phi_1)
  
  idx_w <- which(w_all > 0)
  
  delta_w <- as.numeric(delta_all[idx_w] == 1)
  
  X_phi <- cbind(1, log(w_all[idx_w]))
  
  list_phi <- sample_betaPG(phi_current, X_phi,
                                   mean_phi, diag(sigma_phi, nrow = 2),
                                   rep(1, nrow(X_phi)), k = delta_w - .5)
  
  list("phi_0" = list_phi[1],
       "phi_1" = list_phi[2])
  
}

updateSigmaY <- function(Ct, w_pcr, delta, 
                         C_star, w_star, delta_star,
                         alpha1, alpha2,
                         sigma_y, sigma_alpha, 
                         P){
  
  w_all <- c(w_pcr, w_star)
  delta_all <- c(delta, delta_star)
  Ct_all <- c(Ct, C_star)
  P_all <- c(P, P_star)
  
  w_i <- w_all[delta_all == 1]
  Ct_i <- Ct_all[delta_all == 1]
  
  X_t <- cbind(1, log(w_i))
  alpha_i <- cbind(alpha1[P_all[delta_all == 1]], alpha2[P_all[delta_all == 1]])
  
  Xtalpha <- sapply(1:length(w_i), function(i){
    X_t[i,] %*% alpha_i[i,]
  })
  
  Ct_res <- (Ct_i - Xtalpha)^2
  
  sumsq <- sum(Ct_res)
  n_samples <- length(Ct_i)
  
  sigma_y <- sqrt(rinvgamma(a_sigma_y + n_samples / 2, b_sigma_y + sumsq / 2))
  
  sigma_y
  
}

updateAlpha <- function(Ct, w_pcr, delta, 
                        C_star, w_star, delta_star,
                        sigma_y, sigma_alpha, 
                        P, n_P, alpha_0){
  
  alpha <- matrix(NA, n_P, 2)
  
  w_all <- c(w_pcr, w_star)
  delta_all <- c(delta, delta_star)
  Ct_all <- c(Ct, C_star)
  P_all <- c(P, P_star)
  
  for (p in 1:n_P) {
    
    w_i <- w_all[delta_all == 1 & P_all == p]
    Ct_i <- Ct_all[delta_all == 1 & P_all == p]
    
    X_t <- cbind(1, log(w_i))
    
    Lambda_alpha <- (t(X_t) %*% X_t / sigma_y^2 + diag(1 / sigma_alpha^2, nrow = 2))
    
    mu_alpha <- matrix(as.numeric(sapply(seq_along(w_i), function(i){
      X_t[i,] * Ct_i[i]  / sigma_y^2
    })) , 2, length(w_i))
    # mu_alpha_star <- matrix(as.numeric(sapply(seq_along(w_tilde_i_star), function(i){
      # X_t[length(w_tilde_i) + i,] * C_tilde_i_star[i]  / sigma_y^2
    # })) , 2,length(w_tilde_i_star))
    # mu_alpha <- as.matrix(cbind(mu_alpha, mu_alpha_star))
    
    mu_alpha <- apply(mu_alpha, 1, sum) #+ matrix(alpha_0, 2, 1) / sigma_alpha^2
    
    mean_alpha <- solve(Lambda_alpha) %*% mu_alpha
    
    Cov_alpha <- solve(Lambda_alpha)
    
    alpha_p <- mvrnorm(1, mean_alpha, Cov_alpha)
    
    alpha[p,] <- c(alpha_p[1], alpha_p[2])
    
  }
  
  alpha
}

updateSigma <- function(l, v, X_w, beta_w, beta_w0,
                        a_sigma, b_sigma, idx_site){
  
  Xwbeta <- beta_w0 + X_w %*% beta_w
  
  l_samples <- l[idx_site]
  
  wtilde_res <- (v - (l_samples + Xwbeta ))^2
  
  wtilde_res_l0 <- wtilde_res[l_samples > 0]
  
  sumsq <- sum(wtilde_res_l0)
  n_samples <- length(wtilde_res_l0)
  
  sigma <- sqrt(rinvgamma(a_sigma + n_samples / 2, b_sigma + sumsq / 2))

  sigma
}

updatep0 <- function(delta, delta_star,
                     gamma, gamma_star,
                     a_p0, b_p0){
  
  delta_all <- c(delta, delta_star)
  gamma_all <- c(gamma, gamma_star)
  
  numTrials <- sum(delta_all == 0)
  
  numSuccesses <- sum(gamma_all == 1 & delta_all == 0)
  
  p0 <- rbeta(1, a_p0 + numSuccesses, 
              b_p0 + numTrials - numSuccesses)
  
  p0
}
