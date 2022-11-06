logistic <- function(x){
  1 / (1 + exp(-x))
}

erf <- function(x){
  ifelse(x > 0,
     pnorm(x, sd = 1 / sqrt(2)) - pnorm(-x, sd = 1 / sqrt(2)),
     - (pnorm(-x, sd = 1 / sqrt(2)) - pnorm(x, sd = 1 / sqrt(2))))
  
}

int_1 <- function(b, tau, mu){
  sqrt(pi)*tau*erf(sqrt(2)*b/(2*tau) - mu*sqrt(2)/(2*tau))/(2*sqrt(pi*tau^2))
}

int_2 <- function(x, mu, tau, a_b, b_b, c_b,M_i){
  sqrt(2)*sqrt(pi) * 
    exp(-mu^2/(2*tau^2) - c_b/(2*sigma^2) - (mu/tau^2 + b_b/sigma^2)^2/(4*(-1/(2*tau^2) - a_b/(2*sigma^2)))) * 
    erf(sqrt(2/tau^2 + 2*a_b/sigma^2)*x/2 - (mu/tau^2 + b_b/sigma^2)/sqrt(2/tau^2 + 2*a_b/sigma^2)) / 
    (2*sqrt(pi*tau^2)*(sqrt(2)*sqrt(pi*sigma^2))^M_i*sqrt(2/tau^2 + 2*a_b/sigma^2))
}

tnorm <- function(mu, sigma){
  
  x <- rnorm(1, mu, sigma)
  
  while(x > 0){
    x <- rnorm(1, mu, sigma)
  }
  
  x
}

updateBtilde <- function(wtilde, sigma, tau, X_b, beta_b, beta_b0, 
                         X_w, beta_w, beta_w0, nu, sigma_nu){
  
  n <- nrow(X_b)
  btilde <- rep(NA, n)
  Xbeta_b <- X_b %*% beta_b
  Xbeta_w <- X_w %*% beta_w
  
  for (i in 1:n) {
    
    wtilde_i <- wtilde[1:M[i] + sum(M[seq_len(i-1)])] - beta_w0 - (Xbeta_w)[1:M[i] + sum(M[seq_len(i-1)])]
    
    priorbtilde <- beta_b0 + Xbeta_b[i]
    p_0 <- prod(dnorm(wtilde_i, nu, sigma_nu)) * (int_1(0, tau, priorbtilde) - int_1(-1000, tau, priorbtilde))
    
    a_b <- M[i]
    b_b <- sum(wtilde_i)
    c_b <- sum(wtilde_i^2)
    p_1 <- int_2(1000, priorbtilde, tau, a_b, b_b, c_b, M[i]) - int_2(0, priorbtilde, tau, a_b, b_b, c_b, M[i])
    
    p_b0 <- p_1 / (p_1 + p_0)
    
    b0 <- rbinom(1, 1, p_b0) 
    
    if(b0 == 1){
      
      # b greater than 0
      
      prior_var <- 1 / tau^2 + M[i] / sigma^2
      
      sum_lik <- 0
      
      for(m in 1:M[i]){
        
        idx_m <- m + sum(M[seq_len(i-1)])
        sum_lik <- sum_lik + (wtilde[idx_m] - beta_0w - Xbeta_w[idx_m]) 
        
      }
      
      mu_btilde <- priorbtilde / tau^2 + sum_lik / sigma^2
      
      mean_btilde <- mu_btilde / prior_var
      
      btilde[i] <- rnorm(1, mean_btilde, sqrt(1 / prior_var))
      
        
    } else {
      
      # b less than 0
      btilde[i] <- tnorm(priorbtilde, tau)
      
    }
    
  }
  
  btilde
  
}

updateBetaB <- function(l, X_b, tau, sigma_beta){
  
  X_b1 <- cbind(1, X_b) 
  
  Lambda_beta <- t(X_b1) %*% X_b1 / tau^2 + diag(1 / sigma_beta^2, nrow = ncov_b + 1)
  
  mu_beta <- sapply(seq_along(l), function(i){
    X_b1[i,] * l[i] / tau^2
  })
  
  mu_beta <- apply(mu_beta, 1, sum)
  
  mean_beta <- solve(Lambda_beta) %*% mu_beta
  
  Cov_beta <- solve(Lambda_beta)
  
  beta0beta <- mvrnorm(1, mean_beta, Cov_beta)
  
  list("beta_b0" = beta0beta[1],
    "beta_b" = beta0beta[-1])
}

updateBetaW <- function(l, v, X_w, sigma, sigma_beta, idx_site,
                        isBetaW0){
  
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
  
  
  list("beta_w0" = beta_w0,
       "beta_w" = beta_w)
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
        
        log(dbinom(1, 1, 
               logistic(phi_0 + phi_1 * log(w_i)), 
               log = F) *
          dnorm(Ct_i[i], alpha1_i[i] + alpha2_i[i] * log(w_i), 
                sd = sigma_y) + 
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
        
        log(dbinom(1, 1, 
                   logistic(phi_0 + phi_1 * log(w_tilde_i)), 
                   log = F) *
              dnorm(Ct_i[i], alpha1_i[i] + alpha2_i[i] * log(w_tilde_i), 
                    sd = sigma_y, log = F) + 
              dbinom(0, 1, 
                     logistic(phi_0 + phi_1 * log(w_tilde_i)), 
                     log = F) *
              p0 * dnorm(Ct_i[i], lambda, sigma_lambda))
        
        
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

updateV <- function(v, Ct, l, beta_w0, X_w, beta_w,
                    alpha1, alpha2, P, w_tilde,
                    phi_0, phi_1, lambda, sigma_lambda){
  
  Xwbetaw <- beta_w0 + X_w %*% beta_w
  
  for (i in 1:n) {
    
    l_i <- l[i]
    
    for (m in 1:M[i]) {
      
      idx_current <- m + sum(M[seq_len(i-1)])
      v_mean <- l_i + Xwbetaw[idx_current]
      
      v_current <- v[idx_current]
      
      Ct_i <- Ct[idx_sample_K == idx_current]
      
      alpha1_i <- alpha1[P[idx_sample_K == idx_current]]
      alpha2_i <- alpha2[P[idx_sample_K == idx_current]]
      
      w_tilde_i <- w_tilde[idx_current]
      
      # v_grid <- seq(v_true[idx_current] - 3, v_true[idx_current] + 3, length.out = 100)
      # y_grid <- sapply(v_grid, function(x){
      #   logf_v(x, Ct_i, v_mean, sigma,
      #          alpha1_i, alpha2_i, w_tilde_i,
      #          sigma_y, p0, phi_0, phi_1,
      #          lambda, sigma_lambda)
      # })
      # 
      # qplot(v_grid, y_grid) +
      #   geom_vline(aes(xintercept = v_true[idx_current]))
      # qplot(v_grid, exp(y_grid) / sum(exp(y_grid))) +
      #   geom_vline(aes(xintercept = v_current))
      
      v_star <- rnorm(1, v_current, sd = .2)
      
      logposterior_current <- logf_v(v_current, Ct_i, v_mean, sigma, 
                                     alpha1_i, alpha2_i, w_tilde_i, 
                                     sigma_y, p0, phi_0, phi_1,
                                     lambda, sigma_lambda)
      
      logposterior_star <- logf_v(v_star, Ct_i, v_mean, sigma, 
                                     alpha1_i, alpha2_i, w_tilde_i, 
                                     sigma_y, p0, phi_0, phi_1,
                                     lambda, sigma_lambda)
      exp(logposterior_star - logposterior_current)
      if(runif(1) < exp(logposterior_star - logposterior_current)){
        
        v[idx_current] <- v_star
        
      }
        
        
      
      
    }
    
  }
  
  v
  
}

updateV_rjmcmc <- function(v, Ct, l, beta_w0, X_w, beta_w,
                           alpha1, alpha2, P, w_tilde,
                           phi_0, phi_1, lambda, sigma_lambda,
                           q_vtilde){
  
  Xwbetaw <- beta_w0 + X_w %*% beta_w
  
  for (i in 1:n) {
    
    l_i <- l[i]
    
    for (m in 1:M[i]) {
      
      idx_current <- m + sum(M[seq_len(i-1)])
      v_mean <- l_i + Xwbetaw[idx_current]
      
      v_current <- v[idx_current]
      
      Ct_i <- Ct[idx_sample_K == idx_current]
      
      alpha1_i <- alpha1[P[idx_sample_K == idx_current]]
      alpha2_i <- alpha2[P[idx_sample_K == idx_current]]
      
      w_tilde_i <- w_tilde[idx_current]
      
      # v_grid <- seq(v_true[idx_current] - 3, v_true[idx_current] + 3, length.out = 100)
      # y_grid <- sapply(v_grid, function(x){
      #   logf_v(x, Ct_i, v_mean, sigma,
      #          alpha1_i, alpha2_i, w_tilde_i,
      #          sigma_y, p0, phi_0, phi_1,
      #          lambda, sigma_lambda)
      # })
      # 
      # qplot(v_grid, y_grid) +
      #   geom_vline(aes(xintercept = v_true[idx_current]))
      # qplot(v_grid, exp(y_grid) / sum(exp(y_grid))) +
      #   geom_vline(aes(xintercept = v_current))
      
      v_star <- rnorm(1, v_current, sd = .2)
      
      # if(v_current > 0 & v_star < 0){
      #   
      #   w_tilde_i
      #   
      # }
      
      logposterior_current <- logf_v(v_current, Ct_i, v_mean, sigma, 
                                     alpha1_i, alpha2_i, w_tilde_i, 
                                     sigma_y, p0, phi_0, phi_1,
                                     lambda, sigma_lambda)
      
      logposterior_star <- logf_v(v_star, Ct_i, v_mean, sigma, 
                                     alpha1_i, alpha2_i, w_tilde_i, 
                                     sigma_y, p0, phi_0, phi_1,
                                     lambda, sigma_lambda)
      exp(logposterior_star - logposterior_current)
      if(runif(1) < exp(logposterior_star - logposterior_current)){
        
        v[idx_current] <- v_star
        
      }
        
        
      
      
    }
    
  }
  
  v
  
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
  
  Xwbetaw <- beta_w0 + X_w %*% beta_w
  
  for (i in 1:n) {
    
    l_i <- l[i]
    
    for (m in 1:M[i]) {
      
      idx_current <- m + sum(M[seq_len(i-1)])
      
      if(v[idx_current] > 0){
        
        v_tilde_current <- v_tilde[idx_current]
        
        Ct_i <- Ct[idx_sample_K == idx_current]
        
        alpha1_i <- alpha1[P[idx_sample_K == idx_current]]
        alpha2_i <- alpha2[P[idx_sample_K == idx_current]]
        
        # v_grid <- seq(v_current - 6, v_current + 5, length.out = 100)
        # y_grid <- sapply(v_grid, function(x){
        #   logf_v(x, Ct_i, v_mean, sigma,
        #          alpha1_i, alpha2_i, w_tilde_i,
        #          sigma_y, p0, phi_0, phi_1,
        #          lambda, sigma_lambda)
        # })
        # 
        # qplot(v_grid, y_grid) +
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
      
      p_delta1 <- exp(log_delta1 - log_gamma1) / (1 + exp(log_delta1 - log_gamma1))
      
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

updatePhi <- function(Ct, C_star, w, w_star, 
                      delta, delta_star,
                      phi_0, phi_1, 
                      mean_phi, sigma_phi,
                      Sigma_prop){
  
  
  w_all <- c(w_pcr, w_star)
  delta_all <- c(delta, delta_star)
  
  phi_current <- c(phi_0, phi_1)
  
  idx_w <- which(w_all > 0)
  
  delta_w <- delta_all[idx_w] == 1
  
  X_phi <- cbind(1, log(w_all[idx_w]))
  
  # phi_01_star <- mvrnorm(1, c(phi_0, phi_1), Sigma = Sigma_prop)
  
  phi_star <- optim(phi_current, fn = function(phi){
    - loglik_phi_cpp(phi, X_phi[,2], delta_w)
  }, method = "BFGS", control = list(reltol = .001))$par
  
  Hess_phi_star <- hes_loglik_phi_cpp(phi_star, X_phi[,2], delta_w)
  Sigma_star <- solve(- Hess_phi_star)
  
  phi_new <- mvrnorm(mu = phi_star, Sigma = Sigma_star)
  
  logproposal_current <- dmvnorm_cpp(phi_current, phi_star, Sigma_star, 1)
  logproposal_star <- dmvnorm_cpp(phi_new, phi_star, Sigma_star, 1)
  
  loglik_star <- loglik_phi(phi_new, X_phi, delta_w)
  loglik_current <- loglik_phi(phi_current, X_phi, delta_w)
  
  logratio_star <- loglik_star - logproposal_star
  logratio_current <- loglik_current - logproposal_current
  
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

