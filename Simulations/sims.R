
setwd("~/qPCR_2/Simulations")

beta_b0_settings <- c(0, 1, 2)

for(idxdataSettings in seq_along(beta_b0_settings)){
  
  beta_b0_current <- beta_b0_settings[idxdataSettings]
  
  {
    n = 100
    M = 1 + rep(3, n)
    K = 1 + rep(3, sum(M))
    n_P = n
    idx_site <- rep(1:n, M)
    idx_site_K <- rep(idx_site, K)
  }
  
  dataParams <- list(
    n = n,
    M = M,
    K = K,
    n_P = n,
    P = idx_site_K,
    ncov_b = 2,
    ncov_w = 2,
    ncov_wt = 1,
    K_standards = 3,
    qty = c(3e+07, 3e+06, 3e+05, 3e+04, 3e+03, 3e+02, 3e+01),
    numCycles = 45
  )
  
  hyperParams <- list(
    tau = 1,
    sigma = 2,
    sigma_y = 1,
    beta_b0 = beta_b0_current,
    beta_w0 = 1,
    theta0 = .02,
    p0 = .05,
    nu = -1,
    sigma_nu = 1,
    alpha1_0 = 45,
    alpha2_0 = -1.7,
    lambda = 20,
    sigma_lambda = 10,
    phi_0 = -26.643,
    phi_1 = 8.3
  )
  
  list_out <- simulateData(dataParams, hyperParams)
  
  dataParams <- list_out$dataParams
  trueParams <- list_out$trueParams
  
  list2env(dataParams,globalenv())
  list2env(trueParams,globalenv())
  
  sites <- unique(data_real$Site)
  n <- length(sites)
  M <- sapply(sites, function(x){
    length(unique(data_real$Sample[data_real$Site == x]))
  })
  
  samples <- unique(data_real$Sample)
  K <- sapply(samples, function(x){
    sum(data_real$Sample == x)
  })
  idx_sample_K <- rep(1:sum(M), K)
  
  idx_site <- rep(1:n, M)
  
  Ct <- data_real$Ct
  
  P <- data_real$Plate
  
  # w_star <- data_standard$Qty
  # C_star <- data_standard$Ct
  # P_star <- data_standard$P
  w_star <- data_standard$Qty
  C_star <- data_standard$Ct
  P_star <- data_standard$P
  # Ct ordered by sample!!
  numSampleV <- as.numeric(table(idx_sample_K))
  
  ncov_b <- ncol(X_b)
  ncov_w <- ncol(X_w)
  ncov_wt <- ncol(X_wt)
  
  n_P <- length(unique(data_real$Plate))
  
  
  # MCMC  -----------
  
  trueStartingValues <- F
  
  nchain <- 1
  niter <- 1000
  nburn <- 1000
  nthin <- 1
  
  # output
  {
    l_output <- array(NA, dim = c(nchain, niter, n) )
    beta_output <- array(NA, dim = c(nchain, niter, ncov_b + 1))
    tau_output <- matrix(NA, nchain, niter)
    v1_output <- array(NA, dim = c(nchain, niter, sum(M)) )
    v2_output <- array(NA, dim = c(nchain, niter, sum(M)) )
    delta_output <- array(NA, dim = c(nchain, niter, sum(M)) )
    gamma_output <- array(NA, dim = c(nchain, niter, sum(M)) )
    beta_w_output <- array(NA, c(nchain, niter, ncov_w + 1))
    nu_output <- array(NA, c(nchain, niter))
    beta_wt_output <- array(NA, c(nchain, niter, ncov_wt))
    sigma_output <- matrix(NA, nchain, niter)
    alpha_output <- array(NA, dim = c(nchain, niter, n_P, 2))
    phi_output <- array(NA, dim = c(nchain, niter, 2))
    p0_output <- array(NA, dim = c(nchain, niter))
    sigma_y_output <- array(NA, dim = c(nchain,  niter, n_P))
  }
  
  # params to update
  {
    updBetaB <- T
    updTau <- T
    updL <- T
    updV1 <- T
    updV2 <- T
    updV1V2 <- T
    updBetaW <- T
    updBetaWt <- T
    updSigma <- T
    updSigmaNu <- T
    updDelta <- T
    updPhi <- T
    updp0 <- T
    updSigmaY <- T
    updAlpha <- T
  }
  
  for (chain in 1:nchain) {
    
    # chain output
    {
      l_output_chain <- matrix(NA, niter, n) 
      beta_output_chain <- matrix(NA, niter, ncov_b + 1)
      tau_output_chain <- rep(NA, niter)
      v1_output_chain <- matrix(NA, niter, sum(M)) 
      v2_output_chain <- matrix(NA, niter, sum(M)) 
      delta_output <- array(NA, dim = c(niter, sum(M)) )
      gamma_output <- array(NA, dim = c(niter, sum(M)) )
      beta_w_output_chain <- matrix(NA, niter, ncov_w + 1)
      nu_output_chain <- rep(NA, niter)
      beta_wt_output_chain <- matrix(NA, niter, ncov_wt)
      sigma_output_chain <- rep(NA, niter)
      alpha_output_chain <- array(NA, dim = c(niter, n_P, 2))
      phi_output_chain <- array(NA, dim = c(niter, 2))
      p0_output_chain <- rep(NA, niter)
      sigma_y_output_chain <- matrix(NA, niter, n_P)
      
      v1_values <- matrix(NA, niter + nburn, sum(M)) 
      v2_values <- matrix(NA, niter + nburn, sum(M)) 
      phi_values <- matrix(NA, niter + nburn, 2) 
    }
    
    # starting values
    {
      
      if(trueStartingValues){
        
        l <- l_true
        beta_b <- beta_b_true
        beta_b0 <- beta_b0_true
        tau <- tau_true
        
        beta_w0 <- beta_w0_true
        beta_w <- beta_w_true
        
        v1 <- v1_true
        v2 <- v2_true
        sigma <- sigma_true
        
        phi_0 <- phi_0_true
        phi_1 <- phi_1_true
        
        delta <- as.numeric(delta_true)
        gamma <- as.numeric(gamma_true)
        delta_star <- as.numeric(delta_star_true)
        gamma_star <- as.numeric(delta_star_true)
        
        nu <- nu_true
        sigma_nu <- sigma_nu_true
        
        p0 <- p0_true
        
        sigma_y <- sigma_y_true
        alpha <- cbind(alpha1_true, alpha2_true)
        
      } else {
        
        l <- rep(10, n)
        beta_b <- rep(0, ncov_b)
        tau <- b_tau / (a_tau - 1)
        
        beta_w0 <- 0
        beta_w <- rep(0, ncov_w)
        
        beta_wt <- rep(0, ncov_wt)
        
        v1 <- rep(10, sum(M))
        v2 <- rep(0, sum(M))
        sigma <- b_sigma / (a_sigma - 1)
        
        nu <- nu0
        sigma_nu <- b_sigma_nu / (a_sigma_nu - 1)
        
        phi_0 <- mean_phi[1]
        phi_1 <- mean_phi[2]
        
        delta <- Ct > 0
        gamma <- rep(0, length(delta))
        delta_star <- C_star > 0
        gamma_star <- rep(0, length(delta_star))
        
        p0 <- a_p0 / (a_p0 + b_p0)
        
        sigma_y <- rep(b_sigma_y / (a_sigma_y - 1), n_P)
        alpha <- matrix(alpha_0, n_P, 2, byrow = T)
        alpha1 <- alpha[,1]
        alpha2 <- alpha[,2]
      }
      
      beta_w0 <- 0
      
      w1 <- convertPositive(v1)
      w2 <- convertPositive(v2)
      # w1 <- (v1_tilde > 0) * (exp(v) - 1)
      # w2 <- (v_tilde > 0) * (exp(v_tilde) - 1)
      
      w <- w1 + w2
      w_pcr <- w[idx_sample_K]
      # v1 <- ifelse(v1_tilde > a, v1_tilde, NA)
      # v2 <- ifelse(v2_tilde > a, v2_tilde, NA)
      
      # v_tilde <- ifelse(v1_tilde > a, v1_tilde, v2_tilde)
      # v <- ifelse(!is.na(v1_tilde), v1, v2)
      
      v2_pcr <- v2[idx_sample_K]
      v1_pcr <- v1[idx_sample_K]
      v1_pcr <- ifelse(!is.na(v1_pcr), v1_pcr, -Inf)
      # Ctilde <- data_real$Ct
      # Ctilde[Ctilde > numCycles] <- Ctilde[Ctilde > numCycles] + rnorm(length(Ctilde > numCycles))
      
      # C_tilde_star[C_tilde_star > numCycles] <- C_tilde_star[C_tilde_star > numCycles] + rnorm(length(C_tilde_star > numCycles))
      
      # alpha_0 <- c(45, -1.5)
      # P_all <- c(P, P_star)
      
      # w_tilde_all <- c(w_tilde, w_tilde_star)
      
      # wtilde <-
      # gamma <- c(Ctilde[idx_real] < numCycles, rep(1, length(idx_standards)))
      # idx_real <- 1:length(data_real$Ct)
      # idx_standards <- length(data_real$Ct) + 1:length(data_standard$Ct)
      
      # tau <- 1
      # sigma <- 1
      # sigma_y <- .05
      
      sigmas_v1_1 <- rep(sigma_v1, sum(M))
      sigmas_v1_2 <- rep(sigma_v2, sum(M)) # larger set of proposal to allow larger moves
      sigmas_v2 <- rep(sigma_v1, sum(M))
      sigmas_vv2 <- rep(sigma_v_v2, sum(M))
      
      Sigma_prop <- Sigma_prop_0
    }
    
    for (iter in 1:(nburn + nthin * niter)) {
      
      if(iter %% 100 == 0){
        if(iter <= nburn){
          print(paste0("Chain = ",chain," - Burn-in Iteration = ",iter))
        } else {
          print(paste0("Chain = ",chain," - Iteration = ",iter - nburn))
        }  
        
      }
      
      # UPDATE BETA B -------
      
      if(updBetaB){
        list_beta <- updateBetaB(l, X_b, tau, 
                                 prior_betab0, 
                                 sigma_beta0, sigma_beta)
        beta_b0 <- list_beta$beta_b0
        beta_b <- list_beta$beta_b
      }
      
      # UPDATE TAU -----------
      
      if(updTau){
        tau <- updateTau(l, X_b, beta_b0, 
                         beta_b, a_tau, b_tau)
      }
      
      # UPDATE L -------
      
      if(updL){
        l <- updateL(v1, X_b, beta_b0, beta_b,
                     beta_w0, X_w, beta_w, tau)
      }
      
      # UPDATE V ----------
      
      if(updV1){
        list_v <- updateV_cpp(v1, v2, Ct, 
                              l, M, beta_w0, X_w, beta_w,
                              sigma, nu, X_wt, beta_wt, 
                              sigma_nu,
                              alpha1, alpha2, P, w2,
                              numSampleV, sigma_y, p0,
                              phi_0, phi_1, sigmas_v1_1,
                              lambda, sigma_lambda)
        # list_v <- updateV(v, v_tilde, Ct, l, beta_w0, X_w, beta_w,
        #                   alpha1, alpha2, P,
        #                   nu, X_wt, beta_wt, sigma_nu, sigmas_v,
        #                   phi_0, phi_1, lambda, sigma_lambda)
        v1 <- list_v$v
        v2 <- list_v$v_tilde
        
        list_v <- updateV_cpp(v1, v2, Ct, l, M, beta_w0, 
                              X_w, beta_w,
                              sigma, nu, X_wt, beta_wt, sigma_nu,
                              alpha1, alpha2, P, w2,
                              numSampleV, sigma_y, p0,
                              phi_0, phi_1, sigmas_v1_2,
                              lambda, sigma_lambda)
        # list_v <- updateV(v, v_tilde, Ct, l, beta_w0, X_w, beta_w,
        #                   alpha1, alpha2, P,
        #                   nu, X_wt, beta_wt, sigma_nu, sigmas_v2,
        #                   phi_0, phi_1, lambda, sigma_lambda)
        v1 <- list_v$v
        v2 <- list_v$v_tilde
        
        w1 <- convertPositive(v1) # (v > 0) * (exp(v) - 1)
        w2 <- convertPositive(v2) # (v > 0) * (exp(v) - 1)
        
        w <- w1 + w2
        
        w_pcr <- w[idx_sample_K]
        
        v1_values[iter,] <- v1
        
        if(iter %% 100 == 0){
          
          hat_sigmas_v <- apply(v1_values[1:(iter - 1),], 2, var)
          
          sigmas_v1_1 <- sqrt(proposal_ratio * hat_sigmas_v + 
                                (1 - proposal_ratio) * sigma_v1^2)
          
        }
        
      }
      
      # UPDATE V TILDE ----
      
      if(updV2){
        v2 <- updateVtilde_cpp(v2, v1, Ct, n, M,
                               alpha1, alpha2, P, nu, X_wt, 
                               beta_wt, sigma_nu,
                               numSampleV, sigma_y, p0,
                               phi_0, phi_1, sigmas_v2, 
                               lambda, sigma_lambda)
        
        w2 <- convertPositive(v2)
        # w_tilde <- (v_tilde > 0) * (exp(v_tilde) - 1)
        w2[is.na(v2)] <- 0
        
        w <- w1 + w2
        
        w_pcr <- w[idx_sample_K]
        
        v2_values[iter,] <- v2
        
        if(iter %% 100 == 0){
          
          hat_sigmas_v <- apply(v2_values[1:(iter - 1),], 2, var)
          
          sigmas_v2 <- sqrt(proposal_ratio * hat_sigmas_v +
                              (1 - proposal_ratio) * sigma_v1^2)
          
        }
        
      }
      
      # UPDATE V VTILDE -------
      
      if(updV1V2){
        list_v1v2 <- updateV_Vtilde_cpp(v1, v2, Ct, 
                                        l, M, beta_w0, X_w, 
                                        beta_w, sigma, 
                                        nu, sigma_nu, X_wt, 
                                        beta_wt, alpha1, alpha2, 
                                        P, numSampleV, 
                                        sigma_y, sigmas_vv2, 
                                        p0, phi_0, phi_1, 
                                        lambda, sigma_lambda)
        # list_vvtilde <- updateV_Vtilde(v, v_tilde, Ct, l, beta_w0, X_w, beta_w,
        #                                alpha1, alpha2, P, nu, sigma_nu, sigmas_vvtilde,
        #                                phi_0, phi_1, lambda, sigma_lambda)
        v1 <- list_v1v2$v
        v2 <- list_v1v2$v_tilde
        
        w <- w1 + w2
        
        w_pcr <- w[idx_sample_K]
      }
      
      # UPDATE DELTA GAMMA ----------
      
      # list_deltagamma <- updateDeltaGamma(Ct, C_star, w_pcr, w_star,
      #                                     phi_0, phi_1, p0)
      # delta <- list_deltagamma$delta
      # gamma <- list_deltagamma$gamma
      # delta_star <- list_deltagamma$delta_star
      # gamma_star <- list_deltagamma$gamma_star
      
      if(updDelta){
        list_deltagamma <- updateDeltaGammaCpp(Ct, C_star, w_pcr, w_star,
                                               phi_0, phi_1, p0, alpha1,
                                               alpha2, sigma_y, P, 
                                               P_star, lambda,
                                               sigma_lambda)
        delta <- list_deltagamma$delta
        gamma <- list_deltagamma$gamma
        delta_star <- list_deltagamma$delta_star
        gamma_star <- list_deltagamma$gamma_star
        
        # if(iter > 1 & any(C_star > 0 & delta_star != 1)){
        #   browser()
        # }
      }
      
      # UPDATE SIGMA -----------
      
      if(updSigma){
        sigma <- updateSigma(l, v1, X_w, beta_w, beta_w0,
                             a_sigma, b_sigma, idx_site)
        
      }
      
      # UPDATE BETA W -------
      
      if(updBetaW){
        list_betaw <- updateBetaW(l, v1, X_w, sigma, sigma_beta, 
                                  idx_site, isBetaW0 = T)
        beta_w0 <- list_betaw$beta_w0
        beta_w <- list_betaw$beta_w
      }
      
      # UPDATE BETA Wt -------
      
      if(updBetaWt){
        list_betawt <- updateBetaWt(v2, v1, X_wt, 
                                    nu0, sigma_nu, sigma_beta, 
                                    ncov_wt)
        nu <- list_betawt$nu
        beta_wt <- list_betawt$beta_wt
      }
      
      # UPDATE SIGMA NU -----------
      
      if(updSigmaNu){
        sigma_nu <- updateSigmaNu(v2, v1, X_wt, nu, 
                                  a_sigma_nu, b_sigma_nu)
      }
      
      # UPDATE PHI --------
      
      if(updPhi){
        # list_phi <- updatePhi(Ct, C_star, w, w_star,
        #                       delta, delta_star,
        #                       phi_0, phi_1,
        #                       mean_phi, sigma_phi,
        #                       Sigma_prop)
        list_phi <- updatePhi_PG(Ct, C_star, w, w_star,
                                 delta, delta_star,
                                 phi_0, phi_1,
                                 mean_phi, sigma_phi,
                                 Sigma_prop)
        phi_0 <- list_phi$phi_0
        phi_1 <- list_phi$phi_1
        
        phi_values[iter,] <- c(phi_0, phi_1)
        
        if(iter %% 100 == 0){
          
          hat_Sigma_phi <- cov(phi_values[1:(iter - 1),])
          
          #2.4 / sqrt(2)
          Sigma_prop <- proposal_ratio * hat_Sigma_phi + (1 - proposal_ratio) * Sigma_prop_0
          
        }
      }
      
      # UPDATE p0 --------
      
      if(updp0){
        
        p0 <- updatep0(delta, delta_star,
                       gamma, gamma_star,
                       a_p0, b_p0)
        
      }
      
      # UPDATE SIGMA Y ---------
      
      if(updSigmaY){
        # sigma_y <- updateSigmaY(Ct, w_pcr, delta,
        #                         C_star, w_star, delta_star,
        #                         alpha1, alpha2,
        #                         sigma_y, sigma_alpha,
        #                         P)
        
        sigma_y <- updateSigmaY_cpp(Ct, w_pcr, delta, C_star, w_star,
                                    delta_star, alpha1,
                                    alpha2, sigma_alpha, P, P_star, n_P, 
                                    a_sigma_y, b_sigma_y)
        
      }
      
      # UPDATE ALPHA ----------
      
      if(updAlpha){
        # alpha <- updateAlpha(Ct, w_pcr, delta, 
        #                      C_star, w_star, delta_star,
        #                      sigma_y, sigma_alpha, 
        #                      P, n_P, alpha_0)
        
        alpha <- updateAlpha_cpp(Ct, w_pcr, delta, P, 
                                 C_star, w_star, delta_star, P_star, 
                                 sigma_y, sigma_alpha, n_P, alpha_0)
        alpha1 <- alpha[,1]
        alpha2 <- alpha[,2]
        
      }
      
      # OUTPUT -------
      
      {
        if(iter > nburn){
          trueIter <- iter - nburn
          l_output_chain[trueIter,] <- l
          beta_output_chain[trueIter,] <- c(beta_b0, beta_b)
          tau_output_chain[trueIter] <- tau  
          v1_output_chain[trueIter,] <- v1
          v2_output_chain[trueIter,] <- v2
          sigma_output_chain[trueIter] <- sigma
          beta_w_output_chain[trueIter,] <- c(beta_w0, beta_w)
          nu_output_chain[trueIter] <- nu
          beta_wt_output_chain[trueIter,] <- beta_wt
          alpha_output_chain[trueIter,,] <- alpha
          phi_output_chain[trueIter,] <- c(phi_0, phi_1)
          p0_output_chain[trueIter] <- p0
          sigma_y_output_chain[trueIter,] <- sigma_y
          
        }
        
      }
      
    }
    
    {
      l_output[chain,,] <- l_output_chain
      beta_output[chain,,] <- beta_output_chain
      tau_output[chain,] <- tau_output_chain
      v1_output[chain,,] <- v1_output_chain
      v2_output[chain,,] <- v2_output_chain
      beta_w_output[chain,,] <- beta_w_output_chain
      nu_output[chain,] <- nu_output_chain
      beta_wt_output[chain,,] <- beta_wt_output_chain
      sigma_output[chain,] <- sigma_output_chain
      alpha_output[chain,,,] <- alpha_output_chain
      phi_output[chain,,] <- phi_output_chain
      p0_output[chain,] <- p0_output_chain
      sigma_y_output[chain,,] <- sigma_y_output_chain
    }
    
  }
  
  # generate plots 
  
  beta_b_plot <- checkBeta(beta_b_true,
            beta_output,
            dropFirst = T)
  
  beta_w_plot <- checkBeta(beta_w_true,
            beta_w_output,
            dropFirst = T)
  
  beta_wt_plot <- checkBeta(beta_wt_true,
            beta_wt_output,
            dropFirst = F)
  
  var_plot <- checkVariances(tau_output,
                             tau_true,
                             sigma_output,
                             sigma_true,
                             sigma_y_output,
                             sigma_y_true)
  
  v1_plot <- checkElements(v1_output,
                           v1_true,
                           c(1 + 1:5))
  
  l_plot <- checkElements(l_output,
                           l_true,
                           c(1 + 1:5))
  
  ggsave(paste0("beta_b_plot_",idxdataSettings,".png"), beta_b_plot)
  ggsave(paste0("beta_w_plot",idxdataSettings,".png"), beta_w_plot)
  ggsave(paste0("beta_wt_plot",idxdataSettings,".png"), beta_wt_plot)
  ggsave(paste0("var_plot",idxdataSettings,".png"), var_plot)
  
}
