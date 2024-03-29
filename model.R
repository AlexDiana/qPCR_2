library(MASS); library(ars); library(ggplot2)
library(Rcpp); library(RcppArmadillo)
library(here)
sourceCpp(here("Functions/functions.cpp"))
source("~/qPCR_2/functions.R", echo=TRUE)
source("~/qPCR_2/diagnostics.R", echo=TRUE)

# prior -----

sigma_alpha <- 10
prior_betab0 <- 5;
sigma_beta0 <- 2; sigma_beta <- 1

a_tau <- 2; b_tau <- 2
a_sigma <- 2; b_sigma <- 1


a_sigma_y <- 2
b_sigma_y <- (.25) * (a_sigma_y - 1)

mean_phi <- c(-26, 8)
sigma_phi <- 3

alpha_0 <- c(45, -1.7)

nu0 <- -3
a_sigma_nu <- 2; b_sigma_nu <- 2

a_p0 <- 1; b_p0 <- 10

lambda <- 20
sigma_lambda <- 10

# proposal values
{
  sigma_v1 <- .5
  sigma_v2 <- 5
  sigma_v_v2 <- .01
  proposal_ratio <- .95
  Sigma_prop_0 <- diag(c(8, 1), nrow = 2)
}


# LOAD DATA ----

realData <- T

if(realData) {
  load("~/qPCR_2/Data/sara.rda")
  samples <- unique(data_real$Sample)
  X_w <- as.matrix(X_w[match(X_w[,1], samples),-1])
  sites <- unique(data_real$Site)
  X_b <- as.matrix(X_b[match(X_b[,1], sites),-1])
  X_wt <- matrix(0, nrow(X_w), 0)
}

# DATA CLEANING --------

numCycles <- 45

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
if(realData){
  w_star <- data_standard$Quantity
  C_star <- data_standard$C.
  P_star <- data_standard$Plate
}
# Ct ordered by sample!!
numSampleV <- as.numeric(table(idx_sample_K))

ncov_b <- ncol(X_b)
ncov_w <- ncol(X_w)
ncov_wt <- ncol(X_wt)

n_P <- length(unique(data_real$Plate))

# MCMC  -----------

trueStartingValues <- F

nchain <- 1
niter <- 5000
nburn <- 5000
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
      if(ncov_wt > 0){
        list_betawt <- updateBetaWt(v2, v1, X_wt, 
                                    nu0, sigma_nu, sigma_beta, 
                                    ncov_wt)
        nu <- list_betawt$nu
        beta_wt <- list_betawt$beta_wt  
      }
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

# OUTPUT TRUE ------

setwd("~/qPCR_2/Output")

library(coda)

# ess
{
  diagnosticsCheck <- function(params_output){
    
    dims <- dim(params_output)
    
    if(length(dims) > 2){
      params_output <- apply(params_output, setdiff(1:length(dims), 1:2), c)
      
      dims <- dim(params_output)
      dims_margin <- setdiff(1:length(dims), 1)
      return(apply(params_output, dims_margin, function(x){
        
        if(any(!is.na(x))){
          x <- x[!is.na(x)]
          x_current <- mcmc(x)
          return(as.numeric(effectiveSize(x_current)))
          
        } else {
          return(NA)
        }
      }))
      
    } else {
      
      # wrong
      x_current <- mcmc(params_output[1,])
      return(effectiveSize(x_current))
      
    }
    
    
  }
  
  ess_l <- diagnosticsCheck(l_output)
  print(min(ess_l))
  ess_beta <- diagnosticsCheck(beta_output)
  print(min(ess_beta))
  ess_beta_w <- diagnosticsCheck(beta_w_output)
  print(min(ess_beta_w[-1]))
  ess_alpha <- diagnosticsCheck(alpha_output)
  print(min(ess_alpha))
  ess_v1 <- diagnosticsCheck(v1_output)
  print(min(ess_v1))
  ess_v2 <- diagnosticsCheck(v2_output)
  print(min(ess_v2, na.rm = T))
  ess_nu <- diagnosticsCheck(nu_output)
  print(min(ess_nu))
  ess_beta_wt <- diagnosticsCheck(beta_wt_output)
  print(min(ess_beta_wt))
  ess_tau <-  diagnosticsCheck(tau_output)
  print(min(ess_tau))
  ess_phi <- diagnosticsCheck(phi_output)
  print(min(ess_phi))
  ess_p0 <- diagnosticsCheck(p0_output)
  print(min(ess_p0))
}

# traceplot
{
  qplot(1:niter, beta_w_output[,,3])
  qplot(1:niter, beta_output[,,2])
  qplot(1:niter, beta_w_output[,,4])
  qplot(1:niter, p0_output[1,])
  qplot(1:niter, alpha_output[1,,10,2])
  qplot(1:niter, v_output[1,,1])
  qplot(1:niter, tau_output[1,])
  qplot(1:niter, sigma_output[1,])
}

# covariates
{
  covariatePlot <- function(beta_output, covNames){
    
    beta_CI <- apply(beta_output[,,-1,drop = F], 3, function(x){
      quantile(x, probs = c(.025,.5,.975))
    })
    beta_mean <- beta_CI[2,]
    colnames(beta_CI) <- covNames
    
    idxCovs <- order(beta_mean)
    
    # covnames_order <- covnmaes[idxCovs]
    # 
    # OTUnames <- as.character(1:length(OTUnames))
    # factorSub <- factor(idxSpecies, levels = idxSpecies)
    # namesSpecies <- OTUnames
    
    data_plot <- data.frame(Covariates = covNames,
                            y = beta_CI[2,],
                            ymin = beta_CI[1,],
                            ymax = beta_CI[3,])
    
    (covplot_all <- 
        ggplot2::ggplot(data = data_plot, ggplot2::aes(x = Covariates,
                                                       y = y,
                                                       ymin = ymin,
                                                       ymax = ymax)) + 
        geom_errorbar() + 
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 20),
                       axis.title = ggplot2::element_text(size = 20, face = "bold"),
                       axis.text = ggplot2::element_text(size = 13, face = "bold", angle = 0),
                       panel.grid.major = ggplot2::element_line(colour="grey", size=0.15),
                       legend.title = element_text(size=14), #change legend title font size
                       legend.text = element_text(size=10),
                       panel.background = ggplot2::element_rect(fill = "white", color = "black")) +
        ggplot2::geom_hline(aes(yintercept = 0), color = "red") + 
        ylab("Value") + 
        ggplot2::scale_y_continuous(breaks = (-10):10) + ggplot2::coord_flip())
    
    
    
  }
  
  cov_b_plot <- covariatePlot(beta_output, colnames(X_b))
  
  ggsave("Cov_concentration_plot.jpeg", cov_b_plot)
  
  cov_w_plot <- covariatePlot(beta_w_output, colnames(X_w))
  
  ggsave("Cov_detection_plot.jpeg", cov_w_plot)
  
}

# site plot
{
  sitetype <- ifelse(X_b[,"typelake"] == 1, "lake",
                     ifelse(X_b[,"typereservoir"] == 1, "reservoir",
                            ifelse(X_b[,"typeriver"] == 1, "river","pond")))
  
  sitesPlot <- function(l_output, sites, sitetype){
    
    l_CI <- apply(l_output, 3, function(x){
      quantile(x, probs = c(.025,.5,.975))
    })
    l_mean <- l_CI[2,]
    colnames(l_CI) <- sites
    
    idxCovs <- order(l_mean)
    
    # covnames_order <- covnmaes[idxCovs]
    # 
    # OTUnames <- as.character(1:length(OTUnames))
    # factorSub <- factor(idxSpecies, levels = idxSpecies)
    # namesSpecies <- OTUnames
    
    data_plot <- data.frame(Sites = sites,
                            Type = sitetype,
                            y = l_CI[2,],
                            ymin = l_CI[1,],
                            ymax = l_CI[3,])
    
    (sitesplot_all <- 
        ggplot2::ggplot(data = data_plot, ggplot2::aes(x = Sites,
                                                       y = y,
                                                       ymin = ymin,
                                                       ymax = ymax)) + 
        geom_errorbar() + facet_grid(rows = vars(Type), scale = "free") + 
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 20),
                       axis.title = ggplot2::element_text(size = 20, face = "bold"),
                       axis.text = ggplot2::element_text(size = 13, face = "bold", angle = 0),
                       panel.grid.major = ggplot2::element_line(colour="grey", size=0.15),
                       legend.title = element_text(size=14), #change legend title font size
                       legend.text = element_text(size=10),
                       panel.background = ggplot2::element_rect(fill = "white", color = "black")) +
        # ggplot2::geom_hline(aes(yintercept = 0), color = "red") + 
        ylab("Site concentration") + ggplot2::coord_flip()
      # ggplot2::scale_y_continuous(breaks = (-10):10) + ggplot2::coord_flip()
    )
    
    
  }
  
  sitesplot <- sitesPlot(l_output, sites, sitetype)
  
  ggsave("Sites_concentration_plot.jpeg", sitesplot)
}


# diagnostics plot ---------

qplot(1:niter, l_output[,,1])

qplot(1:niter, beta_output_chain[,1])
qplot(1:niter, beta_w_output[,,3])
qplot(1:niter, beta_output_chain[,3])

# DIAGNOSTICS PLOT -----

tracePlotParameters <- function(param_output, idx = NULL) {
  
  nchain <- dim(param_output)[1]
  niter <- dim(param_output)[2]
  
  if(!is.null(idx)){
    if(length(idx) == 1){
      param_output <- param_output[,,idx,drop=F]
      param_output <-
        matrix(param_output[, , 1], nrow = nchain, ncol = niter)
    } else {
      param_output <- param_output[,,idx[1],idx[2],drop=F]
      param_output <-
        matrix(param_output[, , 1, 1], nrow = nchain, ncol = niter)
    }
  }
  
  param_output_long <- reshape2::melt(param_output)
  
  diagnosticsPlot <-
    ggplot2::ggplot(data = param_output_long, ggplot2::aes(
      x = Var2,
      y = value,
      group = Var1,
      color = factor(Var1)
    )) + ggplot2::geom_line() +
    ggplot2::xlab("Iterations") + ggplot2::ylab("Value") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 17),
      axis.title = ggplot2::element_text(size = 16, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 11, face = "bold"),
      axis.text.x = ggplot2::element_text(
        size = 11,
        face = "bold",
        hjust = 1
      ),
      axis.line = ggplot2::element_line(colour = "black", size = 0.15),
      # panel.grid.minor = element_line(colour="grey", size=0.15),
      panel.grid.major = ggplot2::element_line(colour = "grey", size = 0.15),
      panel.background = ggplot2::element_rect(fill = "white", color = "black"),
      legend.position = "none"
    )
  
  return(diagnosticsPlot)
  
}

# OUTPUT ------

idx_l <- 61
qplot(1:niter, v_output_chain[,idx_l])  + 
  geom_hline(aes(yintercept = v_true[idx_l]))

qplot(1:niter, v_tilde_output_chain[,1])  + geom_hline(aes(yintercept = vtilde_true[1]))

l <- 3
qplot(1:niter, v_output_chain[,l])  + geom_hline(aes(yintercept = v_true[l]))

qplot(1:niter, l_output_chain[,1], geom = "line")  + geom_hline(aes(yintercept = l_true[1])) #+ 
qplot(1:niter, l_output_chain[,1], geom = "line")  + geom_hline(aes(yintercept = l_true[1])) #+ 
  # ylim(c(-2,2))

# beta b

qplot(1:niter, beta_output_chain[,1]) + geom_hline(aes(yintercept = beta_0b_true))

l <- 2
qplot(1:niter, beta_output_chain[,l + 1], geom = "line") + 
  geom_hline(aes(yintercept = beta_b_true[l])) + ylim(c(-2,2))

# beta w
qplot(1:niter, beta_w_output_chain[,1], geom = "line") + 
  geom_hline(aes(yintercept = beta_w0_true))

l <- 2
qplot(1:niter, beta_w_output_chain[,l + 1], geom = "line") + 
  geom_hline(aes(yintercept = beta_w_true[l])) + ylim(c(-2,2))

l <- 1
qplot(1:niter, beta_wt_output_chain[,l], geom = "line") + 
  geom_hline(aes(yintercept = beta_wt_true[l])) + ylim(c(-2,2))

# alpha

p <- 3
qplot(1:niter, alpha_output_chain[,p,1], geom = "line") + 
  geom_hline(aes(yintercept = alpha1_true[p])) 

p <- 1
qplot(1:niter, alpha_output_chain[,p,2], geom = "line") + geom_hline(aes(yintercept = alpha2_true[p])) 

qplot(1:niter, phi_output_chain[,1], geom = "line") +
  geom_hline(aes(yintercept = phi_0_true)) 
qplot(1:niter, phi_output_chain[,2], geom = "line") +
  geom_hline(aes(yintercept = phi_1_true)) 

#

tracePlotParameters(beta_w_output, idx = c(2))
tracePlotParameters(beta_w_output, idx = c(2))

# 