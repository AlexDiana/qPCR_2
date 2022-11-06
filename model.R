library(MASS); library(ars); library(ggplot2)
library(Rcpp); library(RcppArmadillo)
library(here)
sourceCpp(here("Functions/functions.cpp"))
source("~/qPCR_2/functions.R", echo=TRUE)

# prior -----

sigma_alpha <- 10
sigma_beta <- 1

a_tau <- 2; b_tau <- 1
a_sigma <- 2; b_sigma <- 1

a_sigma_y <- 2
b_sigma_y <- (.25) * (a_sigma_y - 1)

mean_phi <- c(-26, 8)
sigma_phi <- 3

alpha_0 <- c(45, -1.7)

nu <- -2
sigma_nu <- 3

p0 <- .05

lambda <- 20
sigma_lambda <- 10

# proposal values
{
  sigma_v <- .5
  proposal_ratio <- .95
  Sigma_prop_0 <- diag(c(8, 1), nrow = 2)
}

# DATA CLEANING --------

numCycles <- 45

load("~/qPCR_2/Data/sara.rda")

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

w_star <- data_standard$Quantity
C_star <- data_standard$C.
P_star <- data_standard$Plate
# Ct ordered by sample!!
numSampleV <- as.numeric(table(idx_sample_K))

ncov_b <- ncol(X_b)
ncov_w <- ncol(X_w)

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
  v_output <- array(NA, dim = c(nchain, niter, sum(M)) )
  v_tilde_output <- array(NA, dim = c(nchain, niter, sum(M)) )
  delta_output <- array(NA, dim = c(nchain, niter, sum(M)) )
  gamma_output <- array(NA, dim = c(nchain, niter, sum(M)) )
  beta_w_output <- array(NA, c(nchain, niter, ncov_w + 1))
  sigma_output <- matrix(NA, nchain, niter)
  alpha_output <- array(NA, dim = c(nchain, niter, n_P, 2))
  phi_output <- array(NA, dim = c(nchain, niter, 2))
  sigma_y_output <- array(NA, dim = c(nchain,  niter, n_P))
}

# params to update
{
  updBetaB <- T
  updTau <- T
  updL <- T
  updV <- T
  updVTilde <- T
  updBetaW <- T
  updSigma <- T
  updDelta <- T
  updPhi <- T
  updSigmaY <- T
  updAlpha <- T
}

for (chain in 1:nchain) {
  
  # chain output
  {
    l_output_chain <- matrix(NA, niter, n) 
    beta_output_chain <- matrix(NA, niter, ncov_b + 1)
    tau_output_chain <- matrix(NA, niter)
    v_output_chain <- matrix(NA, niter, sum(M)) 
    v_tilde_output_chain <- matrix(NA, niter, sum(M)) 
    delta_output <- array(NA, dim = c(niter, sum(M)) )
    gamma_output <- array(NA, dim = c(niter, sum(M)) )
    beta_w_output_chain <- matrix(NA, niter, ncov_w + 1)
    sigma_output_chain <- rep(NA, niter)
    alpha_output_chain <- array(NA, dim = c(niter, n_P, 2))
    phi_output_chain <- array(NA, dim = c(niter, 2))
    sigma_y_output_chain <- matrix(NA, nchain, niter)
    
    v_values <- matrix(NA, niter + nburn, sum(M)) 
    vtilde_values <- matrix(NA, niter + nburn, sum(M)) 
    phi_values <- matrix(NA, niter + nburn, 2) 
  }
  
  # starting values
  {
    
    if(trueStartingValues){
      
      l <- l_true
      beta_b <- beta_b_true
      tau <- tau_true
      
      beta_w0 <- beta_w0_true
      beta_w <- beta_w_true
      
      v <- v_true
      v_tilde <- vtilde_true
      sigma <- sigma_true
      
      phi_0 <- phi_0_true
      phi_1 <- phi_1_true
      
      delta <- as.numeric(delta_true == 1)
      gamma <- as.numeric(delta_true == 2)
      delta_star <- as.numeric(delta_star_true == 1)
      gamma_star <- as.numeric(delta_star_true == 2)
      
      sigma_y <- sigma_y_true
      alpha <- cbind(alpha1_true, alpha2_true)
     
    } else {
      
      l <- rep(0, n)
      beta_b <- rep(0, ncov_b)
      tau <- b_tau / (a_tau - 1)
      
      beta_w0 <- 0
      beta_w <- rep(0, ncov_w)
      
      v <- rep(1, sum(M))
      v_tilde <- rep(0, sum(M))
      sigma <- b_sigma / (a_sigma - 1)
      
      phi_0 <- mean_phi[1]
      phi_1 <- mean_phi[2]
      
      delta <- Ct > 0
      gamma <- rep(0, length(delta))
      delta_star <- C_star > 0
      gamma_star <- rep(0, length(delta_star))
      
      sigma_y <- rep(b_sigma_y / (a_sigma_y - 1), n_P)
      alpha <- matrix(alpha_0, n_P, 2, byrow = T)
      alpha1 <- alpha[,1]
      alpha2 <- alpha[,2]
    }
    
    beta_w0 <- 0
    
    w_hat <- (v > 0) * (exp(v) - 1)
    w_tilde <- (v_tilde > 0) * (exp(v_tilde) - 1)
    
    w <- w_hat + w_tilde
    
    w_pcr <- w[idx_sample_K]
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
    
    sigmas_v <- rep(sigma_v, sum(M))
    sigmas_vtilde <- rep(sigma_v, sum(M))
    
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
      list_beta <- updateBetaB(l, X_b, tau, sigma_beta)
      beta_b0 <- list_beta$beta_b0
      beta_b <- list_beta$beta_b
    }
    
    # UPDATE TAU -----------
    
    if(updTau){
      tau <- updateTau(l, X_b, beta_b0, beta_b, a_tau, b_tau)
    }
    
    # UPDATE L -------
    
    if(updL){
      l <- updateL(v, X_b, beta_b0, beta_b,
                   beta_w0, X_w, beta_w, tau)
    }
    
    # UPDATE V ----------
    
    if(updV){
      v <- updateV_cpp(v, Ct, l, M, beta_w0, X_w, beta_w,
                       alpha1, alpha2, P, w_tilde,
                       numSampleV, sigma, sigma_y, p0,
                       phi_0, phi_1, sigmas_v,
                       lambda, sigma_lambda)
      
      w_hat <- (v > 0) * (exp(v) - 1)
      
      w <- w_hat + w_tilde
      
      w_pcr <- w[idx_sample_K]
      
      v_values[iter,] <- v
      
      if(iter %% 100 == 0){
        
        hat_sigmas_v <- apply(v_values[1:(iter - 1),], 2, var)
        
        sigmas_v <- sqrt(proposal_ratio * hat_sigmas_v + (1 - proposal_ratio) * sigma_v^2)
        
      }
      
    }
    
    # v <- updateV(v, Ct, l, beta_w0, X_w, beta_w,
    #              alpha1, alpha2, P, w_tilde,
    #              phi_0, phi_1, lambda, sigma_lambda)
    
    # UPDATE V TILDE ----
    
    if(updVTilde){
      v_tilde <- updateVtilde_cpp(v_tilde, v, Ct, l, M,
                                  beta_w0, X_w, beta_w, alpha1,
                                  alpha2, P, nu, sigma_nu,
                                  numSampleV, sigma, sigma_y, p0,
                                  phi_0, phi_1, sigmas_vtilde, lambda, sigma_lambda)
      
      w_tilde <- (v_tilde > 0) * (exp(v_tilde) - 1)
      
      w <- w_hat + w_tilde
      
      w_pcr <- w[idx_sample_K]
      
      vtilde_values[iter,] <- v_tilde
      
      if(iter %% 100 == 0){
        
        hat_sigmas_v <- apply(vtilde_values[1:(iter - 1),], 2, var)
        
        sigmas_vtilde <- sqrt(proposal_ratio * hat_sigmas_v +
                                (1 - proposal_ratio) * sigma_v^2)
        
      }
      
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
                                             alpha2, sigma_y, P, P_star, lambda,
                                             sigma_lambda)
      delta <- list_deltagamma$delta
      gamma <- list_deltagamma$gamma
      delta_star <- list_deltagamma$delta_star
      gamma_star <- list_deltagamma$gamma_star
      
      if(iter > 1 & any(C_star > 0 & delta_star != 1)){
        browser()
      }
    }
    
    # UPDATE SIGMA -----------
    
    if(updSigma){
      sigma <- updateSigma(l, v, X_w, beta_w, beta_w0,
                           a_sigma, b_sigma, idx_site)
      
    }
    
    # UPDATE BETA W -------
    
    if(updBetaW){
      list_betaw <- updateBetaW(l, v, X_w, sigma, sigma_beta, 
                                idx_site, isBetaW0 = T)
      beta_w0 <- list_betaw$beta_w0
      beta_w <- list_betaw$beta_w
    }
    
    # UPDATE PHI --------
    
    if(updPhi){
      list_phi <- updatePhi(Ct, C_star, w, w_star,
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
        v_output_chain[trueIter,] <- v
        v_tilde_output_chain[trueIter,] <- v_tilde
        sigma_output_chain[trueIter] <- sigma
        beta_w_output_chain[trueIter,] <- c(beta_w0, beta_w)
        alpha_output_chain[trueIter,,] <- alpha
        phi_output_chain[trueIter,] <- c(phi_0, phi_1)
        sigma_y_output_chain[trueIter,] <- sigma_y
        
      }
      
    }
    
  }
  
  {
    l_output[chain,,] <- l_output_chain
    beta_output[chain,,] <- beta_output_chain
    tau_output[chain,] <- tau_output_chain
    v_output[chain,,] <- v_output_chain
    v_tilde_output[chain,,] <- v_tilde_output_chain
    beta_w_output[chain,,] <- beta_w_output_chain
    sigma_output[chain,] <- sigma_output_chain
    alpha_output[chain,,,] <- alpha_output_chain
    phi_output[chain,,] <- phi_output_chain
    sigma_y_output[chain,,] <- sigma_y_output_chain
  }
  
}

# OUTPUT TRUE ------

library(coda)

diagnosticsCheck <- function(params_output){
  
  dims <- dim(params_output)
  
  params_output <- apply(params_output, setdiff(1:length(dims), 1:2), c)
  
  dims <- dim(params_output)
  dims_margin <- setdiff(1:length(dims), 1)
  apply(params_output, dims_margin, function(x){
    x_current <- mcmc(x)
    effectiveSize(x_current)
  })
  
}

ess_l <- diagnosticsCheck(l_output)
ess_beta <- diagnosticsCheck(beta_output)
ess_beta_w <- diagnosticsCheck(beta_w_output)
ess_alpha <- diagnosticsCheck(alpha_output)
ess_v <- diagnosticsCheck(v_output)
ess_vtilde <- diagnosticsCheck(v_tilde_output)
ess_tau <-  diagnosticsCheck(tau_output)
ess_phi <- diagnosticsCheck(phi_output)

qplot(1:niter, alpha_output[,,19,1])

covariatePlot <- function(beta_output, covNames){
  
  beta_CI <- apply(beta_output[,,-1,drop = F], 3, function(x){
    quantile(x, probs = c(.05,.5,.95))
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
      # ggplot2::scale_x_discrete(name = "Species", 
                                # breaks = idxSpecies,
                                # labels = namesSpecies) +
      ggplot2::scale_y_continuous(breaks = (-10):10) + ggplot2::coord_flip())
  
  
  
}

covName <- colnames(X_b)

covariatePlot(beta_output, colnames(X_b))
covariatePlot(beta_w_output, colnames(X_w))

# diagnostics plot

qplot(1:niter, l_output[,,1])

qplot(1:niter, beta_output_chain[,1])
qplot(1:niter, beta_w_output[,,3])
qplot(1:niter, beta_output_chain[,3])

# OUTPUT ------

qplot(1:niter, v_output_chain[,1])  + geom_hline(aes(yintercept = v_true[1]))

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

# alpha

p <- 3
qplot(1:niter, alpha_output_chain[,p,1], geom = "line") + 
  geom_hline(aes(yintercept = alpha1_true[p])) 

p <- 1
qplot(1:niter, alpha_output_chain[,p,2], geom = "line") + geom_hline(aes(yintercept = alpha2_true[p])) 

qplot(1:niter, phi_output_chain[,1], geom = "line") + geom_hline(aes(yintercept = phi_0_true)) 
qplot(1:niter, phi_output_chain[,2], geom = "line") + geom_hline(aes(yintercept = phi_1_true)) 


