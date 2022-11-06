library(microbenchmark)

microbenchmark({
  list_beta <- updateBetaB(l, X_b, tau, sigma_beta)
  beta_b0 <- list_beta$beta_b0
  beta_b <- list_beta$beta_b
},{
  tau <- updateTau(l, X_b, beta_b0, beta_b, a_tau, b_tau)

},{
  l <- updateL(v, X_b, beta_b0, beta_b,
               beta_w0, X_w, beta_w, tau)
},{
  v <- updateV_cpp(v, Ct, l, M, beta_w0, X_w, beta_w,
                   alpha1, alpha2, P, w_tilde,
                   numSampleV, sigma, sigma_y, p0,
                   phi_0, phi_1, sigmas_v,
                   lambda, sigma_lambda)
},{
  v_tilde <- updateVtilde_cpp(v_tilde, v, Ct, l, M, beta_w0, X_w, beta_w, alpha1,
                              alpha2, P, nu, sigma_nu, numSampleV, sigma, sigma_y, p0,
                              phi_0, phi_1, sigmas_v, lambda, sigma_lambda)
},{
  list_deltagamma <- updateDeltaGammaCpp(Ct, C_star, w_pcr, w_star,
                                         phi_0, phi_1, p0, alpha1,
                                         alpha2, sigma_y, P, P_star, lambda,
                                         sigma_lambda)
  delta <- list_deltagamma$delta
  gamma <- list_deltagamma$gamma
  delta_star <- list_deltagamma$delta_star
  gamma_star <- list_deltagamma$gamma_star
},{
  sigma <- updateSigma(l, v, X_w, beta_w, beta_w0,
                       a_sigma, b_sigma)
},{
  list_betaw <- updateBetaW(l, v, X_w, sigma, sigma_beta, idx_site, F)
  beta_w0 <- list_betaw$beta_w0
  beta_w <- list_betaw$beta_w
},{
  list_phi <- updatePhi(Ct, C_star, w, w_star,
                        delta, delta_star,
                        phi_0, phi_1,
                        mean_phi, sigma_phi,
                        Sigma_prop)
  phi_0 <- list_phi$phi_0
  phi_1 <- list_phi$phi_1
},{
  sigma_y <- updateSigmaY_cpp(Ct, w_pcr, delta, C_star, w_star, delta_star, alpha1,
                             alpha2, sigma_alpha, P, P_star, a_sigma_y, b_sigma_y)
},{
  alpha <- updateAlpha_cpp(Ct, w_pcr, delta, P, 
                           C_star, w_star, delta_star, P_star, 
                           sigma_y, sigma_alpha, n_P, alpha_0)
},times = 10)
