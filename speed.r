library(microbenchmark)

microbenchmark({
  list_beta <- updateBetaB(l, X_b, tau, 
                           prior_betab0, 
                           sigma_beta0, sigma_beta)
},{
  tau <- updateTau(l, X_b, beta_b0, beta_b, a_tau, b_tau)

},{
  l <- updateL(v, X_b, beta_b0, beta_b,
               beta_w0, X_w, beta_w, tau)
},{
  list_v <- updateV_cpp(v, v_tilde, Ct, l, M, beta_w0, X_w, beta_w,
                        sigma, nu, X_wt, beta_wt, sigma_nu,
                        alpha1, alpha2, P, w_tilde,
                        numSampleV, sigma_y, p0,
                        phi_0, phi_1, sigmas_v,
                        lambda, sigma_lambda)
},{
  v_tilde <- updateVtilde_cpp(v_tilde, v, Ct, n, M,
                              alpha1, alpha2, P, nu, X_wt, 
                              beta_wt, sigma_nu,
                              numSampleV, sigma_y, p0,
                              phi_0, phi_1, sigmas_vtilde, 
                              lambda, sigma_lambda)
},{
  list_vvtilde <- updateV_Vtilde_cpp(v, v_tilde, Ct, 
                                     l, M, beta_w0, X_w, 
                                     beta_w, sigma, 
                                     nu, sigma_nu, X_wt, 
                                     beta_wt, alpha1, alpha2, 
                                     P, numSampleV, 
                                     sigma_y, sigmas_vvtilde, 
                                     p0, phi_0, phi_1, 
                                     lambda, sigma_lambda)
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
                       a_sigma, b_sigma, idx_site)
},{
  list_betaw <- updateBetaW(l, v, X_w, sigma, sigma_beta, 
                            idx_site, isBetaW0 = T)
},{
  list_betawt <- updateBetaWt(v_tilde, v, X_wt, 
                              nu0, sigma_nu, sigma_beta, 
                              ncov_wt)
},{
  sigma_nu <- updateSigmaNu(v_tilde, v, X_wt, nu, 
                            a_sigma_nu, b_sigma_nu)
},{
  list_phi <- updatePhi(Ct, C_star, w, w_star,
                        delta, delta_star,
                        phi_0, phi_1,
                        mean_phi, sigma_phi,
                        Sigma_prop)
},{
  list_phi <- updatePhi_PG(Ct, C_star, w, w_star,
                           delta, delta_star,
                           phi_0, phi_1,
                           mean_phi, sigma_phi,
                           Sigma_prop)
},{
  p0 <- updatep0(delta, delta_star,
                 gamma, gamma_star,
                 a_p0, b_p0)
},{
  sigma_y <- updateSigmaY_cpp(Ct, w_pcr, delta, C_star, w_star,
                              delta_star, alpha1,
                              alpha2, sigma_alpha, P, P_star, n_P, 
                              a_sigma_y, b_sigma_y)
},{
  alpha <- updateAlpha_cpp(Ct, w_pcr, delta, P, 
                           C_star, w_star, delta_star, P_star, 
                           sigma_y, sigma_alpha, n_P, alpha_0)
},times = 10)
