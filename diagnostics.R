library(reshape2)

checkBeta <- function(beta_true,
                       beta_output,
                      dropFirst){
  
  ncov_z <- length(beta_true)
  
  if(dropFirst){
    beta_output <- beta_output[,,-1,drop=F]  
  }
  
  nchain <- dim(beta_output)[1]
  niter <- dim(beta_output)[2]
  
  dimnames(beta_output) <- list("chain" = NULL,"iter" = NULL,"cov" = NULL)
  
  beta_output_long <- melt(beta_output)
  
  beta_output_long$true <- beta_true[beta_output_long$cov]
  
  beta_output_long$chain <- factor(beta_output_long$chain)
  
  ggplot(data = beta_output_long, aes(x = iter,
                                      y = value,
                                      color = chain)) + 
    geom_line() + facet_grid(rows = vars(cov), scales = "free") + 
    geom_hline(aes(yintercept = true), color = "red", size = 2) + 
    xlab("Iterations") + ylab("Covariates")
  
}

checkVariances <- function(tau_output,
                           tau_true,
                           sigma_output,
                           sigma_true,
                           sigma_y_output,
                           sigma_y_true){
  
  nchain <- dim(tau_output)[1]
  niter <- dim(tau_output)[2]
  
  var_true <- c("tau" = tau_true, "sigma" = sigma_true)
  
  var_output <- array(NA, dim = c(nchain, niter, 2))
  var_output[,,1] <- tau_output
  var_output[,,2] <- sigma_output
  
  dimnames(var_output) <- list("chain" = NULL,"iter" = NULL,"var" = c("tau","sigma"))
  
  var_output_long <- melt(var_output)
  
  var_output_long$true <- var_true[var_output_long$var]
  
  var_output_long$chain <- factor(var_output_long$chain)
  
  ggplot(data = var_output_long, aes(x = iter,
                                      y = value,
                                      color = chain)) + 
    geom_line() + facet_grid(rows = vars(var), scales = "free") + 
    geom_hline(aes(yintercept = true), color = "red", size = 2) + 
    xlab("Iterations") + ylab("Covariates")
  
  
}

checkElements <- function(param_output,
                          param_true,
                          idxes){
  
  nchain <- dim(param_output)[1]
  niter <- dim(param_output)[2]
  
  param_output <- param_output[,,idxes,drop=F]
  param_true <- param_true[idxes]
  
  dimnames(param_output) <- list("chain" = NULL,"iter" = NULL, "elem" = NULL)
  
  param_output_long <- melt(param_output)
  
  param_output_long$true <- param_true[param_output_long$elem]
  
  param_output_long$chain <- factor(param_output_long$chain)
  
  ggplot(data = param_output_long, aes(x = iter,
                                      y = value,
                                      color = chain)) + 
    geom_line() + facet_grid(rows = vars(elem), scales = "free") + 
    geom_hline(aes(yintercept = true), color = "red", size = 1) + 
    xlab("Iterations") + ylab("Covariates")
  
  
}

