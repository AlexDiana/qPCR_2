# simulate data

simulateData <- function()

p_imk_star <- logistic(phi_0 + phi_1 * log(data_standard$Qty))
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

# 


niter <- 1000



for (iter in 1:niter) {
  
  
  
}