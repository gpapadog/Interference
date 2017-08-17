#' Function that calculates the integral in the denominator corresponding to the
#' marginal likelihood of the vector A, given X and the fixed effects without a
#' cluster level random effect.
Denom_woRE <- function(A, X, coef_hat, alpha) {
  
  X <- as.matrix(cbind(1, X))
  coef_hat <- matrix(coef_hat, nrow = length(coef_hat), ncol = 1)
  
  lin_pred <- X %*% coef_hat
  r <- 1
  for (ii in 1:length(A)) {
    prob_trt <- expit(lin_pred[ii])
    r <- r * (prob_trt / alpha) ^ A[ii]
    r <- r * ((1 - prob_trt) / (1 - alpha)) ^ (1 - A[ii])
  }
  return(list(value = r))
}

