#' Function that calculates the integral in the denominator corresponding to the
#' marginal likelihood of the vector A, given X and the fixed effects.
DenomIntegral <- function(A, X, phi_hat, alpha, integral_bound = 10) {
  
  integral_bound <- abs(integral_bound)
  
  X <- as.matrix(cbind(1, X))
  re_sd <- sqrt(phi_hat[[2]])
  phi_hat[[1]] <- matrix(phi_hat[[1]], nrow = length(phi_hat[[1]]), ncol = 1)
  # Creating the function that we will integrate over.
  f_int <- function(b) {
    r <- 1
    lin_pred <- X %*% phi_hat[[1]]
    for (ii in 1:length(A)) {
      prob_trt <- expit(lin_pred[ii] + b)
      r <- r * (prob_trt / alpha) ^ A[ii]
      r <- r * ((1 - prob_trt) / (1 - alpha)) ^ (1 - A[ii])
    }
    r <- r * dnorm(b, mean = 0, sd = re_sd)
    return(r)
  }
  ans <- integrate(f_int, lower = - integral_bound * re_sd,
                   upper = integral_bound * re_sd)
  return(ans)
}

