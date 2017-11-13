#' Function that calculates the integral in the denominator corresponding to the
#' marginal likelihood of the vector A, given X and the fixed effects.
#' 
#' @export
DenomIntegral <- function(A, X, phi_hat, alpha = NULL, integral_bound = 10,
                          include_alpha = TRUE) {
  
  integral_bound <- abs(integral_bound)
  if (is.null(alpha)) {
    if (include_alpha) {
      stop('No alpha provided.')
    } else {
      alpha <- mean(A)
    }
  }
  
  X <- as.matrix(cbind(1, X))
  re_sd <- sqrt(phi_hat[[2]])
  phi_hat[[1]] <- matrix(phi_hat[[1]], nrow = length(phi_hat[[1]]), ncol = 1)
  # Creating the function that we will integrate over.
  f_int <- function(b) {
    r <- 1
    lin_pred <- X %*% phi_hat[[1]]
    for (ii in 1:length(A)) {
      prob_trt <- expit(lin_pred[ii] + b)
      if (include_alpha) {
        success_weight <- prob_trt / alpha
        failure_weight <- (1 - prob_trt) / (1 - alpha)
      } else {
        success_weight <- prob_trt
        failure_weight <- 1 - prob_trt
      }
      r <- r * success_weight ^ A[ii] * failure_weight ^ (1 - A[ii])
    }
    r <- r * dnorm(b, mean = 0, sd = re_sd)
    return(r)
  }
  ans <- integrate(f_int, lower = - integral_bound * re_sd,
                   upper = integral_bound * re_sd)
  
  return(ans)
}

