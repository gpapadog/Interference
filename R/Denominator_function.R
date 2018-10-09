#' Calculating the denominator of the group estimator.
#' 
#' @param A The treatment vector for units in the group.
#' @param X The covariates for units in the group.
#' @param phi_hat A list with two elements. The first one is a vector of
#' coefficients of the ps, and the second one is the random effect variance.
#' If the second element is 0, the propensity score excludes random effects.
#' @param alpha This value of alpha is used to stabilize calculations. If
#' include_alpha is set to TRUE, alpha needs to coincide with the alpha used
#' for calculating the numerator.
#' @param integral_bound If the propensity score includes a random effect, the
#' integral in the denominator is calculated over the normal distribution of
#' the random effects from - integral_bound to integral_bound.
#' @param include_alpha If include_alpha is set to true, the probabilities in
#' the denominator are divided by the specified value of alpha to stabilize the
#' integral calculation.
#' 
#' @export
Denominator <- function(A, X, phi_hat, alpha = NULL, integral_bound = 10,
                        include_alpha = TRUE) {
  
  integral_bound <- abs(integral_bound)
  if (is.null(alpha)) {
    if (include_alpha) {
      stop('No alpha provided.')
    } # else {  # I dont think I need this part of the code.
    #   alpha <- mean(A)
    # }
  }
  
  X <- as.matrix(cbind(1, X))
  re_sd <- sqrt(phi_hat[[2]])
  phi_hat[[1]] <- matrix(phi_hat[[1]], nrow = length(phi_hat[[1]]), ncol = 1)
  
  # Creating the function that we will integrate over.
  f_int <- function(b) {
    r <- 1
    lin_pred <- X %*% phi_hat[[1]]  # Includes intercept.
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
    if (re_sd > 0) { # If re_sd = 0, there is no random effect in the ps model.
      r <- r * dnorm(b, mean = 0, sd = re_sd)
    }
    return(r)
  }
  
  if (re_sd > 0) {
    ans <- integrate(f_int, lower = - integral_bound * re_sd,
                     upper = integral_bound * re_sd)
  } else {
    ans <- f_int(0)
  }

  return(ans)
}

