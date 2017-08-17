#' Function that calculates the random effect that corresponds to a specific average
#' probability of treatment.
#' 
#' @param alpha Target average probability of treatment. Numberic between 0 and 1.
#' @param lin_pred Linear predictor for the probability of treatment among the
#' remaining units. Includes intercept and fixed effects.
#' @param lower The lower end of the values for bi we will look at. Defaults to - 10.
#' @param upper The upper end of the values for bi we will look at. Defaults to 10.
FromAlphaToRE <- function(alpha, lin_pred, lower = - 10, upper = 10) {
  r <- optimise(f = AlphaToBi, lower = lower, upper = upper, alpha = alpha,
                lin_pred = lin_pred)$minimum
  if (abs(lower - r) < 0.1 | abs(upper - r) < 0.1) {
    warning(paste('bi = ', r, ', lower = ', lower, ', upper = ', upper))
  }
  return(r)
}

AlphaToBi <- function(b, alpha, lin_pred) {
  exp_lin_pred <- exp(lin_pred)
  r <- abs(mean(exp_lin_pred / (exp_lin_pred + exp(- b))) - alpha)
  return(r)
}


# QUESTION: Does it make sense to consider the average probability of treatment where
# the treatment of the unit we are considered is fixed to whatever we are considering?