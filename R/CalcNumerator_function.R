#' Calculate the numerator in the estimator Yhat.
#' 
#' @param Ai_j A vector of length n_i - 1 including the treatment assignment of all but
#' the unit we are currently considering.
#' @param Xi_j The matrix of covariates for all units in the cluster but the one. This
#' can be the observed values or a counterfactual set of covariates.
#' @param coef_hat A vector of the ps coefficients starting with the intercept.
#' @param alpha The average probability of treatment among the n_i - 1 units.
#' @param re_alpha The fixed effect bi that gives average propensity of treatment in
#' the group equal to alpha.
CalcNumerator <- function(Ai_j, Xi_j, coef_hat, alpha, re_alpha) {
  
  coef_hat <- matrix(coef_hat, nrow = length(coef_hat), ncol = 1)
  
  lin_pred <- cbind(1, as.matrix(Xi_j)) %*% coef_hat
  lin_pred <- lin_pred + re_alpha
  probs <- expit(lin_pred)
  
  r <- (probs / alpha) ^ Ai_j * ((1 - probs) / (1 - alpha)) ^ (1 - Ai_j)
  return(list(prob = prod(r), re_alpha = re_alpha))
}