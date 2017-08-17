#' @param dta The dataset as a data frame including treatment, outcome and covariates.
#' @param cov_cols The indeces including the covariates of the propensity score model.
#' @param phi_hat A list with two elements. The first one is a vector of coefficients
#' of the propensity score, and the second one is the random effect variance.
#' @param alpha The values of alpha for which we want to estimate the group average
#' potential outcome.
#' @param neigh_ind List. i^{th} element is a vector with the row indeces of dta that
#' are in cluster i. Can be left NULL.
GroupIPW <- function(dta, cov_cols, phi_hat, neigh_ind, alpha) {
  
  n_neigh <- length(neigh_ind)
  
  yhat <- array(NA, dim = c(n_neigh, 2, length(alpha)))
  dimnames(yhat) <- list(neigh = 1:n_neigh, trt = c(0, 1), alpha = alpha)
  
  for (aa in 1:length(alpha)) {
    curr_alpha <- alpha[aa]
    for (it in c(0, 1)) {
      curr_it <- it
      for (nn in 1:n_neigh) {
        
        # We only calculate group average if there are individuals with that treatment.
        if (any(dta$A[neigh_ind[[nn]]] == curr_it)) {
          y_curr <- 0
          
          bern_prob <- curr_alpha ^ curr_it * (1 - curr_alpha) ^ (1 - curr_it)
          
          for (ind in neigh_ind[[nn]]) {
            y_curr <- y_curr + (dta$A[ind] == curr_it) * dta$Y[ind] / bern_prob
          }
          
          denom <- DenomIntegral(A = dta$A[neigh_ind[[nn]]],
                                 X = dta[neigh_ind[[nn]], cov_cols],
                                 phi_hat = phi_hat, alpha = curr_alpha)
          denom <- length(neigh_ind[[nn]]) * denom$value
          
          yhat[nn, it + 1, aa] <- y_curr / denom 
        }  # Otherwise, group average is NA.
        
      }
    }
  }
  return(yhat)
}



