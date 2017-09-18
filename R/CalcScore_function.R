#' @param dta The data set including (at least) the treatment and covaratiates.
#' @param neigh_ind A list including one element for each neighborhood. That
#' element includes the indices of the observations in dta that belong to each
#' neighborhood.
#' @param phi_hat A list with two elements. The first one is a vector of
#' coefficients of the propensity score, and the second one is the random
#' effect variance.
#' @param cov_cols The indeces including the covariates of the propensity score
#' model.
#' @param trt_name The name of the treatment column. If it is 'A', you can
#' leave NULL.
#' 
#' @export
CalcScore <- function(dta, neigh_ind, phi_hat, cov_cols, trt_name = NULL,
                      integral_bound = 10) {
  
  dta <- as.data.frame(dta)
  n_neigh <- length(neigh_ind)
  num_gamma <- length(phi_hat$coefs) + 1
  phi_hat$coefs <- matrix(phi_hat$coefs, ncol = 1)

  if (is.null(trt_name)) {
    trt_name <- 'A'
  }
  trt_col <- which(names(dta) == 'A')

  scores <- matrix(NA, nrow = num_gamma, ncol = n_neigh)
  
  for (nn in 1 : n_neigh) {
    
    Ai <- dta[neigh_ind[[nn]], trt_col]
    Xi <- dta[neigh_ind[[nn]], cov_cols]
    
    hess_function <- function(gamma) {
      
      phi_hat_hess <- NULL
      phi_hat_hess$coefs <- gamma[- num_gamma]
      phi_hat_hess$re_var <- gamma[num_gamma]
      
      likelihood <- DenomIntegral(Ai, Xi, phi_hat_hess, include_alpha = FALSE,
                                  integral_bound = integral_bound)
      return(log(likelihood$value))
    }
    
    scores[, nn] <- numDeriv::grad(hess_function, x = c(phi_hat$coefs, phi_hat$re_var))

  }
  return(scores)  
}
