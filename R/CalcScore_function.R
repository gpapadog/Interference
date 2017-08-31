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
CalcScore <- function(dta, neigh_ind, phi_hat, cov_cols, trt_name = NULL) {
  
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

    for (jj in 1:(num_gamma - 1)) {  # For each PS coefficient.
      
      fn_deriv <- function(gamma) {
        
        phi_hat_jj <- phi_hat
        phi_hat_jj$coefs[jj, 1] <- gamma
        
        denom <- DenomIntegral(A = Ai, X = Xi, phi_hat = phi_hat_jj,
                               alpha = curr_alpha)$value
        return(log(denom))
      }
      scores[nn, jj] <- grad(fn_deriv, phi_hat$coefs[jj])
    }

    fn_deriv <- function(gamma) {
      phi_hat_var <- phi_hat
      phi_hat_var[[2]] <- gamma
      denom <- DenomIntegral(A = Ai, X = Xi, phi_hat = phi_hat_var,
                             alpha = curr_alpha)$value
      return(log(denom))
    }
    scores[nn, num_gamma] <- grad(fn_deriv, phi_hat[[2]])

  }
  return(scores)  
}
