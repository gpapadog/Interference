CalcB11matrix <- function(dta, neigh_ind, cov_cols, phi_hat) {
  
  num_gamma <- length(phi_hat$coefs) + 1
  B11 <- matrix(0, nrow = num_gamma, ncol = num_gamma)
  
  n_neigh <- length(neigh_ind)
  
  for (nn in 1:n_neigh) {
    
    A <- dta$A[neigh_ind[[nn]]]
    X <- dta[neigh_ind[[nn]], cov_cols]
    
    hess_function <- function(gamma) {
      
      phi_hat_hess <- NULL
      phi_hat_hess$coefs <- gamma[- num_gamma]
      phi_hat_hess$re_var <- gamma[num_gamma]
      
      likelihood <- GroupLikelihood(A = A, X = X, phi_hat = phi_hat_hess)$value
      
      return(likelihood)
    }
    
    deriv <- grad(hess_function, x = c(phi_hat$coefs, phi_hat$re_var))
    deriv <- matrix(deriv, ncol = 1)
    B11 <- B11 + deriv %*% t(deriv)
  }
  B11 <- B11 / n_neigh
  
  return(B11)
}
