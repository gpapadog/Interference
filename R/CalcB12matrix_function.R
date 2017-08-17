CalcB12matrix <- function(dta, neigh_ind, cov_cols, phi_hat, group_aver,
                          pop_aver, re_sd = 2) {
  
  alpha <- as.numeric(dimnames(group_aver)[[3]])
  if (any(alpha != dimnames(pop_aver)[[2]])) {
    stop('alpha for group_aver and pop_aver need to correspond.')
  }
  
  if (any(is.na(group_aver))) {
    warning('Some group averages are missing.')
  }
  
  num_gamma <- length(phi_hat$coefs) + 1
  n_neigh <- length(neigh_ind)
  
  B12 <- array(0, dim = c(num_gamma, 2, length(alpha)))
  dimnames(B12) <- list(ps_coef = 1:num_gamma, a = c(0, 1), alpha = alpha)
  
  for (aa in 1:length(alpha)) {
    for (it in c(0, 1)) {
      
      r <- matrix(0, nrow = n_neigh, ncol = num_gamma)
      
      for (nn in 1:n_neigh) {
        
        A <- dta$A[neigh_ind[[nn]]]
        X <- dta[neigh_ind[[nn]], cov_cols]
        
        der_function <- function(gamma) {
          
          phi_hat_hess <- NULL
          phi_hat_hess$coefs <- gamma[- num_gamma]
          phi_hat_hess$re_var <- gamma[num_gamma]
          
          likelihood <- GroupLikelihood(A = A, X = X, phi_hat = phi_hat_hess)
          
          return(likelihood)
        }
        
        deriv <- grad(der_function, x = c(phi_hat$coefs, phi_hat$re_var))
        r[nn, ] <- deriv * (group_aver[nn, it + 1, aa] - pop_aver[it + 1, aa])
      }
      B12[, it + 1, aa] <- apply(r, 2, mean, na.rm = TRUE)
    }
  }
  return(B12)
}
