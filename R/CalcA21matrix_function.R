CalcA21matrix <- function(dta, neigh_ind, phi_hat, alpha, lower, upper) {
  
  dta <- as.data.frame(dta)
  n_neigh <- length(neigh_ind)
  
  phi_hat$coefs <- matrix(phi_hat$coefs, ncol = 1)
  num_gamma <- length(phi_hat$coefs) + 1  # Coefficients & random effect variance.
  
  A21 <- array(0, dim = c(2, num_gamma, length(alpha)))
  dimnames(A21) <- list(a = c(0, 1), gamma = 1:num_gamma, alpha = alpha)
  
  for (aa in 1:length(alpha)) {
    curr_alpha <- alpha[aa]
    
    for (nn in 1:n_neigh) {
      n_i <- length(neigh_ind[[nn]])
      
      # r is the 2 x num_gamma matrix including the entries from cluster nn.
      r <- matrix(0, nrow = 2, ncol = num_gamma)
      
      for (ii in 1:n_i) {
        
        # Each observation contributes only to the it they belong.
        Aij <- dta$A[neigh_ind[[nn]][ii]]
        
        # Information on the neighbors.
        wh_neigh <- setdiff(neigh_ind[[nn]], neigh_ind[[nn]][ii])
        Ai_j <- dta$A[wh_neigh]
        Xi_j <- dta[wh_neigh, cov_cols]
        des_mat <- cbind(1, as.matrix(Xi_j))
        
        # Information on the whole cluster, same as i_j but with individual j included.
        full_neigh <- c(neigh_ind[[nn]][ii], wh_neigh)
        Ai <- c(Aij, Ai_j)
        Xi <- as.matrix(dta[full_neigh, cov_cols])
        des_mat_full <- cbind(1, Xi)

        # Random effect based only on the neighbors.
        lin_pred <- des_mat %*% phi_hat_est$coefs
        re_alpha <- FromAlphaToRE(alpha = curr_alpha, lin_pred = lin_pred,
                                  lower = lower, upper = upper)
        
        for (jj in 1:(num_gamma - 1)) {  # For each PS coefficient.
          
          fn_deriv <- function(gamma) {
            
            phi_hat_jj <- phi_hat
            phi_hat_jj$coefs[jj, 1] <- gamma
            
            numer <- CalcNumerator(Ai_j = Ai_j, Xi_j = Xi_j, phi_hat = phi_hat_jj,
                                   alpha = curr_alpha, lower = lower, upper = upper,
                                   re_alpha = re_alpha)
            denom <- DenomIntegral(A = Ai, X = Xi, phi_hat = phi_hat_jj,
                                   alpha = curr_alpha)$value
            return(numer / denom)
          }
          
          ratio_deriv <- grad(fn_deriv, phi_hat$coefs[jj])
          
          r[Aij + 1, jj] <- r[Aij + 1, jj] + ratio_deriv * dta$Y[neigh_ind[[nn]][ii]]
          
        }
        
        # Also doing it for the variance.
        numer <- CalcNumerator(Ai_j = Ai_j, Xi_j = Xi_j, phi_hat = phi_hat,
                               alpha = curr_alpha, lower = lower, upper = upper)
        
        fn_deriv <- function(gamma) {
          phi_hat_var <- phi_hat
          phi_hat_var[[2]] <- gamma
          denom <- DenomIntegral(A = Ai, X = Xi, phi_hat = phi_hat_var,
                                 alpha = curr_alpha)
          denom <- denom$value
          return(denom)
        }
        
        ratio_deriv <- numer / grad(fn_deriv, phi_hat_est$coefs[[2]])
        r[Aij + 1, num_gamma] <- r[Aij + 1, num_gamma] +
          ratio_deriv * dta$Y[neigh_ind[[nn]][ii]]
      }
      r[1, ] <- r[1, ] / (n_i * (1 - curr_alpha))
      r[2, ] <- r[2, ] / (n_i * curr_alpha)
      A21[, , aa] <- A21[, , aa] + r
    }
    A21[, , aa] <- A21[, , aa] / n_neigh
  }
  
  if (any(A21 == 0)) {
    warning('At least one entry of A21 is equal to 0.')
  }
  
  return(A21)
}


