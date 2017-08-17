CalcA21matrix <- function(dta, ygroup, neigh_ind, phi_hat, alpha, lower, upper) {
  
  dta <- as.data.frame(dta)
  n_neigh <- length(neigh_ind)
  
  phi_hat$coefs <- matrix(phi_hat$coefs, ncol = 1)
  num_gamma <- length(phi_hat$coefs) + 1  # Coefficients & random effect variance.
  
  A21 <- array(0, dim = c(2, num_gamma, length(alpha)))
  dimnames(A21) <- list(a = c(0, 1), gamma = 1:num_gamma, alpha = alpha)
  
  denom_deriv <- array(NA, dim = c(n_neigh, num_gamma, length(alpha)))
  dimnames(denom_deriv) <- list(neigh = 1:n_neigh, gamma = 1:num_gamma,
                                alpha = alpha)
  
  # The derivative of the denominator with respect 
  for (aa in 1:length(alpha)) {
    curr_alpha <- alpha[aa]
    for (nn in 1:n_neigh) {
      Ai <- matrix(dta$A[neigh_ind[[nn]]])
      Xi <- dta[neigh_ind[[nn]], cov_cols]
      for (gg in 1:num_gamma) {
        
        phi_hat_jj <- phi_hat
        enum_at <- ifelse(gg < num_gamma, phi_hat$coefs[gg], phi_hat$re_var)
        
        fn_deriv <- function(gamma) {

          if (gg < num_gamma) {  # coefficient.
            phi_hat_jj$coefs[gg, 1] <- gamma
          } else {  # residual variance
            phi_hat_jj$re_var <- gamma
          }
          denom <- DenomIntegral(A = Ai, X = Xi, phi_hat = phi_hat_jj,
                                 alpha = curr_alpha)$value
          return(denom)
        }
        denom_deriv[nn, gg, aa] <- grad(fn_deriv, enum_at)
      }
    }
  }
  
  denom <- sapply(1:length(alpha), function(aa) {
    sapply(1:n_neigh, function(nn) {
      DenomIntegral(A = dta$A[neigh_ind[[nn]]],
                    X = dta[neigh_ind[[nn]], cov_cols],
                    phi_hat = phi_hat, alpha = alpha[aa])$value
    })
  })
  
  if (any(is.na(ygroup))) {
    warning('ygroup includes missing values.')
  }
  for (it in c(0, 1)) {
    for (gg in 1:num_gamma) {
      for (aa in 1:length(alpha)) {
        A21[it + 1, gg, aa] <- - mean(denom_deriv[, gg, aa] / denom[, aa] *
                                        ygroup[, it + 1, aa], na.rm = TRUE)
      }
    }
  }
  
  # The numerator does not depend on the residual variance.
  
  
  
  
  
  
  for (aa in 1:length(alpha)) {
    curr_alpha <- alpha[aa]
    for (nn in 1:n_neigh) {
      
      # r is the 2 x num_gamma matrix including the entries from cluster nn.
      r <- matrix(0, nrow = 2, ncol = num_gamma)
      
      for (ii in 1:length(neigh_ind[[nn]])) {
        
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
        phi_hat_var <- phi_hat
        numer <- CalcNumerator(Ai_j = Ai_j, Xi_j = Xi_j, phi_hat = phi_hat_var,
                               alpha = curr_alpha, lower = lower, upper = upper)
        
        fn_deriv <- function(gamma) {
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
      r_denom <- length(neigh_ind[[nn]]) * curr_alpha ^ a * (1 - curr_alpha) ^ (1 - a)
      A21[, , aa] <- A21[, , aa] + r / r_denom
    }
    A21[, , aa] <- A21[, , aa] / n_neigh
  }
  
  if (any(A21 == 0)) {
    warning('At least one entry of A21 is equal to 0.')
  }
  
  return(A21)
}


