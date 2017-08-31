VarEstPS <- function(dta, ygroup, ypop, neigh_ind, phi_hat, cov_cols, var_true,
                     trt_name = NULL) {
 
  n_neigh <- length(neigh_ind)
  num_gamma <- length(phi_hat$coefs) + 1
  phi_hat$coefs <- matrix(phi_hat$coefs, ncol = 1)
  alpha <- as.numeric(dimnames(ygroup))[[3]]
  
  scores <- CalcScore(dta, neigh_ind, phi_hat, cov_cols, trt_name)
  
  A21 <- array(0, dim = c(2, length(num_gamma), length(alpha)))
  B11 <- array(0, dim = c(length(num_gamma), length(num_gamma),
                          length(alpha)))
  B12 <- array(0, dim = c(length(num_gamma), 2, length(alpha)))
  dimnames(A21) <- list(it = c(0, 1), gamma = 1:num_gamma, alpha = alpha)
  dimnames(B11) <- list(gamma = 1 : num_gamma, gamma = 1 : num_gamma,
                        alpha = alpha)
  dimnames(B12) <- dimnames(A21)[c(2, 1, 3)]
  
  for (aa in 1 : length(alpha)) {
    for (nn in 1 : n_neigh) {
      
      scores_nn <- scores[, nn, drop = FALSE]
      B11[, , aa] <- B11[, , a] + t(scores_nn) %*% scores_nn
      
      for (it in c(1, 2)) {
        A21[it, , aa] <- A21[it, , aa] - ygroup[nn, it, aa] * scores_nn
        B12[, it, aa] <- B12[it, , aa] + scores_nn * (ygroup[nn, it, aa] - 
                                                        ypop[it, aa])
      }
    }
  }
  B11 <- B11 / n_neigh
  A21 <- A21 / n_neigh
  B12 <- B12 / n_neigh
  
  var_est_ps <- array(NA, dim = c(2, 2, length(alpha)))
  dimnames(var_est_ps) <- list(it = c(0, 1), it = c(0, 1), alpha = alpha)

  for (aa in 1 : length(alpha)) {
    B11_inv <- solve(B11[, , aa])
    var_est_ps[, , aa] <- var_true[, , aa] + A21[, , aa] %*% B11_inv %*%
      (t(A21[, , aa]) + B12) + t(B12[, , aa]) %*% B11_inv %*% t(A21[, , aa])
  }
  
  return(var_est_ps)
}