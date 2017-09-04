#' Variance of the population average potential outcome for a correctly
#' specified propensity score model.
#' 
#' @param 
VarEstPS <- function(dta, ygroup, ypop, neigh_ind, phi_hat, cov_cols, var_true,
                     trt_name = NULL) {
 
  n_neigh <- length(neigh_ind)
  num_gamma <- length(phi_hat$coefs) + 1
  phi_hat$coefs <- matrix(phi_hat$coefs, ncol = 1)
  alpha <- as.numeric(dimnames(ygroup)[[3]])
  
  # Calculating the score of the ps for every cluster.
  scores <- CalcScore(dta, neigh_ind, phi_hat, cov_cols, trt_name)

  # --- Calculating B11, the information matrix of the cluster ps.
  B11 <- matrix(0, nrow = num_gamma, ncol = num_gamma)
  for (nn in 1 : n_neigh) {
    scores_nn <- scores[, nn, drop = FALSE]
    B11 <- B11 + scores_nn %*% t(scores_nn)
  }
  B11 <- B11 / n_neigh
  
  # ---- Calculating A21, B12.
  A21 <- array(0, dim = c(2, num_gamma, length(alpha)))
  B12 <- array(0, dim = c(num_gamma, 2, length(alpha)))
  for (nn in 1 : n_neigh) {
    scores_nn <- scores[, nn]
    for (aa in 1 : length(alpha)) {
      for (it in c(1, 2)) {
        A21[it, , aa] <- A21[it, , aa] - ygroup[nn, it, aa] * scores_nn
        B12[, it, aa] <- B12[it, , aa] + scores_nn * (ygroup[nn, it, aa] - 
                                                        ypop[it, aa])
      }
    }
  }
  A21 <- A21 / n_neigh
  B12 <- B12 / n_neigh
  
  var_est_ps <- array(NA, dim = c(2, 2, length(alpha)))
  dimnames(var_est_ps) <- list(it = c(0, 1), it = c(0, 1), alpha = alpha)

  B11_inv <- chol2inv(chol(B11))
  for (aa in 1 : length(alpha)) {
    var_est_ps[, , aa] <- 
      A21[, , aa] %*% B11_inv %*% (t(A21[, , aa]) + B12[, , aa]) +
      t(B12[, , aa]) %*% B11_inv %*% t(A21[, , aa])
    var_est_ps[, , aa] <- var_est_ps[, , aa] / n_neigh
    var_est_ps[, , aa] <- var_true[, , aa] + var_est_ps[, , aa]
  }
  
  return(var_est_ps)
}