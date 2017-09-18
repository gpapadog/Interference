#' @export
IEvar <- function(ygroup, alpha, ps = c('true', 'estimated'), scores = NULL) {
  
  ps <- match.arg(ps)
  ie_var <- array(NA, dim = c(2, 2, length(alpha), length(alpha)))
  for (a1 in 1 : (length(alpha) - 1)) {
    for (a2 in (a1 + 1) : length(alpha)) {
      ie_var[, , a1, a2] <- cov(cbind(ygroup[, a1], ygroup[, a2]))
      ie_var[, , a2, a1] <- ie_var[, , a1, a2]
    }
  }
  if (ps == 'true') {
    return(ie_var)
  }
  
  # Else, we are based on the estimated propensity score.
  if (is.null(scores)) {
    stop('Provide score matrix.')
  }
  
  var_est_ps <- array(0, dim = dim(ie_var))
  num_gamma <- dim(scores)[1]
  
  # --- Calculating B11, the information matrix of the cluster ps.
  B11 <- matrix(0, nrow = num_gamma, ncol = num_gamma)
  for (nn in 1 : n_neigh) {
    scores_nn <- scores[, nn, drop = FALSE]
    B11 <- B11 + scores_nn %*% t(scores_nn)
  }
  B11 <- B11 / n_neigh
  B11_inv <- chol2inv(chol(B11))
  
  # Calculating C21, and D12.
  ypop <- apply(ygroup, 1, mean)
  
  for (a1 in 1 : (length(alpha) - 1)) {
    for (a2 in (a1 + 1) : length(alpha)) {
      
      C21 <- array(0, dim = c(2, num_gamma))
      D12 <- array(0, dim = c(num_gamma, 2))
      
      for (nn in 1 : n_neigh) {
        C21[1, ] <- C21[1, ] - ygroup[nn, a1] * scores[, nn]
        C21[2, ] <- C21[2, ] - ygroup[nn, a2] * scores[, nn]
        
        D12[, 1] <- D12[, 1] + scores[, nn] * (ygroup[nn, a1] - ypop[a1])
        D12[, 2] <- D12[, 2] + scores[, nn] * (ygroup[nn, a2] - ypop[a2])
      }
      C21 <- C21 / n_neigh
      D12 <- D12 / n_neigh
      
      chol_B11_inv <- chol(B11_inv)
      mat1 <- C21 %*% t(chol_B11_inv) %*% t(C21 %*% t(chol_B11_inv))
      mat <- C21 %*% B11_inv %*% D12
      
      var_est_ps[, , a1, a2] <- mat1 + mat + t(mat)
      var_est_ps[, , a1, a2] <- var_est_ps[, , a1, a2] / n_neigh
      var_est_ps[, , a1, a2] <- ie_var[, , a1, a2] + var_est_ps[, , a1, a2]
      var_est_ps[, , a2, a1] <- var_est_ps[, , a1, a2]
    }
  }
  return(var_est_ps) 
}
