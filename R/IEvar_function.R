#' @export
IEvar <- function(ygroup, ps = c('true', 'estimated'), scores = NULL) {
  
  ps <- match.arg(ps)
  n_neigh <- dim(ygroup)[1]
  alpha <- as.numeric(dimnames(ygroup)[[2]])
  
  ie_var <- cov(ygroup)
  ie_var <- ie_var * (n_neigh - 1) / (n_neigh ^ 2)
  
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
  ypop <- apply(ygroup, 2, mean)
  
  C21 <- array(0, dim = c(length(alpha), num_gamma))
  D12 <- array(0, dim = c(num_gamma, length(alpha)))
  
  for (nn in 1 : n_neigh) {
    C21 <- C21 - t(ygroup[nn, , drop = FALSE]) %*% t(scores[, nn, drop = FALSE])
    D12 <- D12 + scores[, nn, drop = FALSE] %*% (ygroup[nn, , drop = FALSE] - ypop)
  }
  C21 <- C21 / n_neigh
  D12 <- D12 / n_neigh
  
  chol_B11_inv <- chol(B11_inv)
  mat1 <- C21 %*% t(chol_B11_inv) %*% t(C21 %*% t(chol_B11_inv))
  mat <- C21 %*% B11_inv %*% D12
  
  var_est_ps <- mat1 + mat + t(mat)
  var_est_ps <- ie_var + var_est_ps / n_neigh
  
  return(var_est_ps) 
}
