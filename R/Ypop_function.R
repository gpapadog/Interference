#' Estimates and variance of the population average potential outcome for known
#' or correctly specified propensity score model.
#' 
#' @param ygroup An array including the group average potential outcome
#' estimates where the dimensions correspond to group, individual treatment and
#' value of alpha.
#' @param ps String. Can take values 'true', or 'estimated' for known or
#' estimated propensity score. Defaults to 'true'.
#' @param scores A matrix with rows corresponding to the parameters of the
#' propensity score model and columns for groups. Includes the score of the
#' propensity score evaluated for the variables of each group.
#' @param dta The data set including the variable neigh. Defaults to NULL. Can
#' be left NULL when the true propensity score is used.
#' 
#' @export
Ypop <- function(ygroup, ps = c('true', 'estimated'), scores = NULL,
                 dta = NULL, use = 'everything') {
  
  ps <- match.arg(ps)
  n_neigh <- dim(ygroup)[1]
  alpha <- as.numeric(dimnames(ygroup)[[3]])
  ypop <- apply(ygroup, c(2, 3), mean, na.rm = TRUE)
  
  ypop_var <- apply(ygroup, 3, function(x) cov(x, use = use))
  ypop_var <- array(ypop_var, dim = c(2, 2, length(alpha)))
  # In order to get 1 / N, instead of 1 / (N - 1) in the variance estimates.
  ypop_var <- ypop_var * (n_neigh - 1) / n_neigh
  ypop_var <- ypop_var / n_neigh  # Since we have n_neigh clusters.
  dimnames(ypop_var) <- list(a = c(0, 1), a = c(0, 1), alpha = alpha)
  
  if (ps == 'true') {
    return(list(ypop = ypop, ypop_var = ypop_var))
  }
  

  # Else, the propensity score is estimated.

  neigh_ind <- sapply(1 : n_neigh, function(x) which(dta$neigh == x))

  if (is.null(scores)) {
    stop('scores needs to be specified for estimated propensity score.')
  }
  
  num_gamma <- dim(scores)[1]
  


  var_est_ps <- array(NA, dim = c(2, 2, length(alpha)))
  dimnames(var_est_ps) <- list(it = c(0, 1), it = c(0, 1), alpha = alpha)
  
  # --- Calculating B11, the information matrix of the cluster ps.
  B11 <- matrix(0, nrow = num_gamma, ncol = num_gamma)
  for (nn in 1 : n_neigh) {
    scores_nn <- scores[, nn, drop = FALSE]
    B11 <- B11 + scores_nn %*% t(scores_nn)
  }
  B11 <- B11 / n_neigh
  B11_inv <- chol2inv(chol(B11))
  
  
  # ---- Calculating A21, B12 for each alpha.
  for (aa in 1 : length(alpha)) {
    
    A21 <- array(0, dim = c(2, num_gamma))
    B12 <- array(0, dim = c(num_gamma, 2))
    
    for (it in c(1, 2)) {
      for (nn in 1 : n_neigh) {
        A21[it, ] <- A21[it, ] - ygroup[nn, it, aa] * scores[, nn]
        B12[, it] <- B12[, it] + scores[, nn] * (ygroup[nn, it, aa] -
                                                   ypop[it, aa])
      }
    }
    A21 <- A21 / n_neigh
    B12 <- B12 / n_neigh
    
    chol_B11_inv <- chol(B11_inv)
    mat1 <- A21 %*% t(chol_B11_inv) %*% t(A21 %*% t(chol_B11_inv))
    mat <- A21 %*% B11_inv %*% B12
    
    var_est_ps[, , aa] <- mat1 + mat + t(mat)
    var_est_ps[, , aa] <- var_est_ps[, , aa] / n_neigh
    var_est_ps[, , aa] <- ypop_var[, , aa] + var_est_ps[, , aa]
  }
  
  return(list(ypop = ypop, ypop_var = var_est_ps))
  
}
