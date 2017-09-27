#' Variance of the population average potential outcome for a correctly
#' specified propensity score model.
#' 
#' @param dta The data set including (at least) the treatment and covaratiates.
#' @param yhat_group An array including the group average potential outcome
#' estimates where the dimensions correspond to group, individual treatment and
#' value of alpha.
#' @param yhat_pop A matrix including the population average potential outcome
#' estimates where dimensions are individual treatment and alpha. If left NULL,
#' it will be calculated using the means of yhat_group.
#' @param neigh_ind A list including one element for each neighborhood. That
#' element includes the indices of the observations in dta that belong to each
#' neighborhood. Can be left NULL if dta includes a column neigh.
#' @param phi_hat A list with two elements. The first one is a vector of
#' coefficients of the propensity score, and the second one is the random
#' effect variance.
#' @param cov_cols The indeces including the covariates of the propensity score
#' model.
#' @param trt_name The name of the treatment column. If it is 'A', you can
#' leave NULL.
#' 
#' @export
VarEstPS <- function(dta, yhat_group, yhat_pop = NULL, neigh_ind = NULL,
                     phi_hat, cov_cols, var_true = NULL, trt_name = NULL,
                     scores = NULL) {
 
  num_gamma <- length(phi_hat[[1]]) + 1
  
  if (is.null(neigh_ind)) {
    n_neigh <- max(dta$neigh)
    neigh_ind <- sapply(1 : n_neigh, function(x) which(dta$neigh == x))
  }
  
  n_neigh <- length(neigh_ind)
  
  phi_hat[[1]] <- matrix(phi_hat[[1]], ncol = 1)
  alpha <- as.numeric(dimnames(yhat_group)[[3]])
  
  # Calculating the score of the ps for every cluster.
  if (is.null(scores)) {
    scores <- CalcScore(dta, neigh_ind, phi_hat, cov_cols, trt_name)
  }
  
  # Getting the variance if the PS was the true.
  if (is.null(var_true)) {
    var_true <- YpopTruePS(ygroup = yhat_group, alpha = alpha)$ypop_var
  }
  
  if (is.null(yhat_pop)) {
    yhat_pop <- apply(yhat_group, c(2, 3), mean)
  }

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
        A21[it, ] <- A21[it, ] - yhat_group[nn, it, aa] * scores[, nn]
        B12[, it] <- B12[, it] + scores[, nn] * (yhat_group[nn, it, aa] -
                                                   yhat_pop[it, aa])
      }
    }
    A21 <- A21 / n_neigh
    B12 <- B12 / n_neigh
    
    chol_B11_inv <- chol(B11_inv)
    mat1 <- A21 %*% t(chol_B11_inv) %*% t(A21 %*% t(chol_B11_inv))
    mat <- A21 %*% B11_inv %*% B12
    
    var_est_ps[, , aa] <- mat1 + mat + t(mat)
    var_est_ps[, , aa] <- var_est_ps[, , aa] / n_neigh
    var_est_ps[, , aa] <- var_true[, , aa] + var_est_ps[, , aa]
  }
  
  return(var_est_ps)
}
