YpopEstPS <- function(dataset, ygroup, neigh_ind, cov_cols, phi_hat, lower, upper,
                      alpha, use = 'everything') {
  
  n_neigh <- dim(ygroup)[1]
  # Population average potential outcome.
  ypop <- apply(ygroup, c(2, 3), mean, na.rm = TRUE)
  
  # Asymptotic variance.
  
  # B22 has the same form as the asymptotic variance with known ps.
  B22 <- apply(ygroup, 3, function(x) cov(x, use = use))
  B22 <- array(B22, dim = c(2, 2, length(alpha))) * (n_neigh - 1) / n_neigh
  
  A21 <- CalcA21matrix(dta = dataset, neigh_ind = neigh_ind,
                       phi_hat = phi_hat, alpha = alpha, lower = lower,
                       upper = upper)
  
  B12 <- CalcB12matrix(dataset, neigh_ind, cov_cols, phi_hat, ygroup, ypop)
  
  B11 <- CalcB11matrix(dataset, neigh_ind, cov_cols, phi_hat = phi_hat)
  B11_inv <- solve(B11)
  
  ypop_var <- array(NA, dim = c(2, 2, length(alpha)))

  for (aa in 1:length(alpha)) {
    ypop_var[, , aa] <- A21[, , aa] %*% B11_inv %*% (t(A21[, , aa]) + B12[, , aa]) +
      t(B12[, , aa]) %*% B11_inv %*% t(A21[, , aa]) + B22[, , aa]
  }
  ypop_var <- ypop_var / n_neigh

  return(list(ypop = ypop, ypop_var = ypop_var))
}