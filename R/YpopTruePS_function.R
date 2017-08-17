YpopTruePS <- function(ygroup, alpha, use = 'pairwise.complete.obs') {
  
  ypop <- apply(ygroup, c(2, 3), mean, na.rm = TRUE)
  n_neigh <- dim(ygroup)[1]
  
  ypop_var <- apply(ygroup, 3, function(x) cov(x, use = use))
  ypop_var <- array(ypop_var, dim = c(2, 2, length(alpha)))
  # In order to get 1 / N, instead of 1 / (N - 1) in the variance estimates.
  ypop_var <- ypop_var * (n_neigh - 1) / n_neigh
  ypop_var <- ypop_var / n_neigh  # Since we have n_neigh clusters.
  
  return(list(ypop = ypop, ypop_var = ypop_var))
}