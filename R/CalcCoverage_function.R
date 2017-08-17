CalcCoverage <- function(estimates, true_values, covariances = NULL, po_est = TRUE,
                         alpha_level = 0.05) {
  
  norm_quant <- qnorm(1 - alpha_level / 2)
  if (po_est) {
    cover <- array(NA, dim = dim(estimates))
    dimnames(cover) <- dimnames(estimates)
    for (it in c(0, 1)) {
      lb <- estimates[, it + 1, ] - norm_quant * sqrt(covariances[, 1, 1, ])
      ub <- estimates[, it + 1, ] + norm_quant * sqrt(covariances[, 2, 2, ])
      cover[, it + 1, ] <- (lb < true_values[it + 1, ] & ub > true_values[it + 1, ])
    }
    return(list(cover_indicators = cover,
                cover_prob = apply(cover, c(2, 3), mean, na.rm = TRUE)))
  } else {
    lb <- estimates[, 1, ] - norm_quant * sqrt(estimates[, 2, ])
    ub <- estimates[, 1, ] + norm_quant * sqrt(estimates[, 2, ])
    cover <- t(sapply(1:dim(estimates)[1],
                      function(x) lb[x, ] < true_values & ub[x, ] > true_values))
    dimnames(cover) <- dimnames(estimates)[- 2]
    return(list(cover_indicators = cover, cover_prob = colMeans(cover, na.rm = TRUE)))
  }
}