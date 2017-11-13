#' Indirect effect estimates and asymptotic variance.
#' 
#' @param ygroup An matrix including the group average potential outcome
#' estimates where rows correspond to group, and columns to values of alpha.
#' @param ps String. Can take values 'true', or 'estimated' for known or
#' estimated propensity score. Defaults to 'true'.
#' @param scores A matrix with rows corresponding to the parameters of the
#' propensity score model and columns for groups. Includes the score of the
#' propensity score evaluated for the variables of each group. Can be left
#' NULL for ps = 'true'.
#' 
#' @export
IE <- function(ygroup, ps = c('true', 'estimated'), scores = NULL) {
  
  alpha <- as.numeric(dimnames(ygroup)[[2]])
  ypop <- apply(ygroup, 2, mean)
  names(ypop) <- alpha
  
  ie_var <- IEvar(ygroup = ygroup, ps = ps, scores = scores)
  
  ie <- array(NA, dim = c(4, length(alpha), length(alpha)))
  dimnames(ie) <- list(stat = c('est', 'var', 'LB', 'UB'), alpha1 = alpha,
                       alpha2 = alpha)
  
  for (a1 in 1 : length(alpha)) {
    for (a2 in 1 : length(alpha)) {
      ie[1, a1, a2] <- ypop[a2] - ypop[a1]
      ie[2, a1, a2] <- delta_method(ie_var[c(a1, a2), c(a1, a2)])
      ie[c(3, 4), a1, a2] <- ie[1, a1, a2] + 1.96 * c(- 1, 1) * sqrt(ie[2, a1, a2])
    }
  }

  return(ie)
}