#' Indirect effect estimates and asymptotic variance.
#' 
#' @param ypop A vector of the population average potential outcome for the
#' values of alpha.
#' @param ypop_var An array with dimensions 2, 2, alpha, alpha, including the
#' covariance matrix of the average potential outcome under two specifications
#' of alpha.
#' @param alpha The values of alpha we consider.
#' 
#' @return A 3-dimensional array with dimension 1 equal to 2 including the
#' estimate and the asymptotic variance. The remaining two dimensions
#' correspond to values of alpha.
#' 
#' @export
IE <- function(ypop, ypop_var, alpha) {
  
  ie <- array(NA, dim = c(4, length(alpha), length(alpha)))
  dimnames(ie) <- list(stat = c('est', 'var', 'LB', 'UB'), alpha1 = alpha,
                       alpha2 = alpha)
  
  for (a1 in 1 : length(alpha)) {
    for (a2 in 1 : length(alpha)) {
      ie[1, a1, a2] <- ypop[a1] - ypop[a2]
      ie[2, a1, a2] <- delta_method(ie_var[c(a1, a2), c(a1, a2)])
      ie[c(3, 4), a1, a2] <- ie[1, a1, a2] + 1.96 * c(- 1, 1) * sqrt(ie[2, a1, a2])
    }
  }

  return(ie)
}