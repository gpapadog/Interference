#' Direct effect estimates and asymptotic variance.
#' 
#' @param ypop A matrix with rows corresponding to the potential outcome under
#' control and treatment, and columns corresponding to the cluster-average
#' propensity of treatment.
#' @param ypop_var An 3-dimensional array, where the first two dimensions are
#' equal to 2 and include the variance covariance matrix of the population
#' average potential outcome for each alpha. Dimension 3 is alpha.
#' @param alpha The values of alpha we consider. If ypop has column names,
#' alpha can be left null.
#' 
#' @return A matrix with rows including the estimate and variance of the direct
#' effect and columns corresponding to alpha.
#' 
#' @export
DE <- function(ypop, ypop_var, alpha = NULL) {
  
  if (is.null(alpha)) {
    if (is.null(colnames(ypop))) {
      stop('Specify alpha.')
    }
    alpha <- as.numeric(colnames(ypop))
  }
  
  de <- array(NA, dim = c(2, length(alpha)))
  dimnames(de) <- list(stat = c('est', 'var'), alpha = alpha)

  de[1, ] <- ypop[2, ] - ypop[1, ]
  de[2, ] <- apply(ypop_var, 3, delta_method)
  
  return(de)
}