#' Direct effect estimates and asymptotic variance.
#' 
#' @param ypop A matrix with rows corresponding to the potential outcome under
#' control and treatment, and columns corresponding to the cluster-average
#' propensity of treatment.
#' @param ypop_var An 3-dimensional array, where the first two dimensions are
#' equal to 2 and include the variance covariance matrix of the population
#' average potential outcome for each alpha. Dimension 3 is alpha.
#' @param boots The results of BootVar() function including estimates of the
#' potential outcomes from the bootstrap samples.
#' @param alpha The values of alpha we consider. If ypop has column names,
#' alpha can be left null.
#' 
#' @return A matrix with rows including the estimate and variance of the direct
#' effect and columns corresponding to alpha.
#' 
#' @export
DE <- function(ypop, ypop_var, boots = NULL, alpha = NULL) {
  
  if (is.null(alpha)) {
    if (is.null(colnames(ypop))) {
      stop('Specify alpha.')
    }
    alpha <- as.numeric(colnames(ypop))
  }
  
  dim_names <- c('est', 'var', 'low_int', 'high_int')
  if (!is.null(boots)) {
    dim_names <- c(dim_names, 'boot_var', 'boot_low_int', 'boot_high_int')
  }
  
  de <- array(NA, dim = c(length(dim_names), length(alpha)))
  dimnames(de) <- list(stat = dim_names, alpha = alpha)

  de[1, ] <- ypop[2, ] - ypop[1, ]
  de[2, ] <- apply(ypop_var, 3, delta_method)
  de[3, ] <- de[1, ] - 1.96 * sqrt(de[2, ])
  de[4, ] <- de[1, ] + 1.96 * sqrt(de[2, ])
  
  if (!is.null(boots)) {
    de[5, ] <- apply(boots[2, , ] - boots[1, , ], 1, var, na.rm = TRUE)
    de[6, ] <- de[1, ] - 1.96 * sqrt(de[5, ])
    de[7, ] <- de[1, ] + 1.96 * sqrt(de[5, ])
  }
  
  return(de)
}