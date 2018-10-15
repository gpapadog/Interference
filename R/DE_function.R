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
#' @param alpha_level Numeric. The alpha level of the confidence intervals
#' based on the quantiles of the bootstrap estimates.
#' 
#' @return A matrix with rows including the estimate and variance of the direct
#' effect and columns corresponding to alpha.
#' 
#' @export
DE <- function(ypop, ypop_var, boots = NULL, alpha = NULL,
               alpha_level = 0.05) {
  
  quants <- c(0, 1) + c(1, - 1) * alpha_level / 2
  norm_quant <- - qnorm(alpha_level / 2)
  
  if (is.null(alpha)) {
    if (is.null(colnames(ypop))) {
      stop('Specify alpha.')
    }
    alpha <- as.numeric(colnames(ypop))
  }
  
  dim_names <- c('est', 'var', 'low_int', 'high_int')
  if (!is.null(boots)) {
    dim_names <- c(dim_names, 'boot_var', 'boot_var_LB', 'boot_var_UB',
                   'boot_low_quant', 'boot_high_quant')
  }
  
  de <- array(NA, dim = c(length(dim_names), length(alpha)))
  dimnames(de) <- list(stat = dim_names, alpha = alpha)

  de[1, ] <- ypop[2, ] - ypop[1, ]
  de[2, ] <- apply(ypop_var, 3, delta_method)
  de[3, ] <- de[1, ] - norm_quant * sqrt(de[2, ])
  de[4, ] <- de[1, ] + norm_quant * sqrt(de[2, ])
  
  if (!is.null(boots)) {
    de[5, ] <- apply(boots[2, , ] - boots[1, , ], 1, var)
    de[6, ] <- de[1, ] - norm_quant * sqrt(de[5, ])
    de[7, ] <- de[1, ] + norm_quant * sqrt(de[5, ])
    de[8 : 9, ] <- apply(boots[2, , ] - boots[1, , ], 1, quantile,
                         probs = quants)
  }
  
  return(de)
}