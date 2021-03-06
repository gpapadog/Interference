#' Indirect effect estimates and asymptotic variance.
#' 
#' @param ygroup An matrix including the group average potential outcome
#' estimates where rows correspond to group, and columns to values of alpha.
#' #' @param boots The results of BootVar() function including estimates of the
#' potential outcomes from the bootstrap samples.
#' @param ps String. Can take values 'true', or 'estimated' for known or
#' estimated propensity score. Defaults to 'true'.
#' @param scores A matrix with rows corresponding to the parameters of the
#' propensity score model and columns for groups. Includes the score of the
#' propensity score evaluated for the variables of each group. Can be left
#' NULL for ps = 'true'.
#' @param alpha_level Numeric. The alpha level of the confidence intervals
#' based on the quantiles of the bootstrap estimates.
#' 
#' @export
IE <- function(ygroup, boots = NULL, ps = c('true', 'estimated'),
               scores = NULL, alpha_level = 0.05) {
  
  alpha <- as.numeric(dimnames(ygroup)[[2]])
  ypop <- apply(ygroup, 2, mean)
  names(ypop) <- alpha
  quants <- c(0, 1) + c(1, - 1) * alpha_level / 2
  norm_quant <- - qnorm(alpha_level / 2)
  
  ie_var <- IEvar(ygroup = ygroup, ps = ps, scores = scores)
  dim_names <- c('est', 'var', 'LB', 'UB')
  if (!is.null(boots)) {
    dim_names <- c(dim_names, 'boot_var', 'boot_var_LB', 'boot_var_UB',
                   'boot_low_quant', 'boot_high_quant')
  }
  
  ie <- array(NA, dim = c(length(dim_names), length(alpha), length(alpha)))
  dimnames(ie) <- list(stat = dim_names, alpha1 = alpha, alpha2 = alpha)
  
  for (a1 in 1 : length(alpha)) {
    for (a2 in 1 : length(alpha)) {
      ie[1, a1, a2] <- ypop[a2] - ypop[a1]
      ie[2, a1, a2] <- delta_method(ie_var[c(a1, a2), c(a1, a2)])
      ie_sd <- sqrt(ie[2, a1, a2])
      ie[c(3, 4), a1, a2] <- ie[1, a1, a2] + norm_quant * c(- 1, 1) * ie_sd
    }
  }

  if (!is.null(boots)) {
    ie_var_boots <- array(NA, dim = c(length(alpha), length(alpha)))
    for (a1 in 1 : length(alpha)) {
      for (a2 in 1 : length(alpha)) {
        ie[5, a1, a2] <- var(boots[1, a1, ] - boots[1, a2, ])
        ie_sd <- sqrt(ie[5, a1, a2])
        ie[c(6, 7), a1, a2] <- ie[1, a1, a2] + norm_quant * c(- 1, 1) * ie_sd
        ie[c(8, 9), a1, a2] <- quantile(boots[1, a2, ] - boots[1, a1, ],
                                        probs = quants)
      }
    }
  }
  
  return(ie)
}