#' Function that estimates the group average potential outcome using IPW.
#'
#' @param dta The dataset as a data frame including treatment, outcome and covariates.
#' @param cov_cols The indeces including the covariates of the propensity score model.
#' @param phi_hat A list with two elements. The first one is a vector of coefficients
#' of the propensity score, and the second one is the random effect variance.
#' @param gamma_numer The coefficients of the propensity score model in the numerator.
#' If left NULL, the coefficients in phi_hat will be used instead.
#' @param alpha The values of alpha for which we want to estimate the group average
#' potential outcome.
#' @param neigh_ind List. i^{th} element is a vector with the row indeces of dta that
#' are in cluster i. Can be left NULL.
#' @param trt_col If the treatment is not named 'A' in dta, specify the treatment
#' column index.
#' @param out_col If the outcome is not named 'Y', specify the outcome column index.
#' @param lower The lower end of the values for bi we will look at. Defaults to - 10.
#' @param upper The upper end of the values for bi we will look at. Defaults to 10.
#' @param integral_bound The number of standard deviations of the random effect that
#' will be used as the lower and upper limit.
#' @param keep_re_alpha Logical. If set to TRUE the "random" effect that makes the
#' average probability of treatment equal to alpha will be returned along with the
#' estimated group average potential outcome.
GroupIPW <- function(dta, cov_cols, phi_hat, gamma_numer = NULL, alpha,
                     neigh_ind = NULL, trt_col = NULL, out_col = NULL, lower = - 10,
                     upper = 10, integral_bound = 10, keep_re_alpha = FALSE) {
  
  integral_bound <- abs(integral_bound)
  
  dta <- as.data.frame(dta)
  n_neigh <- length(neigh_ind)
  
  yhat_group <- array(NA, dim = c(n_neigh, 2, length(alpha)))
  dimnames(yhat_group) <- list(neigh = 1:n_neigh, trt = c(0, 1), alpha = alpha)
  
  # Names of treatment and outcome column.
  if (!is.null(trt_col)) {
    names(dta)[trt_col] <- 'A'
  }
  if (!is.null(trt_col)) {
    names(dta)[out_col] <- 'Y'
  }
  if (is.null(gamma_numer)) {
    gamma_numer <- matrix(phi_hat[[1]], ncol = 1)
  }
  
  # Specifyling neigh_ind will avoid re-running the following lines.
  if (is.null(neigh_ind)) {
    for (nn in 1:n_neigh) {
      neigh_ind[[nn]] <- which(dta$neigh == nn)
    }
  }
  
  if (keep_re_alpha) {
    re_alphas <- lapply(1 : n_neigh, function(nn) {
      A <- array(NA, dim = c(length(neigh_ind[[nn]]), 2, length(alpha)))
      dimnames(A) <- list(unit = 1:(dim(A)[1]), a = c(0, 1), alpha = alpha)
      return(A)
    })
  }
  
  for (aa in 1:length(alpha)) {
    curr_alpha <- alpha[[aa]]

    for (it in c(0, 1)) {
      curr_it <- it
      bern_prob <- curr_alpha ^ curr_it * (1 - curr_alpha) ^ (1 - curr_it)

      for (nn in 1:n_neigh) {
        # If no individuals have treatement it we cannot estimate the group average.
        y_curr <- ifelse(sum(dta$A[neigh_ind[[nn]]] == curr_it) == 0, NA, 0)
        for (ind in neigh_ind[[nn]]) {
          if (dta$A[ind] == curr_it) {
            
            wh_others <- setdiff(neigh_ind[[nn]], ind)
            Ai_j <- dta$A[wh_others]
            Xi_j <- dta[wh_others, cov_cols]
            
            # For the observed values of the covariates.
            prob_ind <- CalcNumerator(Ai_j = Ai_j, Xi_j = Xi_j,
                                      coef_hat = gamma_numer,
                                      alpha = curr_alpha,
                                      lower = lower, upper = upper)
            if (keep_re_alpha) {
              obs <- which(neigh_ind[[nn]] == ind)
              re_alphas[[nn]][obs, curr_it + 1, aa] <- prob_ind$re_alpha
            }
            
            y_curr <- y_curr + dta$Y[ind] * prob_ind$prob
          }
        }
        
        denom <- DenomIntegral(A = dta$A[neigh_ind[[nn]]],
                               X = dta[neigh_ind[[nn]], cov_cols],
                               phi_hat = phi_hat, alpha = curr_alpha,
                               integral_bound = integral_bound)
        denom <- length(neigh_ind[[nn]]) * denom$value * bern_prob
        
        yhat_group[nn, it + 1, aa] <- y_curr / denom
      }
    }
  }
  if (keep_re_alpha) {
    return(list(yhat_group = yhat_group, re_alpha = re_alphas))
  }
  return(list(yhat_group = yhat_group))
}