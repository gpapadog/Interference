#' Function that estimates the group average potential outcome using IPW.
#'
#' @param dta Data frame including treatment, outcome and covariates.
#' @param cov_cols The indices including the covariates of the ps model.
#' @param phi_hat A list with two elements. The first one is a vector of
#' coefficients of the ps, and the second one is the random effect variance.
#' @param gamma_numer The coefficients of the ps model in the numerator.
#' If left NULL and estimand is 1, the coefficients in phi_hat will be used
#' instead.
#' @param alpha The values of alpha for which we want to estimate the group
#' average potential outcome.
#' @param neigh_ind List. i^{th} element is a vector with the row indices of
#' dta that are in cluster i. Can be left NULL.
#' @param trt_col If the treatment is not named 'A' in dta, specify the
#' treatment column index.
#' @param out_col If the outcome is not named 'Y', specify the outcome column
#' index.
#' @param alpha_re_bound The lower and upper end of the values for bi we will
#' look at. Defaults to 10, meaning we will look between - 10 and 10.
#' @param integral_bound The number of standard deviations of the random effect
#' that will be used as the lower and upper limit.
#' @param keep_re_alpha Logical. If set to TRUE the "random" effect that makes
#' the average probability of treatment equal to alpha will be returned along
#' with the estimated group average potential outcome.
#' @param estimand Character, either '1' or '2.' If 1 is specified, then the
#' estimand with numerator depending on covariates is estimated. If estimand
#' is set equal to 2, the numerator considered is the product of Bernoulli.
#' 
#' @export
GroupIPW <- function(dta, cov_cols, phi_hat, gamma_numer = NULL, alpha,
                     neigh_ind = NULL, trt_col = NULL, out_col = NULL, 
                     alpha_re_bound = 10, integral_bound = 10,
                     keep_re_alpha = FALSE, estimand = c('1', '2')) {
  
  estimand <- match.arg(estimand)
  integral_bound <- abs(integral_bound)
  alpha_re_bound <- abs(alpha_re_bound)
  phi_hat[[1]] <- matrix(phi_hat[[1]], ncol = 1)
  dta <- as.data.frame(dta)
  
  # We only return the ksi's if we are estimating estimand 1.
  keep_re_alpha <- keep_re_alpha & (estimand == '1')
  
  # Specifyling neigh_ind will avoid re-running the following lines.
  if (is.null(neigh_ind)) {
    neigh_ind <- sapply(1 : max(dta$neigh), function(x) which(dta$neigh == x))
  }
  
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
  
  # If we want to return the ksis that make average propensity alpha.
  if (keep_re_alpha) {
    re_alphas <- matrix(NA, nrow = n_neigh, ncol = length(alpha))
    dimnames(re_alphas) <- list(neigh = 1 : n_neigh, alpha = alpha)
  }
  
  for (aa in 1 : length(alpha)) {
    print(paste('alpha =', round(alpha[aa], 3)))
    curr_alpha <- alpha[[aa]]
    
    for (nn in 1 : n_neigh) {
      
      # For estimand 1, we need to calculate numerator depending on covariates.
      if (estimand == '1') {
        
        # Calculating the random effect that gives alpha.
        Xi <- dta[neigh_ind[[nn]], cov_cols]
        lin_pred <- cbind(1, as.matrix(Xi)) %*% phi_hat[[1]]
        re_alpha <- FromAlphaToRE(alpha = curr_alpha, lin_pred = lin_pred,
                                  alpha_re_bound = alpha_re_bound)
        
        # Keeping the intercept that makes cluster average propensity alpha.
        if (keep_re_alpha) {
          re_alphas[nn, aa] <- re_alpha
        }
        
      }

      for (curr_it in c(0, 1)) {
        
        bern_prob <- curr_alpha ^ curr_it * (1 - curr_alpha) ^ (1 - curr_it)
        prob_ind <- list(prob = 1)  # For estimand 2.
        
        # If no individuals have treatement it we cannot estimate the group average.
        y_curr <- ifelse(sum(dta$A[neigh_ind[[nn]]] == curr_it) == 0, NA, 0)
        
        for (ind in neigh_ind[[nn]]) {
          if (dta$A[ind] == curr_it) {
            
            wh_others <- setdiff(neigh_ind[[nn]], ind)
            Ai_j <- dta$A[wh_others]
            Xi_j <- dta[wh_others, cov_cols]
            
            if (estimand == '1') {
              prob_ind <- CalcNumerator(Ai_j = Ai_j, Xi_j = Xi_j,
                                        coef_hat = gamma_numer,
                                        alpha = curr_alpha, re_alpha = re_alpha)
            }
            y_curr <- y_curr + dta$Y[ind] * prob_ind$prob
          }
        }
        
        denom <- DenomIntegral(A = dta$A[neigh_ind[[nn]]],
                               X = dta[neigh_ind[[nn]], cov_cols],
                               phi_hat = phi_hat, alpha = curr_alpha,
                               integral_bound = integral_bound)
        denom <- length(neigh_ind[[nn]]) * denom$value * bern_prob
        
        yhat_group[nn, curr_it + 1, aa] <- y_curr / denom
      }
    }
  }
  if (keep_re_alpha) {
    return(list(yhat_group = yhat_group, re_alpha = re_alphas))
  }
  return(list(yhat_group = yhat_group))
}