#' Estimates the group average potential outcome using IPW with a PS model without RE.
#'
#' @param dta The dataset as a data frame including treatment, outcome and covariates.
#' @param cov_cols The indeces including the covariates of the propensity score model.
#' @param coef_hat A vector of coefficients of the propensity score.
#' @param alpha The values of alpha for which we want to estimate the group average
#' potential outcome.
#' @param neigh_ind List. i^{th} element is a vector with the row indeces of dta that
#' are in cluster i. Can be left NULL.
#' @param trt_col If the treatment is not named 'A' in dta, specify the treatment
#' column index.
#' @param out_col If the outcome is not named 'Y', specify the outcome column index.
#' @param lower The lower end of the values for bi we will look at. Defaults to - 10.
#' @param upper The upper end of the values for bi we will look at. Defaults to 10.
#' 
GroupIPW_woRE <- function(dta, cov_cols, coef_hat, alpha, neigh_ind, trt_col = NULL,
                          out_col = NULL, lower = - 10, upper = 10) {
  
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
  
  # Specifyling neigh_ind will avoid re-running the following lines.
  if (is.null(neigh_ind)) {
    for (nn in 1:n_neigh) {
      neigh_ind[[nn]] <- which(dta$neigh == nn)
    }
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
                                      coef_hat = coef_hat,
                                      alpha = curr_alpha,
                                      lower = lower, upper = upper)
            y_curr <- y_curr + dta$Y[ind] * prob_ind
          }
        }
        
        denom <- Denom_woRE(A = dta$A[neigh_ind[[nn]]],
                            X = dta[neigh_ind[[nn]], cov_cols],
                            coef_hat = coef_hat, alpha = curr_alpha)
        denom <- length(neigh_ind[[nn]]) * denom$value * bern_prob
        
        yhat_group[nn, it + 1, aa] <- y_curr / denom
      }
    }
  }
  return(yhat_group)
}