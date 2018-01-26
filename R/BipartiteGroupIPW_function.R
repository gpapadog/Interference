#' Function that estimates the cluster average potential outcome in settings of
#' bipartite interference.
#'
#' @param int_dta Data on the interventional units including treatment,
#' neighborhood index and covariates. Neighborhood variable should be named
#' 'neigh' and only include continuous integers from 1 to the number of
#' neighborhoods.
#' @param out_dta Data on the outcome units including outcome, neighborhood,
#' and index of the closest interventional unit as 'closest_int'.
#' Neighborhood variable should be named 'neigh' and only include continuous
#' integers from 1 to the number of neighborhoods.
#' @param cov_cols The indices including the covariates of the ps model.
#' @param phi_hat A list with two elements. The first one is a vector of
#' coefficients of the ps, and the second one is the random effect variance.
#' @param alpha The values of alpha for which we want to estimate the group
#' average potential outcome.
#' @param trt_col If the treatment is not named 'A' in int_dta, specify the
#' treatment column index.
#' @param out_col If the outcome is not named 'Y' in out_dta, specify the
#' outcome column index.
#' @param integral_bound The number of standard deviations of the random effect
#' that will be used as the lower and upper limit. Defaults to 10.
#' 
#' @export
BipartiteGroupIPW <- function(int_dta, out_dta, cov_cols, phi_hat, alpha,
                              trt_col = NULL, out_col = NULL,
                              integral_bound = 10) {
  
  integral_bound <- abs(integral_bound)
  phi_hat[[1]] <- matrix(phi_hat[[1]], ncol = 1)
  int_dta <- as.data.frame(int_dta)
  out_dta <- as.data.frame(out_dta)
  n_neigh <- max(int_dta$neigh)
  
  # Names of treatment and outcome column.
  if (!is.null(trt_col)) {
    names(int_dta)[trt_col] <- 'A'
  }
  if (!is.null(trt_col)) {
    names(out_dta)[out_col] <- 'Y'
  }
  
  # Getting lists of which observations belong to which cluster.  
  int_neigh_ind <- sapply(1 : n_neigh, function(x) which(int_dta$neigh == x))
  out_neigh_ind <- sapply(1 : n_neigh, function(x) which(out_dta$neigh == x))
  
  # Where the results will be saved.
  yhat_group <- array(NA, dim = c(n_neigh, 2, length(alpha)))
  dimnames(yhat_group) <- list(neigh = 1 : n_neigh, trt = c(0, 1),
                               alpha = alpha)
  
  # Adding treatment of closest interventional in out_dta.
  out_dta$closest_trt <- int_dta$A[as.numeric(out_dta$closest_int)]
  
  for (aa in 1 : length(alpha)) {
    
    print(paste('alpha =', round(alpha[aa], 3)))
    curr_alpha <- alpha[[aa]]
    
    for (nn in 1 : n_neigh) {
      for (curr_it in c(0, 1)) {
        
        bern_prob <- curr_alpha ^ curr_it * (1 - curr_alpha) ^ (1 - curr_it)
        
        # How many closest power plants have the treatment.
        num_curr_it <- sum(out_dta$closest_trt[out_neigh_ind[[nn]]] == curr_it)
        # If none of the closest units have the treatment, don't estimate.
        y_curr <- ifelse(num_curr_it == 0, NA, 0)
        
        for (ind in out_neigh_ind[[nn]]) {
          if (out_dta$closest_trt[ind] == curr_it) {
            y_curr <- y_curr + out_dta$Y[ind]
          }
        }
        
        denom <- DenomIntegral(A = int_dta$A[int_neigh_ind[[nn]]],
                               X = int_dta[int_neigh_ind[[nn]], cov_cols],
                               phi_hat = phi_hat, alpha = curr_alpha,
                               integral_bound = integral_bound)
        denom <- length(out_neigh_ind[[nn]]) * denom$value * bern_prob
        
        yhat_group[nn, curr_it + 1, aa] <- y_curr / denom
      }
    }
  }
  return(yhat_group)
}