#' @param dta A data frame including the following variables: neigh varying from 1 to
#' the number of neighborhoods, and two covariates X1 and X2.
#' @param pot_out A list of length equal to the number of clusters. Each element is an
#' array of 3 dimensions. Element ijk of the array is the potential outcome for unit i
#' in the cluster if the treatment of unit is j - 1 and k - 1 of neighbors are treated.
#' @param neigh_ind A list of length equal to the number of custers including the row
#' indices corresponding to each cluster.
#' @param re_sd The standard deviation of the cluster specific random effect for the
#' propensity score model.
#' @param trt_coef Coefficients of the propensity score model corresponding to the
#' intercept, and the two covariates in dta.
#' 
#' @export
GetSimData <- function(dta, pot_out, neigh_ind, re_sd, trt_coef) {
  
  n_obs <- nrow(dta)
  n_neigh <- length(unique(dta$neigh))
  num_cov <- length(trt_coef) - 1
  cov_cols <- which(names(dta) %in% paste0('X', 1 : num_cov))
  
  # Generating the treatment indicators.
  re <- rnorm(n_neigh, mean = 0, sd = re_sd)
  des_mat <- as.matrix(cbind(1, dta[, cov_cols]))
  probs <- expit(des_mat %*% trt_coef + re[dta$neigh])
  print(paste('Range of treatment probability:',
              paste(round(range(probs), 3), collapse = '-')))
  A <- rbinom(n_obs, 1, prob = probs)
  
  # Choosing the observed Y from the matrix of potential outcomes.
  Y <- rep(NA, n_obs)
  for (ii in 1:n_obs) {
    wh_neigh <- dta$neigh[ii]
    wh_obs <- which(neigh_ind[[wh_neigh]] == ii)
    sum_others_trt <- sum(A[setdiff(neigh_ind[[wh_neigh]], ii)])  # Treated neighbors.
    Y[ii] <- pot_out[[wh_neigh]][wh_obs, A[ii] + 1, sum_others_trt + 1]
  }
  
  sim_dta <- cbind(dta, A = A, Y = Y, trt_prob = probs)
  return(list(data = sim_dta, cluster_re = re))
}