#' @param dta A data frame including the following variables: neigh varying from 1 to
#' the number of neighborhoods, and two covariates X1 and X2.
#' @param pot_out A list of length equal to the number of clusters. Each element is an
#' array of 3 dimensions. Element ijk of the array is the potential outcome for unit i
#' in the cluster if the treatment of unit is j - 1 and k - 1 of neighbors are treated.
#' @param neigh_ind A list of length equal to the number of custers including the row
#' indices corresponding to each cluster.
#' @param trt_coef Coefficients of the propensity score model corresponding to the
#' intercept, and the two covariates in dta.
GetSimData_woRE <- function(dta, pot_out, neigh_ind, trt_coef) {
  
  n_obs <- nrow(dta)
  n_neigh <- length(unique(dta$neigh))
  
  # Generating the treatment indicators.
  probs <- expit(cbind(1, dta$X1, dta$X2) %*% trt_coef)
  A <- rbinom(n_obs, 1, prob = probs)
  
  # Choosing the observed Y from the matrix of potential outcomes.
  Y <- rep(NA, n_obs)
  for (ii in 1:n_obs) {
    wh_neigh <- dta$neigh[ii]
    wh_obs <- which(neigh_ind[[wh_neigh]] == ii)
    sum_others_trt <- sum(A[setdiff(neigh_ind[[wh_neigh]], ii)])  # Treated neighbors.
    Y[ii] <- pot_out[[wh_neigh]][wh_obs, A[ii] + 1, sum_others_trt + 1]
  }
  
  sim_dta <- cbind(dta, A = A, Y = Y)
  return(sim_dta)
}