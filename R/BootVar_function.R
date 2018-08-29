#' Bootstrap variance of potential outcomes.
#' 
#' Using re-sampling of clusters to acquire an estimate of the potential
#' outcome estimator variance.
#' 
#' @param dta The data frame including the observed data set.
#' @param B Number of bootstrap samples.
#' @param alpha The values of alpha where the potential outcomes are estimated.
#' @param ps Character. Whether the propensity score is known or is estimated.
#' Options include 'true', 'est'. Defaults to 'true'.
#' @param cov_cols Vector of column indices of the covariates used in the
#' propensity score model.
#' @param phi_hat_true List. Specify if the propensity score is known (ps set
#' to 'true'). Elements of the list are trt_coef and re_var including the
#' coefficients of the propensity score and random effect variance.
#' @param ps_info_est List of elements for acquiring estimates based on the
#' estimated propensity score. The list includes 1) glm_form including an
#' element of formula class including the random intercept formula of the
#' propensity score model. 2) trt_coef The coefficients of the covariates in
#' the counterfactual treatment allocation model.
#' @param verbose Logical. Whether progress is printed. Defaults to TRUE.
#' 
BootVar <- function(dta, B = 500, alpha, ps = c('true', 'est'), cov_cols,
                    phi_hat_true = NULL, ps_info_est = NULL, verbose = TRUE) {
  
  num_clus <- max(dta$neigh)
  ps <- match.arg(ps)
  
  boots <- array(NA, dim = c(2, length(alpha), B))
  dimnames(boots) <- list(po = c('y0', 'y1'), alpha = alpha, sample = 1 : B)
  re_var_positive <- rep(FALSE, B)
  
  for (bb in 1 : B) {
    
    if (verbose) {
      if (bb %% 10 == 0) {
        print(paste0('bootstrap sample ', bb))
      }
    }
    
    boot_clusters <- sample(1 : num_clus, num_clus, replace = TRUE)
    
    # Binding data without accidentally merging repeated clusters.
    boot_dta <- NULL
    for (nn in 1 : num_clus) {
      D <- subset(dta, neigh == boot_clusters[nn])
      D$neigh <- nn
      boot_dta <- rbind(boot_dta, D)
    }
    neigh_ind <- lapply(1 : num_clus, function(nn) which(boot_dta$neigh == nn))
    
    if (ps == 'true') {
      
      # Known propensity score.
      ygroup_boot <- GroupIPW(dta = boot_dta, cov_cols = cov_cols,
                              phi_hat = phi_hat_true, gamma_numer = NULL,
                              alpha = alpha, neigh_ind = neigh_ind,
                              keep_re_alpha = FALSE, estimand = '1',
                              verbose = FALSE)$yhat_group
      ygroup_boot[is.na(ygroup_boot)] <- 0
      boots[, , bb] <- apply(ygroup_boot, c(2, 3), mean) 
      
    } else {
      
      glmod <- lme4::glmer(ps_info_est$glm_form, data = boot_dta,
                           family = binomial)
      re_var <- as.numeric(summary(glmod)$varcor)
      
      if (re_var > 0) {
        re_var_positive[bb] <- TRUE
        phi_hat_est <- list(coefs = summary(glmod)$coef[, 1], re_var = re_var)
        ygroup_boot <- GroupIPW(dta = boot_dta, cov_cols = cov_cols,
                                phi_hat = phi_hat_est,
                                gamma_numer = ps_info_est$trt_coef,
                                alpha = alpha, neigh_ind = neigh_ind,
                                keep_re_alpha = FALSE, estimand = '1',
                                verbose = FALSE)$yhat_group
        ygroup_boot[is.na(ygroup_boot)] <- 0
        boots[, , bb] <- apply(ygroup_boot, c(2, 3), mean)
      }
    }
  }
  if (ps == 'true') {
    return(list(boots = boots))
  }
  return(list(boots = boots, re_var_positive = re_var_positive))
}