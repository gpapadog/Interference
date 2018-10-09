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
#' estimated propensity score. The list includes 1) glm_form: Element of 
#' formula class. The formula can be either for a fixed effects model or for a
#' model including random intercepts. 2) ps_with_re: An indicator of whether
#' the propensity score is a mixed model (set to TRUE) or not (set to FALSE).
#' 3) gamma_numer: The coefficients of the covariates in the counterfactual
#' treatment allocation model, and 4) use_control: Set to TRUE or FALSE if
#' you want or do not want additional elements in fitting the mixed model.
#' use_control does not have to be specified.
#' @param verbose Logical. Whether progress is printed. Defaults to TRUE.
#' @param trt_col If the treatment is not named 'A' in dta, specify the
#' treatment column index.
#' @param out_col If the outcome is not named 'Y', specify the outcome column
#' index.
#' 
#' @export
BootVar <- function(dta, B = 500, alpha, ps = c('true', 'est'), cov_cols,
                    phi_hat_true = NULL, ps_info_est = NULL, verbose = TRUE,
                    ps_specs = NULL, trt_col = NULL, out_col = NULL) {
  
  ps <- match.arg(ps)
  
  boots <- array(NA, dim = c(2, length(alpha), B))
  dimnames(boots) <- list(po = c('y0', 'y1'), alpha = alpha, sample = 1 : B)
  
  for (bb in 1 : B) {
    
    if (verbose) {
      if (bb %% 10 == 0) {
        print(paste0('bootstrap sample ', bb))
      }
    }
    
    boot_dta <- GetBootSample(dta)
    neigh_ind <- lapply(1 : max(boot_dta$neigh),
                        function(nn) which(boot_dta$neigh == nn))
    
    if (ps == 'true') { # Known propensity score.
      
      ygroup_boot <- GroupIPW(dta = boot_dta, cov_cols = cov_cols,
                              phi_hat = phi_hat_true, gamma_numer = NULL,
                              alpha = alpha, neigh_ind = neigh_ind,
                              keep_re_alpha = FALSE, estimand = '1',
                              verbose = FALSE, trt_col = trt_col,
                              out_col = out_col)$yhat_group
      boots[, , bb] <- apply(ygroup_boot, c(2, 3), mean) 
      
    } else {  # Estimated propensity score.
      
      if (ps_info_est$ps_with_re) {  # The PS model includes random intercepts.
        
        if (!is.null(ps_info_est$use_control) & ps_info_est$use_control) {
          glm_control <- glmerControl(optimizer = "bobyqa",
                                      optCtrl = list(maxfun = 2e5))
          glmod <- lme4::glmer(ps_info_est$glm_form, data = boot_dta,
                               family = binomial, control = glm_control)
        } else {
          glmod <- lme4::glmer(ps_info_est$glm_form, data = boot_dta,
                               family = binomial)
        }
        re_var <- as.numeric(summary(glmod)$varcor)
        
        
      } else {  # The PS model includes only fixed effects.
        
        glmod <- glm(ps_info_est$glm_form, data = boot_dta, family = binomial)
        re_var <- 0
        
      }
      
      phi_hat_est <- list(coefs = summary(glmod)$coef[, 1], re_var = re_var)
      ygroup_boot <- GroupIPW(dta = boot_dta, cov_cols = cov_cols,
                              phi_hat = phi_hat_est,
                              gamma_numer = ps_info_est$gamma_numer,
                              alpha = alpha, neigh_ind = neigh_ind,
                              keep_re_alpha = FALSE, estimand = '1',
                              verbose = FALSE, trt_col = trt_col,
                              out_col = out_col)$yhat_group
      boots[, , bb] <- apply(ygroup_boot, c(2, 3), mean)
    }
  }
  return(boots)
  
}

