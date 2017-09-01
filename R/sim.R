# --------------------
# Change when running on the cluster.
setwd('~/Github/Interference/R/')
source('~/Documents/Functions/expit_function.R')
source('DenomIntegral_function.R')
source('CalcNumerator_function.R')
source('FromAlphaToRE_function.R')
source('GroupIPW_function.R')
source('GroupLikelihood_function.R')
source('GetSimData_function.R')
source('YpopTruePS_function.R')
source('VarEstPS_function.R')
source('CalcScore_function.R')
source('CalcB11matrix_function.R')

setwd('~/Documents/Interference/Simulations/')
load_path <- 'Population_quantities/Data/Data3/'
out_path <- NULL
library(lme4)
library(numDeriv)
# --------------------


# --------------------
#  Estimation - alpha.
alpha <- seq(0.2, 0.8, by = 0.05)
a <- c(0, 1)
lower <- - 10
upper <- 10
# --------------------


# Loading the covariates.
load(paste0(load_path, 'covariates.dat'))
load(paste0(load_path, 'potential_outcomes.dat'))
load(paste0(load_path, 'data_specs.dat'))

# --------------------
# Setting up data generation.
re_sd <- data_specs$re_sd  # Can be changed. Data generation does not depend on it.
trt_coef <- data_specs$trt_fixed_effects
# --------------------

n_neigh <- max(dta$neigh)
neigh_ind <- lapply(1 : n_neigh, function(nn) which(dta$neigh == nn))

sim_dta <- GetSimData(dta = dta, pot_out = y, neigh_ind = neigh_ind, re_sd = re_sd,
                      trt_coef = trt_coef)$data
cov_cols <- c(which(names(sim_dta) == 'X1'), which(names(sim_dta) == 'X2'))

glmod <- glmer(A ~ X1 + X2 + (1 | neigh), data = sim_dta, family = binomial)

# What is the observed proportion of treated in each cluster?
obs_alpha <- sapply(1 : n_neigh, function(nn) {
D <- sim_dta[neigh_ind[[nn]], ]
return(mean(D[, which(names(D) == 'A')]))
})

sim_dta <- subset(sim_dta, ! (neigh %in% which(obs_alpha %in% c(0, 1))))
sim_dta$neigh <- as.numeric(as.factor(sim_dta$neigh))
n_neigh <- max(sim_dta$neigh)
neigh_ind <- lapply(1 : n_neigh, function(nn) which(sim_dta$neigh == nn))

#  -------------------------------------------
#   ----  FOR THE TRUE PROPENSITY SCORE ----
#  -------------------------------------------

# ---- Calculating the ipw, asymptotic variance.

phi_hat_true <- list(coefs = trt_coef, re_var = re_sd ^ 2)
ygroup_true <- GroupIPW(dta = sim_dta, cov_cols = cov_cols, phi_hat = phi_hat_true,
                        gamma_numer = NULL, alpha = alpha, neigh_ind = neigh_ind,
                        keep_re_alpha = TRUE)
re_alpha_true <- ygroup_true$re_alpha
ygroup_true <- ygroup_true$yhat_group

ypop <- YpopTruePS(ygroup_true, alpha, use = 'pairwise.complete.obs')
ypop_true <- ypop$ypop
ypop_true_var <- ypop$ypop_var


#  -------------------------------------------
#  ---  FOR THE ESTIMATED PROPENSITY SCORE ---
#  -------------------------------------------

# Calculating the ipw.

phi_hat_est <- list(coefs = summary(glmod)$coef[, 1],
                    re_var = as.numeric(summary(glmod)$varcor))
ygroup_est <- GroupIPW(dta = sim_dta, cov_cols = cov_cols, phi_hat = phi_hat_est,
                       gamma_numer = trt_coef, alpha = alpha, neigh_ind = neigh_ind,
                       keep_re_alpha = TRUE)$yhat_group

ypop <- YpopTruePS(ygroup_est, alpha, use = 'pairwise.complete.obs')
ypop_est <- ypop$ypop
ypop_est_var <- VarEstPS(sim_dta, ygroup_est, ypop_est, neigh_ind, phi_hat_est,
                         cov_cols, ypop$ypop_var)


par(mfrow = c(1, 2))
plot(ypop_true, ypop_est)
abline(a = 0, b = 1)

plot(ypop_true_var, ypop_est_var)
abline(a = 0, b = 1)
