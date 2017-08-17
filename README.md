# Interference
## IPW estimators of causal effect estimands in the presence of interference

The estimands  of Hudgens and Halloran describe average potential outcomes (in groups or population)
under counterfactual treatment allocation that have a two-stage randomization design. In those
counterfactual treatment allocation, the direct effect describes changes in the average potential
outcome under treatment and control when treatment in the remaining of the neighbors is assigned
as independent Bernoulli trials with some fixed probability $\alpha$. Hudgens & Halloran provided
estimators in the situation of randomized trials.

Tchetgen Tchegen & Vanderweele provided IPW estimators of the estimands described above in
observational studies.

We define new estimands that have a practical interpretation in the context of air pollution and
other situations. In these estimands the counterfactual treatment allocation considered take the
covariates into consideration. Units within a cluster with a specific covariate set are more likely
to accept treatment compared to others. The new estimands take this information into consideration.

New IPW estimators are derived that are consistent and asymptotically normal under some regularity
conditions. The code in this repository calculates the IPW estimators for an assumed treatment
assignment that follows a logistic model with a random cluster intercept.
