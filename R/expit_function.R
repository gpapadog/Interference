#' Expit.
#' 
#' Calculating the expit (inverse logit) of a number.
#' 
#' @export
expit <- function(x) {
  return(exp(x) / (1 + exp(x)))
}