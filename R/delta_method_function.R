delta_method <- function(x, vec = c(1, - 1))  {
  
  if (length(vec) != dim(x)[1] | length(vec) != dim(x)[2]) {
    stop('Wrong dimensions.')
  }
  
  vec <- matrix(vec, ncol = 1)
  return(t(vec) %*% x %*% vec)
}