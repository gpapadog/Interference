delta_method <- function(x, vec = c(1, - 1))  {
  vec <- matrix(vec, ncol = 1)
  return(t(vec) %*% x %*% vec)
}