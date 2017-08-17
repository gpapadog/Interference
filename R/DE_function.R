DE <- function(ypop, ypop_var, alpha) {
  
  de <- array(NA, dim = c(2, length(alpha)))
  dimnames(de) <- list(stat = c('est', 'var'), alpha = alpha)

  de[1, ] <- ypop[1, ] - ypop[2, ]
  de[2, ] <- apply(ypop_var, 3, delta_method)
  
  return(de)
}