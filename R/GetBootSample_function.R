GetBootSample <- function(dta) {
  
  num_clus <- max(dta$neigh)
  boot_clusters <- sample(1 : num_clus, num_clus, replace = TRUE)
  
  # Binding data without accidentally merging repeated clusters.
  boot_dta <- NULL
  for (nn in 1 : num_clus) {
    D <- subset(dta, neigh == boot_clusters[nn])
    D$neigh <- nn
    boot_dta <- rbind(boot_dta, D)
  }
  
  return(list(boot_dta = boot_dta, chosen_clusters = boot_clusters))
}