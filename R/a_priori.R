a_priori <- function(data){
  # get number of each allel in our two allels
  num_allel <- t(table(c(as.vector(data$allel1),
                         as.vector(data$allel2))))
  # get a total number of all allels
  num_all_allels <- ncol(data) * nrow(data)
  # compute
  a_priori_prob <- as.vector(num_allel / num_all_allels)
  
  return(a_priori_prob)
}