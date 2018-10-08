exp_frequency <- function(a_priori_prob, data){
  # calculate expected relative frequencies
  expected_prob <- as.vector(a_priori_prob %*% t(a_priori_prob))
  prob_ex <- as.data.frame(expected_prob)
  # add possible allele combination to a new column
  col <- colnames(t(table(c(as.vector(data$allel1),
                            as.vector(data$allel2)))))
  n <- length(col)
  prob_ex$allel1 <- rep(col, each = n)
  prob_ex$allel2 <- col
  prob_ex$allel_comb <- paste(prob_ex$allel1, prob_ex$allel2)
  
  return(prob_ex)
}