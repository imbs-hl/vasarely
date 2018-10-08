obs_freq_heterozygotes <- function(prob){
  # two copies of prob, one is ordered by allel1 the other by allel2
  prob1 <- prob[order(prob$allel1),]
  prob2 <- prob[order(prob$allel2),]
  # add prob1 and prob2 so heterozygotes have the same real relative frequency
  prob3 <- as.data.frame(prob1$real_prob + prob2$real_prob)
  colnames(prob3)[1] <- "real_prob"
  # real probs for homozygotes were doubled and have to be corrected
  # real probs for heterozygotes have to be split for both possible combinations
  prob3$real_prob <- prob3$real_prob / 2

  return(prob3$real_prob)
}
