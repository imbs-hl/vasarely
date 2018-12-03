observed_frequency <- function(data){
  # new column with allel combination
  data$allel_combination <- paste(data$allel1,data$allel2)
  # compute
  real_probability <- as.data.frame(table(data$allel_combination) /
                                      nrow(data))
  colnames(real_probability) <- c("allel_comb", "real_prob")

  return(real_probability)
}
