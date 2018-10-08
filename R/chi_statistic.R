chi_statistic <- function(prob){
  # calculate
  chi_statistic <- sum((prob$expected_prob - prob$real_prob)^2 /
  prob$expected_prob)

  return (chi_statistic)
}
