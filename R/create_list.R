create_list <- function(prob, data, 
                        chi_p_value, chi_statistic, 
                        num_observation, p){
  # delete columns in dataframe prob,
  # they are not needed for returning values
  prob$allel1 <- NULL
  prob$allel2 <- NULL
  # create list with calculated data and plot
  colnames(prob) <- c("allel_comb", "expected_freq", "real_freq")
  list <- list(statistic_values = list(p_value = chi_p_value,
                                       statistic = chi_statistic),
               rel_frequencies = prob,
               number_observations = nrow(data),
               plot = p)
  return(list)
  }