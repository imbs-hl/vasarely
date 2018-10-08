check_prob <- function(prob, lower_color_value, upper_color_value){
  # find minimum and maximum prob
  minimum <- min(min(prob$real_prob), min(prob$expected_prob))
  maximum <- max(max(prob$real_prob), max(prob$expected_prob))
  # check if corresponding to color values
  if((minimum > lower_color_value &&
      minimum > upper_color_value) ||
     (maximum < lower_color_value &&
      maximum < upper_color_value)){
    message("Chosen limits for color_values do not correspond to calculated relative frequencies!")
    message(paste0("Your minimum relative frequency is: ", minimum))
    message(paste0("Your maximum relative frequency is: ", maximum))
    return()
  }
}