check_parameters <- function(data, color, 
                             name_xaxis, name_yaxis, 
                             lower_color_value,  upper_color_value){
  library(assertive)
  ## check input data
  # check if input data have to columns
  if(ncol(data) != 2){
    message("Input data must have exactly two columns!")
    return()
    # check if color vector contains character
  } else if(!is.character(color) && !is.null(color)){
    message("Color data must be a character vector!")
    return()
    # check if color vector has more than one element
  } else if(length(color) < 2 && !is.null(color)){
    message("Color vector must have more than one element!")
    return()
    # check if names of x- and y-axis are characters
  } else if((!is.character(name_xaxis)) ||
            !is.character(name_yaxis) ||
            length(name_xaxis) > 1 ||
            length(name_yaxis) > 1){
    message("Name_xaxis and name_yaxis must be characters!")
    return()
    # check if lower_color_value is a number
  } else if(!is_a_number(lower_color_value)){
    message("lower_color_value must be a number!")
    return()
    # check if lower_color_value is between 0 and 1
  } else if(lower_color_value > 1 || lower_color_value < 0){
    message("lower_color_value must be between 0 and 1!")
    return()
    # check if upper_color_value is a number
  } else if(!is_a_number(upper_color_value)){
    message("upper_color_value must be a number!")
    return()
    # check if upper_color_value is between 0 and 1
  } else if(upper_color_value > 1 || upper_color_value < 0){
    message("upper_color_value must be between 0 and 1!")
    return()
    # check if lower_color_value is smaller than upper_color_value
  } else if(!is.null(lower_color_value) &&
            !is.null(upper_color_value) &&
            lower_color_value >= upper_color_value){
    message("lower_color_value must be smaller than upper_color_ value!")
    return()
  }
}