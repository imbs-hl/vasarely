#' vasarely chart
#' @description vasarely produces a so called vasarely chart.
#' E.g. for some (genetic) input data expected and real
#' probabilty of (allel) frequencies can be calculated
#' and plotted. For all possible combinations (of allels)
#' the expected probability is described by a geom_raster
#' plot. Each cell of the plot has a color which depends
#' on the value of the expected probability. The real
#' probability is described by a dot in the corresponding
#' cell which also has a color that depends on the
#' calculated probability. The chart and its colors for
#' the two probabilities helps to check how expected
#' and real probability differ from each other.
#' @param data The input data to create the chart,
#' e.g. genetic snp data. The input data must have
#' exactly to columns.
#' @param color The optional colors which should be used
#' for the chart, e.g. blues9. Color must be a character
#' vector with at least two colors. If color is NULL grey
#' values will be taken.
#' @param name_xaxis The optional title for the x-axis of
#' the plot. Otherwise the x-axis is called "allele 2"
#' as a genetic input is expected.
#' @param name_yaxis The optional title for the y-axis of
#' the plot. Otherwise the y-axis is called "allele 1"
#' as a genetic input is expected.
#' @param lower_color_value The lower limit for
#' spreading the color over the probability values.
#' It must be a number between 0 and 1.
#' @param upper_color_value The upper limit for
#' spreading the color over the probability values.
#' It must be a number between 0 and 1.
#' @export
#' @import ggplot2 assertive forcats
#' @return returns a list of the vasarely chart and
#' a list of the p-value and statistic of
#' a chi-squared test.
#'
#' @examples
#' # create data
#' a1 <- c(rep("A", each = 25), rep("B", each = 75))
#' a2 <- c(rep("A", each = 50), rep("B", each = 50))
#' data <- data.frame(a1, a2)
#'
#' # use function
#' vasarely(data = data, color = c("yellow", "red"),
#' name_xaxis = "a1", name_yaxis = "a2")
#'


vasarely <- function(data,
                     color = grey.colors(256, start = 1,
                                         end = 0),
                     name_xaxis = "allele 2",
                     name_yaxis = "allele 1",
                     lower_color_value = 0,
                     upper_color_value = 1){

  ## save input
  data <- as.data.frame(data)

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



  # name columns to work with
  colnames(data) <- c("allel1", "allel2")


  ## compute a priori probability
  # get number of each allel in our two allels
  num_allel <- t(table(c(as.vector(data$allel1),
                         as.vector(data$allel2))))
  # get a total number of all allels
  num_all_allels <- ncol(data) * nrow(data)
  # compute
  a_priori_prob <- as.vector(num_allel / num_all_allels)



  ## compute expected probability
  expected_prob <- as.vector(a_priori_prob %*% t(a_priori_prob))
  prob_ex <- as.data.frame(expected_prob)
  # add possible allel combination to a new column
  col <- colnames(num_allel)
  n <- length(col)
  prob_ex$allel1 <- rep(col, each = n)
  prob_ex$allel2 <- col
  prob_ex$allel_comb <- paste(prob_ex$allel1, prob_ex$allel2)



  ## compute real probability
  # new column with allel combination
  data$allel_combination <- paste(data$allel1,data$allel2)
  #compute
  real_probability <- as.data.frame(table(data$allel_combination) /
                                      nrow(data))
  colnames(real_probability) <- c("allel_comb", "real_prob")

  # check if there are missing allel combinations in real probability
  # compared to expected, set values to 0
  prob <- merge(x = prob_ex,
                y = real_probability,
                by.x = "allel_comb",
                by.y = "allel_comb",
                all = TRUE)
  prob[is.na(prob)] <- 0


 # check if probability values correspond to the chosen limits for
 # spreading the colors
  min_real <- min(prob$real_prob)
  min_expected <- min(prob$expected_prob)
  max_real <- max(prob$real_prob)
  max_expected <- max(prob$expected_prob)
  minimum <- min(min_real, min_expected)
  maximum <- max(max_real, max_expected)

  if((minimum > lower_color_value &&
      minimum > upper_color_value) ||
     (maximum < lower_color_value &&
      maximum < upper_color_value)){
    message("Chosen limits for color_values do not correspond to calculated probabilities!")
    message(paste0("Your minimum probability is: ", minimum))
    message(paste0("Your maximum probability is: ", maximum))
    return()
  }

  ## create plot for expected and real probability
  ## with geom_reaster and geom_dotplot

  # libraries needed for plotting
  library("ggplot2")
  library("forcats")


  # create plot and turn arount y-axis, so the origin
  # of the plot is top left
  # changed real_probability to prob_real
  p <- ggplot(prob, aes(x = prob_ex$allel2,
                        y = forcats::fct_rev(prob_ex$allel1))) +

          # create plot with squares for the expected probability,
          geom_raster(aes(fill = prob$expected_prob),
                      hjust = 0.5,
                      vjust = 0.5) +

          # create plot with dots for the real probability,
          # binwidth: regulate size of dots,
          # binaxis: direction to group dots
          # stackdir: dots in the center of squares#
          # color: regulate the lines around the dots
          geom_dotplot(aes(fill = prob$real_prob),
                       binwidth = 0.90,
                       binaxis = "y",
                       stackdir = 'center',
                       color = 0.001) +
          # take color and spread it from lower to upper color value
          scale_fill_gradientn(colours = color,
                               limits = c(lower_color_value,
                                          upper_color_value)) +

          # title of legend, title of y-axis
          labs(fill = "probability",
               x = name_xaxis,
               y = name_yaxis) +

          # put x-axis to the top of the plot
          scale_x_discrete(position = "top",
                           expand = c(0,0)) +
          scale_y_discrete(expand = c(0,0)) +
          # x-axis and y-axis should always have the same length
          # so the fields are squares
          coord_fixed(ratio = 1/1) +

          # fix size of legend text and size of legend title
          # legend title in bold
          theme(legend.text = element_text(size = 8),
                legend.title = element_text(size = 10,
                                            face = "bold")) +
          # put title of the legend to the top, so it is readable
          guides(fill = guide_legend(title.position = "top"))


  ## calculate Chi-Squared-Test and put result into plot:
  chi <- chisq.test(x = data$allel1,
                    y = data$allel2,
                    correct = FALSE)
  p <- p + labs(caption = paste0( "Chi-squared test: p-value ",
                                  round(x = chi$p.value, digits = 8),
                                  ", statistic ",
                                  round(x=chi$statistic, digits = 4)))
  # create list with calculated data and plot
  list <- list(statistic_values = list(p_value = chi$p.value,
                                       statistic = chi$statistic),
               plot = p)

  ## call plotting function and return list
  plot_vasarely(p)

  return(list)
}
