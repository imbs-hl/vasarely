#' vasarely chart
#' @description vasarely produces a so called vasarely chart.
#' E.g. for some (genetic) input data expected and observed
#' relative frequencies of allele combinations can be calculated
#' and plotted. For all possible combinations
#' the expected relative frequency is described by a geom_raster
#' plot. Each cell of the plot has a color which depends
#' on the value of the expected relative frequcency. The real
#' relative frequency is described by a dot in the corresponding
#' cell which also has a color that depends on the
#' calculated relative frequency. The chart and its colors for
#' the two relative frequencies help to check how expected
#' and observed relative frequencies differ from each other.
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
#' spreading the color over the relative frequencies values.
#' It must be a number between 0 and 1.
#' @param upper_color_value The upper limit for
#' spreading the color over the relative frequencies values.
#' It must be a number between 0 and 1.
#' @export
#' @import ggplot2 assertive forcats
#' @return returns a list of the calculated relative
#' frequencies, the number of observations, the vasarely
#' chart and a list of the p-value and statistic of
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
                     color = grey.colors(256,
                                         start = 1,
                                         end = 0),
                     name_xaxis = "allele 2",
                     name_yaxis = "allele 1",
                     lower_color_value = 0,
                     upper_color_value = 1){

  # save input and prepare data
  data <- as.data.frame(data)
  colnames(data) <- c("allel1", "allel2")
  num_observation <- nrow(data)

  # check input parameters
  check_parameters(data, color,
                   name_xaxis, name_yaxis,
                   lower_color_value,  upper_color_value)

  # compute a priori probability
  a_priori_prob <- a_priori(data)

  # compute expected relative frequency
  prob_ex <- exp_frequency(a_priori_prob, data)

  # compute observed relative frequency
  real_probability <- observed_frequency(data)

  # check if there are missing allel combinations in observed
  # relative frequencies compared to expected, set values to 0
  prob <- merge(x = prob_ex,
                y = real_probability,
                by.x = "allel_comb",
                by.y = "allel_comb",
                all = TRUE)
  prob[is.na(prob)] <- 0

  # add observed relative frequencies of heterozygotes
  # so they have the same prob.
  prob$real_prob <- obs_freq_heterozygotes(prob)

  # check if relative frequency values correspond to the
  # chosen limits for spreading the colors
  check_prob(prob, lower_color_value, upper_color_value)

  # calculate chi-squared-test
  chi_statistic <- chi_statistic(prob) * num_observation
  chi_p_value <- pchisq(chi_statistic,
                        df = (sqrt(nrow(prob)) - 1) ^ 2,
                        lower.tail = FALSE)

  # create plot for expected and observed relative frequencies
  p <- plot_vasarely(prob, prob_ex,
                     color, upper_color_value, lower_color_value,
                     name_xaxis, name_yaxis,
                     chi_statistic, chi_p_value)
  p

  # create list with calculated values and return it
  list <- create_list(prob, data,
                      chi_p_value, chi_statistic,
                      num_observation, p)

  return(list)
}
