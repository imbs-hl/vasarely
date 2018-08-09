#' plot_vasarely
#' @description plot_vasarely produces the plot of a vasarely chart and allows to add further functions of ggplot2.
#' @param vasarely_chart vasarely chart which has been calculated by the function vasarely
#' @return returns the plotted vasarely chart
#' @examples
#' # create data
#' a1 <- c(rep("A", each = 25), rep("B", each = 75))
#' a2 <- c(rep("A", each = 50), rep("B", each = 50))
#' data <- data.frame(a1, a2)
#'
#' # use function vasarely to calculate list with chart
#' vasarely <- vasarely(data = data)
#'
#' # extract plot from list
#' vasarely_plot <- vasarely$plot
#'
#' # put plot into plotting function and add further function of ggplot2
#' plot_vasarely(vasarely_plot) + labs(title = "My vasarely chart")

plot_vasarely <- function(vasarely_chart){
  vasarely_chart
}
