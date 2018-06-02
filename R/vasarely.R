#' vasarely chart
#' @description vasarely produces a so called vasarely chart. E.g. for some (genetic) input data expected and real probabilty of (allel) frequencies can be calculated and plotted. For all possible combinations (of allels) the expected probabilty is described by a geom_raster plot. Each cell of the plot has a color which depends on the value of the expected probability. The real probability is described by a dot in the corresponding cell which also has a color that depends on the calculated probability. The chart and its colors for the two probabilities helps to check how expected and real probability differ from each other.
#' @param dat The input data to create the chart, e.g. genetic snp data. The input data must have exactly to columns.
#' @param colour The optional colors which should be used for the chart, e.g. blues9. Colour must be a character vector with at least two colors. If colour is NULL grey values will be taken.
#' @param name_xaxis The optional title for the x-axis of the plot.
#' @param name_yaxis The optional title for the y-axix of the plot
#'
#' @return returns the calculated vasarely chart of the input data
#' @export ??
#'
#' @examples
#' # create data
#' a1 <- c(rep("A", each = 25), rep("B", each = 75))
#' a2 <- c(rep("A", each = 50), rep("B", each = 50))
#' data <- data.frame(a1, a2)
#'
#' # use function
#' vasarely(dat = data, colour = c("yellow", "red"), name_xaxis = "a1", name_yaxis = "a2")


vasarely <- function(dat, colour = NULL, name_xaxis = NULL, name_yaxis = NULL){

  ## save input
  data <- as.data.frame(dat)
  color <- colour
  name_x <- name_xaxis
  name_y <- name_yaxis

  ## check input data
  # check if input data have to columns
  if(ncol(data) != 2){
    print("Input data must have exactly two columns!")
    return()
  # check if color vector contains character
  }else if(!is.character(color) && !is.null(color)){
    print("Colour data must be a character vector!")
    return()
  # check if color vector has more than one element
  }else if(length(color) < 2 && !is.null(color)){
    print("Colour vector must have more than one element!")
    return()
  # check if names of x- and y-axix are characters
  }else if((!is.null(name_x) && !is.character(name_x)) || (!is.null(name_y) && !is.character(name_y)) ){
    print("Name_xaxis and name_yaxis must be characters!")
    return()
  }


  # name columns to work with
  colnames(data) <- c("allel1", "allel2")

  ## compute a priori probability
  # get number of each allel in our two allels
  num_allel <- t(table(data$allel1) + table(data$allel2))
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
  prob_ex$allel1 <- NULL
  prob_ex$allel2 <- NULL




  ## compute real probability
  # new column with allel combination
  data$allel_combination <- paste(data$allel1,data$allel2)
  #compute
  real_probability <- as.data.frame(table(data$allel_combination)/nrow(data))
  colnames(real_probability) <- c("allel_comb", "real_prob")
  # check if there are missing allel combinations in real probability compared to expected, set values to 0
  prob <- merge(x = prob_ex, y = real_probability, by.x = "allel_comb", by.y = "allel_comb", all = TRUE)
  prob[is.na(prob)] <- 0


  # vectors needed for plottig later
  prob_real <- prob$real_prob
  prob_expected <- prob$expected_prob
  allel1 <- substring(prob$allel_comb,1,1)
  allel2 <- substring(prob$allel_comb,2)

  ## create plot for expected and real probability with geom_reaster and geom_dotplot
  # libraries needed for plotting
  library("ggplot2")
  library("forcats")


  # create plot and turn arount y-axis, so the origin of the plot is top left
  # changed real_probability to prob_real
  p <- ggplot(prob, aes(x = allel1, y = forcats::fct_rev(allel2))) +

          # create plot with squares for the expected probability,
          geom_raster(aes(fill = prob_expected), hjust = 0.5, vjust = 0.5) +

          # create plot with dots for the real probability,
          # binwidth: regulate size of dots, binaxis: direction to group dots
          # stackdir: dots in the center of squares, color: regulate the lines around the dots
          geom_dotplot(aes(fill = prob_real), binwidth = 0.90, binaxis = "y", stackdir = 'center', color = 0.001) +

          # title of legend, title of y-axis
          labs(fill = "probability", y = "allel2") +

          # put x-axis to the top of the plot
          scale_x_discrete(position = "top", expand = c(0,0)) +
          scale_y_discrete(expand = c(0,0)) +

          # fix size of legend text and size of legend title, legend title in bold
          theme(legend.text = element_text(size = 8),legend.title = element_text(size = 10, face = "bold"))



 # if no parameter for color take grey values
 if (is.null(color)){
   # grey values: probability = 1 is black, prob = 0 is white
   color <- grey.colors(256, start = 1, end = 0)
   # spread colors from probability 0 to 1
   p <- p + scale_fill_gradientn(colours = color, limits = c(0,1))
 }else{
 # else take color and spread it from probability 0 to 1
   p <- p + scale_fill_gradientn(colours = color, limits = c(0,1))
 }

 # if new name for x-axis is given in parameters
 if(!is.null(name_x)){
   p <- p + labs(x = name_x)
 }

 # if new name for y-axis is given in parameters
 if(!is.null(name_y)){
   p <- p + labs(y = name_y)
 }

  p
 return(p)
}

##################### do some tests######################
# load testdata (without column X) and test plot function
setwd("C:/Users/limes/Documents/Semester6/Bachelorarbeit/vasarely/R")
#testdata <- as.data.frame(read.csv("testdata5"))

## test grey_value parameter
#vasarely(testdata, testdata)

## test data with more than two columns
#testdata3 <- as.data.frame(read.csv("testdata"))
#vasarely(testdata3,blues9)

## test data with less than two columns
#testdata4 <- as.data.frame(read.csv("testdata_new"))
#testdata4$allel2 <- NULL
#vasarely(testdata4, blues9)

# test chart with grey values
#vasarely(testdata,blues9)
#vasarely(testdata)
#p + ggtitle("vasarely")

## example with numbers--> does not work with ten numbers but with 9????
a <- as.character(rep(c(1:10), each = 10))
b <- as.character(sample(c(1:10), 100, replace = TRUE))
data2 <- data.frame(a, b)

vasarely(data2)
