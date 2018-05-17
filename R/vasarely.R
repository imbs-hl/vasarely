#' vasarely chart
#'
#' @param dat the input data. For example, a matrix of genetic snp data
#'
#' @return returns the calculated vasarely chart of the input data
#' @export
#'
#' @examples

### do not forget: STR + L and rm(list=ls())


vasarely <- function(dat, colour = NULL, name_xaxix = NULL, name_yaxis = NULL){
  # save input
  data <- as.data.frame(dat)
  color <- colour
  name_x <- name_xaxix
  name_y <- name_yaxis

  # check input data
  if(ncol(data) != 2){
    print("Input data must have exactly two columns!")
    return()
  }else if(!is.character(color) && !is.null(color)){
    print("Colour data must be a character vector!")
    return()
  }else if(length(color) < 2 && !is.null(color)){
    print("Colour vector must have more than one element!")
    return()
  }else if((!is.null(name_x) && !is.character(name_x)) || (!is.null(name_y) && !is.character(name_y)) ){
    print("Name_xaxis and name_yaxis must be characters!")
    return()
  }


  #name columns to work with
  colnames(data) <- c("allel1", "allel2")

  ## compute a priori probability
  # get number of each allel in our two allels
  num_allel <- t(table(data$allel1) + table(data$allel2))
  # get a total number of all allels
  num_all_allels <- ncol(data) * nrow(data)
  # compute
  a_priori_prob <- as.vector(num_allel / num_all_allels)


  ## compute expected probability
  prob_ex <- as.vector(a_priori_prob %*% t(a_priori_prob))


  ## compute real probability
  # new column with allel combination
  data$allel_combination <- paste(data$allel1,data$allel2)
  #compute
  real_probability <- as.data.frame(table(data$allel_combination)/nrow(data))
  colnames(real_probability) <- c("allel_comb", "real_prob")
  # vectors needed for plottig later
  prob_real <- real_probability$real_prob
  allel1 <- substring(real_probability$allel_comb,1,1)
  allel2 <- substring(real_probability$allel_comb,2)


  ## create plot for expected and real probability with geom_reaster and geom_dotplot
  library("ggplot2")
 p <-  ggplot(real_probability,aes(allel1, allel2)) +
          geom_raster(aes(fill=prob_ex), hjust = 0.5, vjust = 0.5, interpolate = FALSE) +
          geom_dotplot(aes(fill=prob_real), binwidth = 0.98, binaxis = "y", stackdir='center', color = 0.01) +
          labs(fill="probability") +
          scale_x_discrete(position = "top", expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
          theme(legend.text=element_text(size = 8),legend.title=element_text(size = 10, face = "bold"))

 # if no parameter for color take grey values
 if (is.null(color)){
   color <- grey.colors(256, start = 1, end = 0)
   p <- p + scale_fill_gradientn(colours = color, limits = c(0,1))
 }else{
 # take color
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

# load testdata (without column X) and test plot function
setwd("C:/Users/limes/Documents/Semester6/Bachelorarbeit/vasarely/R")
testdata <- as.data.frame(read.csv("testdata5"))

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
vasarely(testdata)
#p + ggtitle("vasarely")

