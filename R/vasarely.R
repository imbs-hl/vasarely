#' vasarely chart
#'
#' @param dat the input data. For example, a matrix of genetic snp data
#'
#' @return returns the calculated vasarely chart of the input data
#' @export
#'
#' @examples

### do not forget: STR + L and rm(list=ls())


vasarely <- function(dat){
  # prepare data
  dat <- dat
  colnames(dat) <- c("allel1", "allel2")

  ## compute a priori probability
  # get number of each allel in our two allels
  num_allel <- t(table(dat$allel1) + table(dat$allel2))
  # get a total number of all allels
  num_all_allels <- ncol(dat) * nrow(testdata)
  # compute
  a_priori_prob <- as.vector(num_allel / num_all_allels)


  ## compute expected probability
  expected_probability <- as.data.frame(a_priori_prob %*% t(a_priori_prob))
  colnames(expected_probability) <- colnames(num_allel)
  rownames(expected_probability) <- colnames(num_allel)


  ## compute real probability
  # new column with allel combination
  dat$allel_combination <- paste(testdata$allel1,testdata$allel2)
  #compute
  real_probability <- as.data.frame(table(dat$allel_combination)/nrow(dat))
  colnames(real_probability) <- c("allel_comb", "real_prob")

  ## put real probability and expected probability together into a dataframe
  probabilities <- real_probability
  probabilities$expected_prob <- as.vector(t(expected_probability))
  probabilities$allel1 <- substring(probabilities$allel_comb,1,1)
  probabilities$allel2 <- substring(probabilities$allel_comb,2)

  ## prepare chart
  library("ggplot2")
  prob_ex <- probabilities$expected_prob
  prob_real <- probabilities$real_prob
  allel1 <- probabilities$allel1
  allel2 <- probabilities$allel2

  # create plot for expected and real probability with geom_reaster and geom_dotplot
  print(ggplot(probabilities,aes(allel1, allel2)) +
          geom_raster(aes(fill=prob_ex), hjust = 0.5, vjust = 0.5, interpolate = FALSE) +
          geom_dotplot(aes(fill=prob_real), binwidth = 0.8, binaxis = "y", stackdir='center', color = 0.01))
}

# load second testdata (without column X) and test plot function
data2 <- testdata <- as.data.frame(read.csv("testdata2"))
vasarely(data2)
