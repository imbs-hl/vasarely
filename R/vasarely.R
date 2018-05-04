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
  # prepare data: save input in dataframe and name columns to work with
  dat <- dat
  colnames(dat) <- c("allel1", "allel2")

  ## compute a priori probability
  # get number of each allel in our two allels
  num_allel <- t(table(dat$allel1) + table(dat$allel2))
  # get a total number of all allels
  num_all_allels <- ncol(dat) * nrow(dat)
  # compute
  a_priori_prob <- as.vector(num_allel / num_all_allels)


  ## compute expected probability
  prob_ex <- as.vector(t(as.data.frame(a_priori_prob %*% t(a_priori_prob))))


  ## compute real probability
  # new column with allel combination
  dat$allel_combination <- paste(dat$allel1,dat$allel2)
  #compute
  real_probability <- as.data.frame(table(dat$allel_combination)/nrow(dat))
  colnames(real_probability) <- c("allel_comb", "real_prob")
  # vectors needed for plottig later
  prob_real <- real_probability$real_prob
  allel1 <- substring(real_probability$allel_comb,1,1)
  allel2 <- substring(real_probability$allel_comb,2)


  ## create plot for expected and real probability with geom_reaster and geom_dotplot
  library("ggplot2")
  print(ggplot(real_probability,aes(allel1, allel2)) +
          geom_raster(aes(fill=prob_ex), hjust = 0.5, vjust = 0.5, interpolate = FALSE) +
          geom_dotplot(aes(fill=prob_real), binwidth = 0.8, binaxis = "y", stackdir='center', color = 0.01))
}

# load second testdata (without column X) and test plot function
testdata2 <- as.data.frame(read.csv("testdata2"))
vasarely(testdata2)
