#' vasarely chart
#'
#' @param dat the input data. For example, a matrix of genetic snp data
#'
#' @return returns the calculated vasarely chart of the input data
#' @export
#'
#' @examples

### do not forget: rm(list=ls())

### now just trying something that works with test data :)

## load our testdata
testdata <- as.data.frame(read.csv("testdata"))

## compute a priori probability
# get number of each allel in our two allels
num_allel <- t(table(testdata$allel1) + table(testdata$allel2))
# get a total number of all allels, -1 because first column is the id of our data
num_all_allels <- (ncol(testdata) - 1) * nrow(testdata)
# compute
a_priori_prob <- as.vector(num_allel / num_all_allels)
a_priori_prob

## compute expected probability
expected_probability <- as.data.frame(a_priori_prob %*% t(a_priori_prob))
colnames(expected_probability) <- colnames(num_allel)
rownames(expected_probability) <- colnames(num_allel)
expected_probability


## compute real probability
# maybe we need variable testdata later..
testdata2 <- testdata
# new column with allel combination
testdata2$allel_combination <- paste(testdata$allel1,testdata$allel2)
# delete both colums with redundant data
testdata2$allel1 <- NULL
testdata2$allel2 <- NULL
#compute
real_probability <- as.data.frame(table(testdata2$allel_combination)/nrow(testdata2))
colnames(real_probability) <- c("allel_comb", "prob")
real_probability



vasarely <- function(dat){

}
