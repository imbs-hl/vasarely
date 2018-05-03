#' vasarely chart
#'
#' @param dat the input data. For example, a matrix of genetic snp data
#'
#' @return returns the calculated vasarely chart of the input data
#' @export
#'
#' @examples

### do not forget: STR + L and rm(list=ls())

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
#compute
real_probability <- as.data.frame(table(testdata2$allel_combination)/nrow(testdata2))
colnames(real_probability) <- c("allel_comb", "real_prob")
real_probability

## put real probability and expected probability together into a dataframe
probabilities <- real_probability
probabilities$expected_prob <- as.vector(t(expected_probability))
probabilities$allel1 <- substring(probabilities$allel_comb,1,1)
probabilities$allel2 <- substring(probabilities$allel_comb,2)
probabilities

## create chart
library("ggplot2")
prob_ex <- probabilities$expected_prob
prob_real <- probabilities$real_prob
allel1 <- probabilities$allel1
allel2 <- probabilities$allel2
# create plot for expected probability with geom_reaster
l <- ggplot(probabilities,aes(allel1, allel2))
p <- l + geom_raster(aes(fill=prob_ex), hjust = 0.5, vjust = 0.5, interpolate = FALSE)
#print(p)
q <- p + geom_dotplot(aes(fill=prob_real),binwidth = 0.5,binaxis = "y")
print(q)

vasarely <- function(dat){

}
