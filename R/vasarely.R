#' vasarely chart
#'
#' @param dat the input data. For example, a matrix of genetic snp data
#'
#' @return returns the calculated vasarely chart of the input data
#' @export
#'
#' @examples

## do not forget: rm(list=ls())
## now just trying something that works with test data :)

# load our testdata
testdata <- as.data.frame(read.csv("testdata"))
num_allel <- t(table(testdata$allel1) + table(testdata$allel2))
names_colums <- colnames(num_allel)

# calculate a priori probability
num_all_allels <- (ncol(testdata) - 1) * nrow(testdata)
a_priori_prob <- as.vector(num_allel / num_all_allels)

# calculate expected probability
expected_probability <- as.data.frame(a_priori_prob %*% t(a_priori_prob))
colnames(expected_probability) <- names_colums
rownames(expected_probability) <- names_colums


# calculate real probability
testdata2 <- testdata
testdata2$allel_combination <- paste(testdata$allel1,testdata$allel2)
testdata2$allel1 <- NULL
testdata2$allel2 <- NULL
real_probability <- as.data.frame(table(testdata2$allel_combination)/nrow(testdata2))
colnames(real_probability) <- c("allel_comb", "prob")



vasarely <- function(dat){

}
