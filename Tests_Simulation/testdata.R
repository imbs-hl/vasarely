# create a vector as basis to create allels
x = c("A", "B", "C", "D", "E")

# use vector to create two random allel vectors with 200 elements each
allel1 <- sample(x, 200, replace = TRUE)
allel2 <- sample(x, 200, replace = TRUE)

# put the two allel vectors into a dataframe
two_allels <- data.frame(allel1, allel2)

# save the dataframe, so we can always use the same data for testing
write.csv(two_allels, file = "testdata5", row.names = FALSE)


# modifiy first testdata for test of the vasarely function: delete first column
#testdata2 <- as.data.frame(read.csv("testdata"))
#testdata2$X <- NULL
#write.csv(testdata2, file = "testdata_new", row.names=F)
