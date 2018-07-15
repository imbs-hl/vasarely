## load different testdata
testdata <- as.data.frame(read.csv("testdata5"))
#testdata with more than one column
testdata_more_columns <- as.data.frame(read.csv("testdata"))
# testdata with less than two columns
testdata_less_columns <- as.data.frame(read.csv("testdata_new"))
testdata_less_columns$allel2 <- NULL
#example
a1 <- c(rep("A", each = 25), rep("B", each = 75))
a2 <- c(rep("A", each = 50), rep("B", each = 50))
example <- data.frame(a1, a2)
#testdata with numbers
a <- rep(c(1:10), each = 10)
b <- sample(c(1:10), 100, replace = TRUE)
numbers <- data.frame(a, b)
#empty testdata
empty_data <- testdata
empty_data$allel1 <- NULL
empty_data$allel2 <- NULL

context("package")
test_that("package is successfully installed", {
  expect_true("vasarely" %in% rownames(installed.packages()))
})

context("ncol")
test_that("testdata have x columns", {
  expect_equal(ncol(testdata), 2)
  expect_equal(ncol(testdata_less_columns), 1)
  expect_equal(ncol(testdata_more_columns), 3)
  expect_equal(ncol(empty_data), 0)
})

context("columns")
test_that("testdata have two columns", {
  expect_message(vasarely(empty_data), "Input data must have exactly two columns!", fixed = TRUE)
  expect_message(vasarely(testdata_less_columns), "Input data must have exactly two columns!", fixed = TRUE)
  expect_message(vasarely(testdata_more_columns), "Input data must have exactly two columns!", fixed = TRUE)
})

context("color")
test_that("color parameter input", {
  expect_message(vasarely(testdata, testdata), "Color data must be a character vector!")
  expect_message(vasarely(testdata, 0), "Color data must be a character vector!")
  expect_message(vasarely(testdata, c(3:6)), "Color data must be a character vector!")
  expect_message((vasarely(testdata, "a")), "Color vector must have more than one element!")
  expect_message(vasarely(testdata, c("yellow")), "Color vector must have more than one element!")
  expect_message(vasarely(testdata, list("yellow", "red")), "Color data must be a character vector!")
  })

context("greyvalueparameter")
test_that("upper and lower grey value parameters work", {
  # value is not a number
  expect_message(vasarely(example, upper_color_value = a), "upper_color_value must be a number!", fixed = TRUE)
  expect_message(vasarely(example, upper_color_value = "2"), "upper_color_value must be a number!", fixed = TRUE)
  expect_message(vasarely(example, lower_color_value = a), "lower_color_value must be a number!", fixed = TRUE)
  expect_message(vasarely(example, lower_color_value = "2"), "lower_color_value must be a number!", fixed = TRUE)
  # value is a vector
  expect_message(vasarely(example, upper_color_value = c(0:1)), "upper_color_value must be a number!", fixed = TRUE)
  expect_message(vasarely(example, lower_color_value = c(0:1)), "lower_color_value must be a number!", fixed = TRUE)
  # value > 1
  expect_message(vasarely(example, lower_color_value = 1.1), "lower_color_value must be between 0 and 1!", fixed = TRUE)
  expect_message(vasarely(example, upper_color_value = 1.1), "upper_color_value must be between 0 and 1!", fixed = TRUE)
  # value < 0
  expect_message(vasarely(example, lower_color_value = -0.1), "lower_color_value must be between 0 and 1!", fixed = TRUE)
  expect_message(vasarely(example, upper_color_value = -0.1), "upper_color_value must be between 0 and 1!", fixed = TRUE)
  # lower value is bigger than upper
  expect_message(vasarely(example, lower_color_value = 1, upper_color_value = 0), "lower_color_value must be smaller than upper_color_ value!", fixed = TRUE)
  # lower color value does not correspond with prob. values for testdata:
  # min expected prob: 0.140625 , min real prob: 0.25,
  # max expected prob: 0.390625 , max real prob: 0.5
  expect_message(vasarely(example, lower_color_value = 0.140625, upper_color_value = 0.4), "Chosen limits for color_values do not correspond to calculated probabilities!")
  expect_message(vasarely(example, lower_color_value = 0.15, upper_color_value = 0.6), "Chosen limits for color_values do not correspond to calculated probabilities!")
  # upper color value does not correspond with prob. values
  expect_message(vasarely(example, lower_color_value = 0.1, upper_color_value = 0.3), "Chosen limits for color_values do not correspond to calculated probabilities!")
  expect_message(vasarely(example, lower_color_value = 0.1, upper_color_value = 0.45), "Chosen limits for color_values do not correspond to calculated probabilities!")
  })

#devtools::test()
