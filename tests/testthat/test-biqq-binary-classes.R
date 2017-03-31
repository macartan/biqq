context("Binary class frequencies")

# XYK clue data
clue_data <-
  data.frame(X = c(rep(1, times = 12), rep(0, times = 8)),
             Y = c(rep(1,11), 0, rep(0:1, each = 4)),
             K = c(rep(1, times = 4), rep(0, times = 8), 1, rep(0, times = 7)))

test_that("creation of binary class frequencies works", {
  # for XYK data
  expect_equivalent(biqq_binary_classes(clue_data, ordering = c(1,5,3,7,2,6,4,8)), c(3,1,4,7,1,0,0,4))
  # for XY data
  expect_equivalent(biqq_binary_classes(clue_data[,c("X","Y")], ordering = c(1,3,2,4)), c(4,1,4,11))
})

test_that("errors in creation of binary class frequencies work", {
  # for duplicates in ordering
  expect_error(biqq_binary_classes(clue_data, ordering = c(1,5,3,7,2,6,8,8)))
  # for ordering indexes out of range
  expect_error(biqq_binary_classes(clue_data[,c("X","Y")], ordering = c(1,3,2,5)))
  # for length of ordering exceeding actual range length
  expect_error(biqq_binary_classes(clue_data[,c("X","Y")], ordering = c(1,3,2,4,5)))
})
