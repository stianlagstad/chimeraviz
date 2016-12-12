context("Test scaling a list of integers to an interval [a, b]")

list123 <- c(1, 2, 3)
list123ScaledTo3 <- .scaleListToInterval(list123, min(list123), max(list123))
list123ScaledTo10 <- .scaleListToInterval(list123, 1, 10)

test_that("scaleListToInterval works", {
  expect_equal(list123ScaledTo3, list123)
  expect_equal(min(list123ScaledTo10), 1)
  expect_equal(max(list123ScaledTo10), 10)

  expect_error(.scaleListToInterval(c(1), 1, 10))
})
