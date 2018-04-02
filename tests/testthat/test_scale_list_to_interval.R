context("Test scaling a list of integers to an interval [a, b]")

list123 <- c(1, 2, 3)
list_123_scaled_to_3 <- .scale_list_to_interval(
  list123,
  min(list123),
  max(list123)
)
list_123_scaled_to_10 <- .scale_list_to_interval(list123, 1, 10)

test_that(".scale_list_to_interval works", {
  expect_equal(list_123_scaled_to_3, list123)
  expect_equal(min(list_123_scaled_to_10), 1)
  expect_equal(max(list_123_scaled_to_10), 10)

  expect_error(.scale_list_to_interval(c(1), 1, 10))
})
