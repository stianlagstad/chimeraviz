context("Test plot_circle()")

# Load data
defuse833ke <- system.file(
  "extdata",
  "defuse_833ke_results.filtered.tsv",
  package = "chimeraviz")
fusions <- import_defuse(defuse833ke, "hg19", 3)
# Temporary file to store the plot
png_filename <- tempfile(
  pattern = "circlePlot",
  fileext = ".png",
  tmpdir = tempdir())

test_that("plot_circle produces a png file", {
  # Open device
  png(png_filename, width = 1000, height = 750)
  # Plot and expect certain message output
  plot_circle(fusions)
  # Close device
  dev.off()
  # See if we actually produced a file with size > 0
  expect_equal(file.info(png_filename)$size > 0, TRUE)
})

test_that("plot_circle() fails with invalid input", {
  # List without fusion objects?
  expect_error(plot_circle(list()))
  # Inputs that are not lists?
  expect_error(plot_circle(42))
  expect_error(plot_circle("hey this should fail"))
})
