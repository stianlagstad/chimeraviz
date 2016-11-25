context("Test plotCircle()")

# Load data
defuse833ke <- system.file(
  "extdata",
  "defuse_833ke_results.filtered.tsv",
  package="chimeraviz")
fusions <- importDefuse(defuse833ke, "hg19", 3)
# Temporary file to store the plot
pngFilename <- tempfile(
  pattern = "circlePlot",
  fileext = ".png",
  tmpdir = tempdir())

test_that("plotCircle produces a png file", {
  # Open device
  png(pngFilename, width = 1000, height = 750)
  # Plot and expect certain message output
  plotCircle(fusions)
  # Close device
  dev.off()
  # See if we actually produced a file with size > 0
  expect_equal(file.info(pngFilename)$size > 0, TRUE)
})

test_that("plotCircle() fails with invalid input", {
  # List without fusion objects?
  expect_error(plotCircle(list()))
  # Inputs that are not lists?
  expect_error(plotCircle(42))
  expect_error(plotCircle("hey this should fail"))
})
