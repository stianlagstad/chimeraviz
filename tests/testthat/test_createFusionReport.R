context("Test creating Fusion reports")

# Load data
defuse833ke <- system.file(
  "extdata",
  "defuse_833ke_results.filtered.tsv",
  package="chimeraviz")
fusions <- importDefuse(defuse833ke, "hg19", 5)
# Temporary file to store the report
outputFilename <- tempfile(
  pattern = "fusionReport",
  fileext = ".html",
  tmpdir = tempdir())

test_that("createFusionReport() produces a html file", {
  createFusionReport(fusions, outputFilename)
  # See if we actually produced a file with size > 0
  expect_equal(file.info(outputFilename)$size > 0, TRUE)
})

test_that("createFusionReport() fails with invalid input", {
  # List without fusion objects?
  expect_error(createFusionReport(list(42)))
})
