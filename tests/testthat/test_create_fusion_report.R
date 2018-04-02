context("Test creating Fusion reports")

# Load data
defuse833ke <- system.file(
  "extdata",
  "defuse_833ke_results.filtered.tsv",
  package = "chimeraviz")
fusions <- import_defuse(defuse833ke, "hg19", 5)
# Temporary file to store the report
output_filename <- tempfile(
  pattern = "fusionReport",
  fileext = ".html",
  tmpdir = tempdir())

test_that("create_fusion_report() produces a html file", {
  create_fusion_report(fusions, output_filename)
  # See if we actually produced a file with size > 0
  expect_equal(file.info(output_filename)$size > 0, TRUE)
})

test_that("create_fusion_report() fails with invalid input", {
  # List without fusion objects?
  expect_error(create_fusion_report(list(42)))
})
