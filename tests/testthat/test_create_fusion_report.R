context("Test creating Fusion reports")

test_that("create_fusion_report() produces a html file", {
  # Load data
  defuse833ke <- system.file(
    "extdata",
    "defuse_833ke_results.filtered.tsv",
    package = "chimeraviz")
  fusions <- import_defuse(defuse833ke, "hg19", 5)
  # Temporary file to store the report
  random_filename <- paste0(
    paste0(sample(LETTERS, 5, replace = TRUE), collapse=''),
    ".png"
  )
  create_fusion_report(fusions, random_filename)
  # See if we actually produced a file with size > 0
  expect_equal(file.info(random_filename)$size > 0, TRUE)
  # Delete the temporary file
  file.remove(random_filename)
})

test_that("create_fusion_report() fails with invalid input", {
  # List without fusion objects?
  expect_error(create_fusion_report(list(42)))
})
