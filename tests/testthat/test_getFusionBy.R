context("Test filter functions that operate on lists of Fusion objects")

defuse833ke <- system.file(
  "extdata",
  "defuse_833ke_results.filtered.tsv",
  package="chimeraviz")
fusions <- importDefuse(defuse833ke, "hg19", 3)

test_that("getFusionById works", {
  expect_equal(getFusionById(fusions, "58")@id, "58")
  expect_equal(getFusionById(fusions, "5267")@id, "5267")

  expect_error(getFusionById(fusions, 9999999))
})

test_that("getFusionByGeneName works", {
  expect_equal(length(getFusionByGeneName(fusions, "NEMF")), 1)
  expect_equal(length(getFusionByGeneName(fusions, "HENMT1")), 1)

  expect_warning(getFusionByGeneName(fusions, "Non-existent-gene-name"))
})

test_that("getFusionByChromosome works", {
  expect_equal(length(getFusionByChromosome(fusions, "chr11")), 1)
  expect_equal(length(getFusionByChromosome(fusions, "chr1")), 1)

  expect_warning(getFusionByChromosome(fusions, "Non-existent-chromosome"))
})
