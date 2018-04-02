context("Test filter functions that operate on lists of Fusion objects")

defuse833ke <- system.file(
  "extdata",
  "defuse_833ke_results.filtered.tsv",
  package = "chimeraviz")
fusions <- import_defuse(defuse833ke, "hg19", 3)

test_that("get_fusion_by_id works", {
  expect_equal(get_fusion_by_id(fusions, "58")@id, "58")
  expect_equal(get_fusion_by_id(fusions, "5267")@id, "5267")

  expect_error(get_fusion_by_id(fusions, 9999999))
})

test_that("get_fusion_by_gene_name works", {
  expect_equal(length(get_fusion_by_gene_name(fusions, "NEMF")), 1)
  expect_equal(length(get_fusion_by_gene_name(fusions, "HENMT1")), 1)

  expect_warning(get_fusion_by_gene_name(fusions, "Non-existent-gene-name"))
})

test_that("get_fusion_by_chromosome works", {
  expect_equal(length(get_fusion_by_chromosome(fusions, "chr11")), 1)
  expect_equal(length(get_fusion_by_chromosome(fusions, "chr1")), 1)

  expect_warning(get_fusion_by_chromosome(fusions, "Non-existent-chromosome"))
})
