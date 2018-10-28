context("Test importing STAR-Fusion data into a list of Fusion objects")

starfusion_data <- system.file(
  "extdata",
  "star-fusion.fusion_candidates.final.abridged.txt",
  package = "chimeraviz")
fusions <- import_starfusion(starfusion_data, "hg19", 3)
fusion <- fusions[[1]]

test_that("import_starfusion imports data into Fusion objects", {
  expect_equal(length(fusions), 3)
  expect_equal(fusion@id, "1")
  expect_equal(fusion@gene_upstream@name, "THRA")
  expect_equal(fusion@gene_upstream@chromosome, "chr17")
  expect_equal(fusion@gene_upstream@ensembl_id, "ENSG00000126351")
})

test_that("import_starfusion fails to import data", {
  expect_error(import_starfusion("Non-existent-file", "hg18"))
  expect_error(import_starfusion(starfusion_data, "Invalid-genome"))
  expect_error(import_starfusion(starfusion_data, "hg19", "invalid limit"))
})

test_that("import_starfusion limits work", {
  expect_equal(length(import_starfusion(starfusion_data, "hg19", 2)), 2)
})
