context("Test importing STAR-Fusion data into a list of Fusion objects")

starfusion_data <- system.file(
  "extdata",
  "star-fusion.fusion_predictions.abridged.annotated.coding_effect.tsv",
  package = "chimeraviz")
fusions <- import_starfusion(starfusion_data, "hg19", 3)
fusion <- fusions[[1]]

test_that("import_starfusion imports data into Fusion objects", {
  expect_equal(length(fusions), 3)
  expect_equal(fusion@id, "1")
  expect_equal(fusion@gene_upstream@name, "RP11-9L18.2")
  expect_equal(fusion@gene_upstream@chromosome, "chr1")
  expect_equal(fusion@gene_upstream@ensembl_id, "ENSG00000215835")
})

test_that("import_starfusion fails to import data", {
  expect_error(import_starfusion("Non-existent-file", "hg18"))
  expect_error(import_starfusion(starfusion_data, "Invalid-genome"))
  expect_error(import_starfusion(starfusion_data, "hg19", "invalid limit"))
})

test_that("import_starfusion limits work", {
  expect_equal(length(import_starfusion(starfusion_data, "hg19", 2)), 2)
})

test_that("import_starfusion imports additional columns when available", {
  expect_equal(fusions[[3]]@gene_upstream@junction_sequence@length, 1431)
  expect_true(!is.null(fusion@fusion_tool_specific_data$annots))
  expect_true(!fusions[[3]]@inframe)
})
