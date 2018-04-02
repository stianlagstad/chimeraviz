context("Test importing FusionMap data into a list of Fusion objects")

fusionmap_data <- system.file(
  "extdata",
  "FusionMap_01_TestDataset_InputFastq.FusionReport.txt",
  package = "chimeraviz")
fusions <- import_fusionmap(fusionmap_data, "hg19", 3)
fusion <- fusions[[1]]

test_that("import_fusionmap imports data into Fusion objects", {
  expect_equal(length(fusions), 3)
  expect_equal(fusion@id, "FUS_11184689_446432191(++)")
  expect_equal(fusion@gene_upstream@name, "FADS3")
  expect_equal(fusion@gene_upstream@chromosome, "chr11")
})

test_that("import_fusionmap fails to import data", {
  expect_error(import_fusionmap("Non-existent-file", "hg18"))
  expect_error(import_fusionmap(fusionmap_data, "Invalid-genome"))
  expect_error(import_fusionmap(fusionmap_data, "hg19", "invalid limit"))
})

test_that("import_fusionmap limits work", {
  expect_equal(length(import_fusionmap(fusionmap_data, "hg19", 2)), 2)
})
