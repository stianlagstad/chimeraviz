context("Test importing FusionMap data into a list of Fusion objects")

fusionmapData <- system.file(
  "extdata",
  "FusionMap_01_TestDataset_InputFastq.FusionReport.txt",
  package = "chimeraviz")
fusions <- importFusionmap(fusionmapData, "hg19", 3)
fusion <- fusions[[1]]

test_that("importFusionmap imports data into Fusion objects", {
  expect_equal(length(fusions), 3)
  expect_equal(fusion@id, "FUS_11184689_446432191(++)")
  expect_equal(fusion@geneA@name, "FADS3")
  expect_equal(fusion@geneA@chromosome, "chr11")
})

test_that("importFusionmap fails to import data", {
  expect_error(importFusionmap("Non-existent-file", "hg18"))
  expect_error(importFusionmap(fusionmapData, "Invalid-genome"))
  expect_error(importFusionmap(fusionmapData, "hg19", "invalid limit"))
})

test_that("importFusionmap limits work", {
  expect_equal(length(importFusionmap(fusionmapData, "hg19", 2)), 2)
})
