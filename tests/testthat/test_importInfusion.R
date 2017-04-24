context("Test importing InFusion data into a list of Fusion objects")

infusionData <- system.file(
  "extdata",
  "infusion_fusions.txt",
  package = "chimeraviz")
fusions <- importInfusion(infusionData, "hg19", 3)
fusion <- fusions[[1]]

test_that("importJaffa imports data into Fusion objects", {
  expect_equal(length(fusions), 3)
  expect_equal(fusion@id, "283813")
  expect_equal(fusion@geneA@name, "TMPRSS2")
  expect_equal(fusion@geneA@chromosome, "chr21")
})

test_that("importInfusion fails to import data", {
  expect_error(importInfusion("Non-existent-file", "hg18"))
  expect_error(importInfusion(infusionData, "Invalid-genome"))
  expect_error(importInfusion(infusionData, "hg19", "invalid limit"))
})

test_that("importInfusion limits work", {
  expect_equal(length(importInfusion(infusionData, "hg19", 2)), 2)
})
