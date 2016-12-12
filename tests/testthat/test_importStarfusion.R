context("Test importing STAR-Fusion data into a list of Fusion objects")

starfusionData <- system.file(
  "extdata",
  "star-fusion.fusion_candidates.final.abridged.txt",
  package = "chimeraviz")
fusions <- importStarfusion(starfusionData, "hg19", 3)
fusion <- fusions[[1]]

test_that("importStarfusion imports data into Fusion objects", {
  expect_equal(length(fusions), 3)
  expect_equal(fusion@id, "1")
  expect_equal(fusion@geneA@name, "CACNG6")
  expect_equal(fusion@geneA@chromosome, "chr19")
})

test_that("importStarfusion fails to import data", {
  expect_error(importStarfusion("Non-existent-file", "hg18"))
  expect_error(importStarfusion(starfusionData, "Invalid-genome"))
  expect_error(importStarfusion(starfusionData, "hg19", "invalid limit"))
})

test_that("importStarfusion limits work", {
  expect_equal(length(importStarfusion(starfusionData, "hg19", 2)), 2)
})
