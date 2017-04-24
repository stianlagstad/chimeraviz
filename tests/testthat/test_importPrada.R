context("Test importing PRADA data into a list of Fusion objects")

pradaData <- system.file(
  "extdata",
  "PRADA.acc.fusion.fq.TAF.tsv",
  package = "chimeraviz")
fusions <- importPrada(pradaData, "hg19", 3)
fusion <- fusions[[1]]

test_that("importPrada imports data into Fusion objects", {
  expect_equal(length(fusions), 3)
  expect_equal(fusion@id, "1")
  expect_equal(fusion@geneA@name, "TBCD")
  expect_equal(fusion@geneA@chromosome, "chr17")
  expect_equal(fusion@geneA@strand, "+")
})

test_that("importPrada fails to import data", {
  expect_error(importPrada("Non-existent-file", "hg18"))
  expect_error(importPrada(pradaData, "Invalid-genome"))
  expect_error(importPrada(pradaData, "hg19", "invalid limit"))
})

test_that("importPrada limits work", {
  expect_equal(length(importPrada(pradaData, "hg19", 2)), 2)
})
