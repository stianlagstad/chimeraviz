context("Test importing EricScript data into a list of Fusion objects")

ericscriptData <- system.file(
  "extdata",
  "ericscript_SRR1657556.results.total.tsv",
  package = "chimeraviz")
fusions <- importEricscript(ericscriptData, "hg19", 3)
fusion <- fusions[[1]]

test_that("importJaffa imports data into Fusion objects", {
  expect_equal(length(fusions), 3)
  expect_equal(fusion@id, "1")
  expect_equal(fusion@geneA@name, "PPIA")
  expect_equal(fusion@geneA@chromosome, "chr7")
})

test_that("importEricscript fails to import data", {
  expect_error(importEricscript("Non-existent-file", "hg18"))
  expect_error(importEricscript(ericscriptData, "Invalid-genome"))
  expect_error(importEricscript(ericscriptData, "hg19", "invalid limit"))
})

test_that("importEricscript limits work", {
  expect_equal(length(importEricscript(ericscriptData, "hg19", 2)), 2)
})
