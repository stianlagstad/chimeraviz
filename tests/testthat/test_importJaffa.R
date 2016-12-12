context("Test importing JAFFA data into a list of Fusion objects")

jaffaData <- system.file(
  "extdata",
  "jaffa_results.csv",
  package = "chimeraviz")
fusions <- importJaffa(jaffaData, "hg19", 3)
fusion <- fusions[[1]]

test_that("importJaffa imports data into Fusion objects", {
  expect_equal(length(fusions), 3)
  expect_equal(fusion@id, "1")
  expect_equal(fusion@geneA@name, "BCAS4")
  expect_equal(fusion@geneA@chromosome, "chr20")
  expect_equal(fusion@geneA@strand, "*")
})

test_that("importJaffa fails to import data", {
  expect_error(importJaffa("Non-existent-file", "hg18"))
  expect_error(importJaffa(jaffaData, "Invalid-genome"))
  expect_error(importJaffa(jaffaData, "hg19", "invalid limit"))
})

test_that("importJaffa limits work", {
  expect_equal(length(importJaffa(jaffaData, "hg19", 2)), 2)
})
