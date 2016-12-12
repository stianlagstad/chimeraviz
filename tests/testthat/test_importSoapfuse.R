context("Test importing soapfuse data into a list of Fusion objects")

soapfuse833ke <- system.file(
  "extdata",
  "soapfuse_833ke_final.Fusion.specific.for.genes",
  package = "chimeraviz")
fusions <- importSoapfuse(soapfuse833ke, "hg19", 3)
fusion <- fusions[[1]]

test_that("importSoapfuse imports data into Fusion objects", {
  expect_equal(length(fusions), 3)
  expect_equal(fusion@id, "1")
  expect_equal(fusion@geneA@name, "APLP2")
  expect_equal(fusion@geneA@chromosome, "chr11")
})

test_that("importSoapfuse fails to import data", {
  expect_error(importSoapfuse("Non-existent-file", "hg38"))
  expect_error(importSoapfuse(soapfuse833ke, "Invalid-genome"))
  expect_error(importSoapfuse(soapfuse833ke, "hg38", "invalid limit"))
})

test_that("importSoapfuse limits work", {
  expect_equal(length(importSoapfuse(soapfuse833ke, "hg38", 2)), 2)
})
