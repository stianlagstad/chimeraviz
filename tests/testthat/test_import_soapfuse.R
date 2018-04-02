context("Test importing soapfuse data into a list of Fusion objects")

soapfuse833ke <- system.file(
  "extdata",
  "soapfuse_833ke_final.Fusion.specific.for.genes",
  package = "chimeraviz")
fusions <- import_soapfuse(soapfuse833ke, "hg19", 3)
fusion <- fusions[[1]]

test_that("import_soapfuse imports data into Fusion objects", {
  expect_equal(length(fusions), 3)
  expect_equal(fusion@id, "1")
  expect_equal(fusion@gene_upstream@name, "APLP2")
  expect_equal(fusion@gene_upstream@chromosome, "chr11")
})

test_that("import_soapfuse fails to import data", {
  expect_error(import_soapfuse("Non-existent-file", "hg38"))
  expect_error(import_soapfuse(soapfuse833ke, "Invalid-genome"))
  expect_error(import_soapfuse(soapfuse833ke, "hg38", "invalid limit"))
})

test_that("import_soapfuse limits work", {
  expect_equal(length(import_soapfuse(soapfuse833ke, "hg38", 2)), 2)
})
